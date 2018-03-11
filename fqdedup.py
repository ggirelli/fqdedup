#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ==============================================================================
# 
# author: Gabriele Girelli
# email: gigi.ga90@gmail.com
# version: 1.1.1dev
# date: 180308
# project: pre-processing sequencing data
# 
# credits:
#   Dr. F. Agostini for the nice chats and for providing an initial prototype.
# 
# aim:
#   Deduplicate FASTQ file: remove duplicate sequence by keeping the higher
#   quality one. Moreover, remove reads with "N" in the initial portion of a
#   read, if requested by the user.
# 
# description:
#   Initially, the records in the FASTQ are quickly counted with bash "wc -l".
#   Then, the full FASTQ file is read and parsed with Bio.SeqIO. Each record is
#   stored in plain text format alongside its quality in a dictionary, with its
#   sequence as key. If "N" (i.e., any nucleotide) is found in the initial
#   portion (user-defined) of a sequence, the sequence is discarded. Each
#   sequence is compared to the encountered ones and replaces it only and only
#   if its quality is higher (either sum or mean). It is also possible to
#   manually set an upper limit of resident memory using the --max-mem option.
# 
# notes:
#   The current implementation requires less RAM than previous ones, and shorter
#   times to compute. Instead of storing each FASTQ record as parsed, it stores
#   them as plain text alongside sequence and its quality (minor redundancy).
#   For a 20 GB plain FASTQ, approx. 15 GB of resident memory are required.
#	
# ==============================================================================

# DEPENDENCIES =================================================================

import argparse
import binascii
from Bio import SeqIO
import gzip
from io import StringIO
import numpy as np
import os
import psutil
from subprocess import check_output
import sys
from tqdm import tqdm


# PARAMETERS ===================================================================

# Add script description
parser = argparse.ArgumentParser(description = '''
author: Gabriele Girelli
email: gigi.ga90@gmail.com
version: 1.1.1dev
date: 180308
project: pre-processing sequencing data

credits:
  Dr. F. Agostini for the nice chats and for providing an initial prototype.

aim:
  Deduplicate FASTQ file: remove duplicate sequence by keeping the higher
  quality one. Moreover, remove reads with "N" in the initial portion of a
  read, if requested by the user.

description:
  Initially, the records in the FASTQ are quickly counted with bash "wc -l".
  Then, the full FASTQ file is read and parsed with Bio.SeqIO. Each record is
  stored in plain text format alongside its quality in a dictionary, with its
  sequence as key. If "N" (i.e., any nucleotide) is found in the initial
  portion (user-defined) of a sequence, the sequence is discarded. Each
  sequence is compared to the encountered ones and replaces it only and only
  if its quality is higher (either sum or mean). It is also possible to
  manually set an upper limit of resident memory using the --max-mem option.

notes:
  The current implementation requires less RAM than previous ones, and shorter
  times to compute. Instead of storing each FASTQ record as parsed, it stores
  them as plain text alongside sequence and its quality (minor redundancy).
  For a 20 GB plain FASTQ, approx. 15 GB of resident memory are required.
''', formatter_class = argparse.RawDescriptionHelpFormatter)

# Add mandatory arguments
parser.add_argument('fastq', type = str, nargs = 1,
	help = '''Path to input FASTQ file.
	Both gzipped and plain FASTQ formats are supported''')

# Add arguments with default value
parser.add_argument('-n', type = int, nargs = 1, metavar = 'nt', default = [0],
	help = """Length [nt] of sequence initial portion to search for N.
	Default: 0.""")
parser.add_argument('--max-mem', type = int, nargs = 1, metavar = 'MB',
	help = """Upper limit (in MB) of resident memory for the deduplication
	process. Use -1 for unlimited. Not compatible with MacOS. Default: -1.""",
	default = [-1])

# Add flags
parser.add_argument('--use-mean-qual',
	action = 'store_const', dest = 'doMean',
	const = True, default = False,
	help = 'Select sequences based on mean quality instead of quality sum.')

# Version flag
version = "1.1.1dev"
parser.add_argument('--version', action = 'version',
	version = '%s v%s' % (sys.argv[0], version,))

# Parse arguments
args = parser.parse_args()

# Assign to in-script variables
ipath = args.fastq[0]
basename = os.path.splitext(os.path.basename(ipath))[0]
linker_length = args.n[0]

max_mem = args.max_mem[0]
if max_mem < 0:
	max_mem = np.inf

doMean = args.doMean
if doMean:
	qCalc = lambda x: np.mean(x)
else:
	qCalc = lambda x: np.sum(x)

# FUNCTIONS ====================================================================

def get_mem():
	# Memory profiling
	# From https://goo.gl/HkfNpu
	import resource
	rusage_denom = 1024.
	if sys.platform == 'darwin':
		rusage_denom = rusage_denom * rusage_denom
	mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
	return mem

def check_mem():
	# Print memory profiling
	print("%f MB" % (get_mem(),))

def set_default(v, default):
	# Set default value for function argument
	if type(v) == type(None): return(default)
	return(v)

def is_gz_file(filepath):
	# https://stackoverflow.com/a/47080739
	with open(filepath, 'rb') as test_f:
		return binascii.hexlify(test_f.read(2)) == b'1f8b'

def write_output(oh, records):
	'''
	Write output after filtering.

	Args:
		oh (file/gzip): output file handle.
		records (dict): records dictionary after filtering.
	'''

	# Convert to list to empty mem
	records = list(records.values())
	for i in tqdm(range(len(records))):
		# Pop to empty mem sooner
		oh.write(records.pop(0)[0])

def cmp_record(rec, records, ncounter, linker_length):
	'''
	Compares a record to stored one, replace previous ones based on quality and
	discard if "N" is present in the initial portion of the sequence.

	Args:
		rec (SeqRecord): single FASTQ record.
		records (dict): records dictionary during filtering.
		ncounter (int): number of records discarded due to "N".
		linker_length (int): Length [nt] of sequence portion to search for N.

	Returns:
		dict: records dictionary after comparison.
	'''

	# Extract record's sequence, let's make it comfy
	seq = str(rec.seq)

	# Skip if N in linker sequence
	if "N" not in seq[:linker_length]:
		# Prepare record for storage
		ltmp = [qCalc(rec.letter_annotations["phred_quality"]),
			rec.format("fastq")]

		if seq not in records.keys():
			# Store record
			records[seq] = ltmp
		elif ltmp[1] > records[seq][1]:
			# Replace stored record
			records[seq] = ltmp
	else:
		ncounter += 1

	return((records, ncounter))

def log_result(ncounter, nrecs):
	print("%d records removed due to presence of 'N'." % ncounter)
	print("%d records after deduplication." % nrecs)
	print("Peaked at %.1f MB of resident memory." % get_mem())

def run(ih, oh, linker_length, nrecs):
	'''
	Run the script on the input file: remove records with Ns in the initial
	portion of the sequence, if linker_length is larger than 0.

	Args:
		ih (file/gzip): input file handle.
		oh (file/gzip): output file handle.
		linker_length (int): Length [nt] of sequence portion to search for N.
		nrecs (int): expected number of records.
	'''

	# Read all records
	print("Reading and filtering...")
	records = {}
	ncounter = 0

	for rec in tqdm(SeqIO.parse(ih, "fastq"), total = nrecs):
		records, ncounter = cmp_record(rec, records, ncounter, linker_length)
	log_result(ncounter, len(records))

	# Remove duplicates and write output
	print("Writing...")
	write_output(oh, records)

def run_mm(ih, oh, linker_length, nrecs, max_mem = None):
	'''
	Run the script on the input file: remove records with Ns in the initial
	portion of the sequence, if linker_length is larger than 0.
	Performs resident memory profiling with upper limit set by the user.

	Args:
		ih (file/gzip): input file handle.
		oh (file/gzip): output file handle.
		linker_length (int): Length [nt] of sequence portion to search for N.
		nrecs (int): expected number of records.
		max_mem (int): upper resident memory limit in MB.
	'''

	# Default memory limit to infinity
	if max_mem < 0: max_mem = None
	max_mem = set_default(max_mem, np.inf)

	# Read all records
	print("Reading and filtering...")
	records = {}
	ncounter = 0

	for rec in tqdm(SeqIO.parse(ih, "fastq"), total = nrecs):
		# Stop when the mem limit is hit
		if get_mem() >= max_mem:
			sys.exit("!ABORTED! Hit resident memory limit of %d MB." % (
				max_mem,))

		# Compare current record with stored ones
		records, ncounter = cmp_record(rec, records, ncounter, linker_length)
	log_result(ncounter, len(records))

	# Remove duplicates and write output
	print("Writing...")
	write_output(oh, records)

# RUN ==========================================================================

# Log input --------------------------------------------------------------------

print("\n# fqdedup.py v%s - Single-end FASTQ deduplication" % version)
print("Input: %s" % (ipath,))

if ( is_gz_file(ipath) ):
	print("! Gzipped FASTQ deduplication style.")
	# Prepare to parse a gzipped FASTQ input
	catter = "zcat"
	opath = "%s/%s.dedup.gz" % (os.path.dirname(ipath), basename)
	oh = gzip.open(opath, "wt")
	ih = gzip.open(ipath, "rt")
else:
	print("! Plain FASTQ deduplication style.")
	# Prepare to parse a plain FASTQ input
	catter = "cat"
	opath = "%s/%s.dedup.fastq" % (os.path.dirname(ipath), basename)
	oh = open(opath, "wt")
	ih = open(ipath, "rt")

if doMean: print("!Using average quality for sequence selection.")
else: print("! Using quality sum for sequence selection.")

if 0 != linker_length:
	print("! Discarding sequences with N in the first %d bases." % (
		linker_length,))

if np.inf != max_mem:
	print("! Upper resident memory limit set to %d MB." % (max_mem,))

print()

# Count records in input -------------------------------------------------------

print("Counting records...")
nrecs = float(check_output(
	["bash", "-c", "%s '%s' | wc -l" % (catter, ipath)])) / 4
print("> Found %d records." % (nrecs,))

# Run --------------------------------------------------------------------------

if np.inf == max_mem:
	run(ih, oh, linker_length, nrecs)
else:
	run_mm(ih, oh, linker_length, nrecs, max_mem)

# END ==========================================================================

################################################################################
