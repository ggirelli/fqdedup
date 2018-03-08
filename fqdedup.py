#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ==============================================================================
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Date: 180308
# Project: pre-processing sequencing data
# 
# Credits:
# 	Dr. F. Agostini for suggestions on the script structure and for providing
# 		an initial prototype.
# 
# Aim:
# 	Deduplicate FASTQ file: remove duplicate sequence by keeping the higher
# 	quality one. Moreover, remove reads with "N" in the initial portion of a
# 	read, if requested by the user.
# 
# Description:
# 	The full FASTQ file is read and parsed with Bio.SeqIO. The records are
# 	stored in plain text format in a list, and in a separate array their
# 	qualities. During the read/parse operation, if "N" (i.e., any nucleotide)
# 	is found in the initial portion (user-defined) of a sequence, the sequence
# 	is discarded. Sorted indexes from the quality array are used to navigate
# 	through the records list and keep only the first occurrence (highest
# 	quality) of each sequence.
# 
# Notes:
# 	For a 20 GB plain FASTQ, approximately 30 GB of free RAM are required.
# 	Designed to work on single-end FASTQ files.
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
from subprocess import check_output
import sys
from tqdm import tqdm

# PARAMETERS ===================================================================

# Add script description
parser = argparse.ArgumentParser(description = '''
author: Gabriele Girelli
email: gigi.ga90@gmail.com
version: 1.0.0
date: 180308
project: pre-processing sequencing data

credits:
  Dr. F. Agostini for suggestions on the script structure and for providing
  an initial prototype.

aim:
  Deduplicate FASTQ file: remove duplicate sequence by keeping the higher
  quality one. Moreover, remove reads with "N" in the initial portion of a
  read, if requested by the user.

description:
  The full FASTQ file is read and parsed with Bio.SeqIO. The records are
  stored in plain text format in a list, and in a separate array their
  qualities. During the read/parse operation, if "N" (i.e., any nucleotide)
  is found in the initial portion (user-defined) of a sequence, the sequence
  is discarded. Sorted indexes from the quality array are used to navigate
  through the records list and keep only the first occurrence (highest
  quality sum) of each sequence. Use the --use-mean-qual option to select the
  records with the highest quality mean instead.

notes:
  The current implementation requires less RAM than previous ones, but takes
  longer times to compute. Instead of storing each FASTQ record as parsed,
  it stores them separately as strings and then parses them again when needed.
  For a 20 GB plain FASTQ, approximatley 30 GB of free RAM are required. Also,
  in this version, quality is calculated in tandem with the input, and sorting
  is performed on the quality list, and matched to the input records.
''', formatter_class = argparse.RawDescriptionHelpFormatter)

# Add mandatory arguments
parser.add_argument('fastq', type = str, nargs = 1,
	help = '''Path to input FASTQ file.
	Both gzipped and plain FASTQ formats are supported''')

# Add arguments with default value
parser.add_argument('-n', type = int, nargs = 1, metavar = 'nt', default = [0],
	help = """Length [nt] of sequence initial portion to search for N.
	Default: 0.""")

# Add flags
parser.add_argument('--use-mean-qual',
	action = 'store_const', dest = 'doMean',
	const = True, default = False,
	help = 'Select sequences based on mean quality instead of quality sum.')

# Version flag
version = "1.0.0"
parser.add_argument('--version', action = 'version',
	version = '%s v%s' % (sys.argv[0], version,))

# Parse arguments
args = parser.parse_args()

# Assign to in-script variables
ipath = args.fastq[0]
basename = os.path.splitext(os.path.basename(ipath))[0]
linker_length = args.n[0]
doMean = args.doMean

if doMean:
	qCalc = lambda x: np.mean(x)
else:
	qCalc = lambda x: np.sum(x)

# FUNCTIONS ====================================================================


def is_gz_file(filepath):
	# https://stackoverflow.com/a/47080739
	with open(filepath, 'rb') as test_f:
		return binascii.hexlify(test_f.read(2)) == b'1f8b'

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
	print("Reading input...")
	qs = []
	records = []

	for record in tqdm(SeqIO.parse(ih, "fastq"), total = nrecs):
		seq = str(record.seq)
		# Skip if N in linker sequence
		if "N" in seq[:linker_length]:
			continue
		else:
			records.append(record.format("fastq"))
			qs.append(qCalc(record.letter_annotations['phred_quality']))

	print("%d records without N in the linker region." % (len(records),))

	# Sort the reads by quality (higher on top)
	print("Sorting by quality...")
	idxs = np.argsort(qs).tolist()[::-1]

	# Remove duplicates and write output
	print("Deduplicating...")
	encountered = set()

	for i in tqdm(idxs):
		record = SeqIO.read(StringIO(records[i]), "fastq")

		if str(record.seq) not in encountered:
			encountered.add(str(record.seq))
			oh.write(record.format("fastq"))

	print("%d records after deduplication." % len(encountered))

# RUN ==========================================================================

# Log input
print("\n# fqdedup3.py - Single-end FASTQ deduplication")
print("Input: %s" % (ipath,))
if doMean: print("!Using average quality for sequence selection.")
else: print("! Using quality sum for sequence selection.")

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

print()

# Count records in input
print("Counting records...")
nrecs = float(check_output(
	["bash", "-c", "%s '%s' | wc -l" % (catter, ipath)])) / 4
print("> Found %d records." % (nrecs,))

# Run
run(ih, oh, linker_length, nrecs)

# END ==========================================================================

################################################################################
