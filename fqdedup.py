#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ==============================================================================
# 
# 180308 - Gabriele Girelli
# Project: pre-processing sequencing data
# 
# Credits:
# 	Dr. F. Agostini for suggestions on the script structure.
# 
# Aim: deduplicate FASTQ file.
# 	- Remove reads with N in a certain window.
# 	- Remove duplicate sequence by keeping the higher quality one.
# 	
# Notes:
# 	Current implementation can be used only if the whole fastq file can be
# 	stored in the machine's RAM.
#	
# ==============================================================================

# DEPENDENCIES =================================================================

import gzip
import os
from subprocess import check_output
import sys
from tqdm import tqdm

from Bio import SeqIO

# PARAMETERS ===================================================================

zip_path = "/media/gire/MiSo/TK93_fastq_dedup/TK93_S5_LALL_R1_001.sample.fastq.gz"
linker_length = 22

# FUNCTIONS ====================================================================


# RUN ==========================================================================

# Log input
print("Input: %s" % (zip_path,))

# Count records in input
print("Counting records...")
nrecs = check_output(["bash", "-c", "zcat '%s' | wc -l" % zip_path])
nrecs = float(nrecs) / 4
print("> Found %d records." % (nrecs,))

# Prepare empty dictionary
dseq = {}

# Open output handle
basename = os.path.splitext(os.path.basename(zip_path))[0]
opath = "%s/%s.dedup.gz" % (os.path.dirname(zip_path), basename)
oh = gzip.open(opath, "wt")

# Read zipped fastq
with gzip.open(zip_path, "rt") as ih:

	# Read all records
	print("Reading input...")
	records = []
	for record in tqdm(SeqIO.parse(ih, "fastq"), total = nrecs):
		seq = str(record.seq)

		# Skip if N in linker sequence
		if "N" in seq[:linker_length]:
			continue
		else:
			records.append(record)
	print("%d records without N in the linker region." % (len(records),))

	# Sort the reads by quality (higher on top)
	print("Sorting by quality...")
	records.sort(key = lambda x: -sum(x.letter_annotations['phred_quality']))

	# Remove duplicates and write output
	print("Deduplicating...")
	encountered = set()
	for i in tqdm(range(len(records))):
		if str(records[i].seq) not in encountered:
			encountered.add(str(records[i].seq))
			oh.write(records[i].format("fastq"))
	print("%d records after deduplication." % len(encountered))

# Close output handle
oh.close()

# END ==========================================================================

################################################################################
