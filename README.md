# fqdedup
Single-ended FASTQ file deduplication


```plain
usage: fqdedup.py [-h] [-n nt] [--use-mean-qual] [--version] fastq

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

positional arguments:
  fastq            Path to input FASTQ file. Both gzipped and plain FASTQ
                   formats are supported

optional arguments:
  -h, --help       show this help message and exit
  -n nt            Length [nt] of sequence initial portion to search for N.
                   Default: 0.
  --use-mean-qual  Select sequences based on mean quality instead of quality
                   sum.
  --version        show program's version number and exit

```