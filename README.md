# fqdedup
Single-ended FASTQ file deduplication


```plain
usage: fqdedup.py [-h] [-n nt] [--max-mem MB] [--use-mean-qual] [--version]
                  fastq

author: Gabriele Girelli
email: gigi.ga90@gmail.com
version: 1.1.0
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

positional arguments:
  fastq            Path to input FASTQ file. Both gzipped and plain FASTQ
                   formats are supported

optional arguments:
  -h, --help       show this help message and exit
  -n nt            Length [nt] of sequence initial portion to search for N.
                   Default: 0.
  --max-mem MB     Upper limit (in MB) of resident memory for the
                   deduplication process. Use -1 for unlimited. Not compatible
                   with MacOS. Default: -1.
  --use-mean-qual  Select sequences based on mean quality instead of quality
                   sum.
  --version        show program's version number and exit
```