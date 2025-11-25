# CorrectR1
Pairwise correction of amplicon MiSeq paired reads

CorrectR1_v1.12.py\
Parallel R1 correction across multiple samples in a directory.

Inputs organised as:\
  input_directory/ \
      sample1_correctedR1_R1.fastq.gz\
      sample1_correctedR1_R2.fastq.gz\
      sample2_correctedR1_R1.fastq.gz\
      sample2_correctedR1_R2.fastq.gz\
      ...\
      (fastq file names must end in R1 or R2, not 001)\
Outputs organised as:\
  R1corrected_<params>_reads/ \
      sample1_correctedR1_R1.fastq.gz\
      sample1_correctedR1_R2.fastq.gz\
      ...\
  R1corrected_<params>_figures/ \
      sample1_correction_counts.png \
      sample1_correction_positions.png\
      sample1_overlap_heatmap.png\
      ...\
  combined_summary.tsv\
  combined_run.log\

Conducts pairwise alignment of R1 and R2. Where there is a mismatch in the alignment, the base with the higher quality score is written to R1. Indels are not currently considered.\
Tested on HTT Exon1 repeat amplicon MiSeq data only.\
Parallelises samples in directory, cannot currently specify a single sample be processed - create a directory with R1 and R2 from this sample.\

Recommend using default alignment parameters (i.e. do not set m, mm, go, ge).\
Recommend trimming R2 reverse complement from -150 to -30. For example, a 200 base R2 (non-reverse complement) will have 30 bases trimmed from the start and 50 bases trimmed from the end; a 300 base R2 will have 30 bases trimmed from the start and 150 bases trimmed from the end.\
--trimstart 150 --trimend 30\
\

usage: CorrectR1_v1.12.py [-h] -d DIR [-t T] [-s S] [--trimstart TRIMSTART] [--trimend TRIMEND] [-m M] [-mm MM] [-go GO] [-ge GE] [--motif MOTIF]

Batch R1 correction using R2 overlap

options:\
  -h, --help            show this help message and exit\
  -d, --dir DIR         Directory containing paired FASTQ files\
  -t T                  Number of threads (samples in parallel)\
  -s S                  Reads to skip from start\
  --trimstart TRIMSTART\
  --trimend TRIMEND\
  -m M\
  -mm MM\
  -go GO\
  -ge GE\
  --motif MOTIF\

Example:
python correct_R1_v1.12.py -d "/path/to/directory" -t 8 --trimstart 150 --trimend 30
