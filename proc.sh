#!/bin/bash

# Sample slurm submission script for the Chimera compute cluster
# Lines beginning with # are comments, and will be ignored by
# the interpreter.Lines beginning with #SBATCH are directives
# to the scheduler.These in turn can be commented out by
# adding a second # (e.g. ##SBATCH lines will not be processed
# by the scheduler).
#
#
# set name of job
#SBATCH --job-name=verru2
#
# set the number of processors/tasks needed
#SBATCH -n 24

#set an account to use
#if not used then default will be used
# for scavenger users, use this format:
#BATCH --account=patrick.kearns
# for contributing users, use this format:
##SBATCH --account=

# set max wallclock timeDD-HH:MM:SS

# the default time will be 1 hour if not set
#SBATCH --time=00-96:00:00

# set a memory request
#SBATCH --mem=64gb

# Set filenames for stdout and stderr.%j can be used for the jobid.
# see "filename patterns" section of the sbatch man page for
# additional options
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out
#

# set the partition where the job will run.Multiple partitions can
# be specified as a comma separated list
# Use command "sinfo" to get the list of partitions
#SBATCH --partition=DGXA100
##SBATCH --partition=Intel6240,Intel6248,DGXA100

#When submitting to the GPU node, these following three lines are needed:

##SBATCH --gres=gpu:1
##SBATCH --export=NONE
#source /etc/profile
 

#Optional
# mail alert at start, end and/or failure of execution
# see the sbatch man page for other options
#SBATCH --mail-type=ALL
# send mail to this address
#SBATCH --mail-user=patrick.kearns@umb.edu

# Put your job commands here, including loading any needed
# modules or diagnostic echos.


conda activate qiime2-2023.5
source activate qiime2-2023.5

#change to directory where sequences are in
cd  /hpcstor6/scratch01/p/patrick.kearns/6-11-25_run

#load raw FASTQ reads into QIIME
#qiime tools import --type EMPSingleEndSequences --input-path ./data --output-path abby_seqs.qza

#demultiplex reads
qiime demux emp-paired\
  --i-seqs seqs.qza \
 --m-barcodes-file abby_map2.txt \
 --m-barcodes-column BarcodeSequence \
  --o-per-sample-sequences abby_demux2.qza \
  --o-error-correction-details  abby_demux-details2.qza \
  --p-no-golay-error-correction 

qiime demux emp-paired \
  --i-seqs seqs.qza \
 --m-barcodes-file abby_map.txt \
 --m-barcodes-column BarcodeSequence \
  --o-per-sample-sequences abby_demux.qza \
  --o-error-correction-details  abby_demux-details.qza \
  --p-no-golay-error-correction 
  
#quality filer
qiime quality-filter q-score \
--i-demux  abby_demux.qza \
--o-filtered-sequences  abby_demux-filtered.qza \
--o-filter-stats  abby_demux-filter-stats.qza \
--p-min-quality 15

qiime quality-filter q-score \
--i-demux  abby_demux2.qza \
--o-filtered-sequences  abby_demux-filtered2.qza \
--o-filter-stats  abby_demux-filter-stats2.qza \
--p-min-quality 15
 
 
  #call ASVs with deblur
  qiime deblur denoise-16S \
  --i-demultiplexed-seqs  abby_demux-filtered.qza \
  --p-trim-length 120 \
  --p-jobs-to-start 24 \
  --o-representative-sequences  abby_rep-seqs-deblur.qza \
  --o-table  abby_table-deblur.qza \
   --o-stats  abby_deblur-stats.qza

  qiime deblur denoise-16S \
  --i-demultiplexed-seqs  abby_demux-filtered2.qza \
  --p-trim-length 120 \
  --p-jobs-to-start 24 \
  --o-representative-sequences  abby_rep-seqs-deblur2.qza \
  --o-table  abby_table-deblur2.qza \
   --o-stats  abby_deblur-stats2.qza
  
 #export rep seqs
 qiime tools export --input-path abby_rep-seqs-deblur.qza --output-path rep_seqs
 qiime tools export --input-path abby_rep-seqs-deblur2.qza --output-path rep_seqs2
 
 #export OTU table to biom then to text file
 qiime tools export --input-path abby_table-deblur.qza --output-path asv_table
 biom convert -i asv_table/feature-table.biom --to-tsv -o asv_table.txt
 
qiime tools export --input-path abby_table-deblur2.qza --output-path asv_table2
 biom convert -i asv_table/feature-table2.biom --to-tsv -o asv_table2.txt
 

 #assign taxonomy with sklearn and silva database
qiime feature-classifier classify-sklearn   --i-classifier silva-138-99-515-806-nb-classifier.qza   --i-reads  abby_rep-seqs-deblur.qza   --o-classification abby_taxonomy.qza 
qiime feature-classifier classify-sklearn   --i-classifier silva-138-99-515-806-nb-classifier.qza   --i-reads  abby_rep-seqs-deblur2.qza   --o-classification abby_taxonomy2.qza 


 #export taxonomy file
 qiime tools export --input-path abby_taxonomy.qzv --output-path taxonomy
 qiime tools export --input-path abby_taxonomy2.qzv --output-path taxonomy2
 
#make stacked bar visualizations
qiime taxa barplot --o-visualization taxa_plot  --m-metadata-file abby_map.txt  --i-taxonomy abby_taxonomy.qza --i-table abby_table-deblur.qza
qiime taxa barplot --o-visualization taxa_plot2  --m-metadata-file abby_map2.txt  --i-taxonomy abby_taxonomy2.qza --i-table abby_table-deblur2.qza

#normalize table by 16S copy # using gcn-norm
qiime gcn-norm copy-num-normalize \
  --i-table abby_table-deblur2.qza \
  --i-taxonomy abby_taxonomy2.qza \
  --o-gcn-norm-table abby_table-normalized2.qza

qiime gcn-norm copy-num-normalize \
  --i-table abby_table-deblur.qza \
  --i-taxonomy abby_taxonomy.qza \
  --o-gcn-norm-table abby_table-normalized.qza

 #export normalized OTU table to biom then to text file
 qiime tools export --input-path abby_table-normalized2.qza --output-path asv_table2_norm
 biom convert -i asv_table2_norm/feature-table.biom --to-tsv -o asv_table2_norm.txt
 
qiime tools export --input-path abby_table-normalized.qza --output-path asv_table_norm
 biom convert -i asv_table_norm/feature-table.biom --to-tsv -o asv_table_norm.txt
 

 #import plasmid sequence to QIIME
#qiime tools import --input-path DF_16S_complete\ ITS\ construct.fasta --output-path plasmid.qza --type 'FeatureData[Sequence]'

#run deblur
  qiime deblur denoise-other \
  --i-demultiplexed-seqs  abby_demux-filtered2.qza \
  --i-reference-seqs plasmid.qza \
  --p-trim-length 120 \
  --o-representative-sequences  plasmid_rep-seqs-deblur2.qza \
  --o-table  plasmid_table-deblur2.qza \
   --o-stats  plasmid_deblur-stats2.qza

  qiime deblur denoise-other \
  --i-demultiplexed-seqs  abby_demux-filtered2.qza \
  --i-reference-seqs plasmid.qza \
  --p-trim-length 120 \
  --o-representative-sequences  plasmid_rep-seqs-deblur.qza \
  --o-table  plasmid_table-deblur.qza \
   --o-stats  plasmid_deblur-stats.qza

#export tables
 qiime tools export --input-path plasmid_table-deblur.qza --output-path plasmid_table
 biom convert -i plasmid_table/feature-table.biom --to-tsv -o plasmid_table.txt
 
qiime tools export --input-path plasmid_table-deblur2.qza --output-path plasmid_table2
 biom convert -i plasmid_table2/feature-table.biom --to-tsv -o plasmid_table2.txt
 


