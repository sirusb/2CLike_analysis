
bpipe run ../../scripts/DropSeq_workflow.groovy ../../fastq/0H/read1.fastq.gz ../../fastq/0H/read2.fastq.gz
mkdir Repeats
cd Repeats

original_r1="/nfs2/nadhir/Projects/201804_xiaoji_Dux/fastq/0H/read1.fastq.gz"
original_r2="/nfs2/nadhir/Projects/201804_xiaoji_Dux/fastq/0H/read2.fastq.gz"

bpipe run ../../../scripts/filter_readNotInGenomicRegions.groovy -p read1=${original_r1}  -p read2=${original_r2} ../merged_exon_tagged.bam
bpipe run ../../../scripts/DropSeq_rpeats_workflow.groovy read1.fastq.gz read2.fastq.gz
