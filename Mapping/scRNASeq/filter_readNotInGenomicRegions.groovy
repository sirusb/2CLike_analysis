
get_nonGenomicReads = {
  
 produce("IntergenicReads.bed"){
 exec """
   python extract_notMappedToGenes.py ${input} ${output}
 """
 }
}


extractReads  = {

  def read1 = "/nfs2/nadhir/Projects/201804_xiaoji_Dux/fastq/0H/read1.fastq.gz"
  def read2 = "/nfs2/nadhir/Projects/201804_xiaoji_Dux/fastq/0H/read2.fastq.gz"
  exec """
    seqtk subseq ${read1} ${input} > read1.fastq &&
    seqtk subseq ${read2} ${input} > read2.fastq
  """
}



run {

  "%.bam" * [get_nonGenomicReads + extractReads]
}

