require(GenomicRanges)
require(rtracklayer)

refFlat <- read.table("repeats.refFlat")


## Prepare genes info
repeats_genes_gtf = GRanges(refFlat$V3,
                            IRanges(refFlat$V5,refFlat$V6),
                            type=rep("gene",nrow(refFlat)),                    
                            strand=refFlat$V4, 
                            
                            # gene info
                            gene_id=refFlat$V1, 
                            gene_name=refFlat$V1,
                            
                            # transcript info
                            transcript_id=rep(NA,nrow(refFlat)),
                            transcript_name= rep(NA, nrow(refFlat)),

                            # exon info
                            exon_number = rep(NA,nrow(refFlat)),
                            exon_id = rep(NA, nrow(refFlat))
                            )


## Pepare transcripts info
repeats_transcript_gtf = GRanges(refFlat$V3,
                            IRanges(refFlat$V5,refFlat$V6),
                            type=rep("transcript",nrow(refFlat)),                    
                            strand=refFlat$V4, 
                            
                            # gene info
                            gene_id=refFlat$V1, 
                            gene_name=refFlat$V1,
                            
                            # transcript info
                            transcript_id=paste0(refFlat$V1,"-001"),
                            transcript_name= paste0(refFlat$V1,"-001"),

                            # exon info
                            exon_number = rep(NA,nrow(refFlat)),
                            exon_id = rep(NA, nrow(refFlat))
                            )



## Pepare exon info
repeats_exons_gtf = GRanges(refFlat$V3,
                            IRanges(refFlat$V5,refFlat$V6),
                            type=rep("exon",nrow(refFlat)),                    
                            strand=refFlat$V4, 
                            
                            # gene info
                            gene_id=refFlat$V1, 
                            gene_name=refFlat$V1,
                            
                            # transcript info
                            transcript_id=paste0(refFlat$V1,"-001"),
                            transcript_name= paste0(refFlat$V1,"-001"),

                            # exon info
                            exon_number = rep(1,nrow(refFlat)),
                            exon_id = paste0("Exon-",refFlat$V1)
                            )


## Add them in order 

repeats_gtf <- GRanges()

for(i in 1:nrow(refFlat)){
    tmp <- c(repeats_genes_gtf[i], repeats_transcript_gtf[i], repeats_exons_gtf[i])
    repeats_gtf <- c(repeats_gtf,tmp)
}





