about title: "Drop-Seq pipeline"

Picard='java -Xmx4g -jar /nfs2/nadhir/tools/Drop-seq_tools-1.13/3rdParty/picard/picard.jar'

FastqToSam = {
	 doc title: "Convert Fastq to a samfile",      
      author: "Nadhir (djek.nad@gmail.com)"

	produce("starting.bam"){
		exec """
			$Picard FastqToSam FASTQ=${input1} FASTQ2=${input2} OUTPUT=$output SAMPLE_NAME=Falong_neg_400cell
		"""
	}
	
}

DropSeqTools="/nfs2/nadhir/tools/Drop-seq_tools-1.13/"
ExtractCellBarCode = {
	produce("${input.prefix}_cell.bam"){
		exec """
		$DropSeqTools/TagBamWithReadSequenceExtended   
				SUMMARY=./unaligned_tagged_Cellular.bam_summary.txt 
				BASE_RANGE=1-12 BASE_QUALITY=10 
				BARCODED_READ=1 DISCARD_READ=false 
				TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 
				INPUT=$input
				OUTPUT=${output}
	"""
	}	
}


ExtractMolecularBarCode = {
	produce("${input.prefix}_molecule.bam"){
		exec """
			$DropSeqTools/TagBamWithReadSequenceExtended SUMMARY=./unaligned_tagged_Cellular_molecular.bam_summary.txt 
				BASE_RANGE=13-20 BASE_QUALITY=10 
				BARCODED_READ=1 DISCARD_READ=true 
				TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 
				INPUT=$input 
				OUTPUT=${output}
		"""
	}
}

FilterBam = {
	produce("${input.prefix}_filtered_unaligned.bam"){
		exec """
		$DropSeqTools/FilterBAM TAG_REJECT=XQ 
				INPUT=$input
				OUTPUT=${output}		
	"""	
	}	
}


TrimAdaptor = {
	def Sequence="AAGCAGTGGTATCAACGCAGAGTGAATGGG"
	produce("${input.prefix}_smart.bam"){
		exec """
			$DropSeqTools/TrimStartingSequence INPUT=${input}
				OUTPUT=${output}
				OUTPUT_SUMMARY=${input_prefix}_adapter_trimming_report.txt
				SEQUENCE=${Sequence}
				MISMATCHES=0 
				NUM_BASES=5
		"""
	}
}

TrimPolyA = {
	produce("${input.prefix}_polyA.bam"){
		exec """
			$DropSeqTools/PolyATrimmer INPUT=${input}
			OUTPUT=${output}
			OUTPUT_SUMMARY=polyA_trimming_report.txt 
			MISMATCHES=0 NUM_BASES=6
		"""
	}
}

SamToFastq = {
	transform(".bam") to (".fastq"){
		exec """
			$Picard SamToFastq INPUT=${input}
			FASTQ=${output}
		"""
	}
}


GENOME="/nfs2/nadhir/genomes/mm9/RepEnrich/rms_Pseudo_out/bowtie_index/rmsk_pseudo"
bowtie1Map = {

  if ( inputs.size() == 1){
        produce("BowtieAligned.out.sam"){
        exec """
                bowtie --best --strata --sam -m 2 -p 10 $GENOME $input $output
        """
        }
  }
  else{
      produce("BowtieAligned.out.sam"){
        exec """
         bowtie --sam -m 2 --best --strata -p10 $GENOME -1 $input1 -2 $input2 $output
       """
      }
  }
}

STAR_align = {
	def GenomeDir="/nfs2/nadhir/genomes/mm9/RepEnrich/rms_Pseudo_out/STAR_index"
	produce("starAligned.out.sam"){
        exec """
		STAR --genomeDir ${GenomeDir} --readFilesIn ${input} 
			--outFileNamePrefix star --outReadsUnmapped Fastq 
			--outFilterIntronMotifs RemoveNoncanonical
         		--alignIntronMax 1 --runThreadN 10 --outSAMprimaryFlag OneBestScore	

	"""
        }	
}


getPrimary = {
   filter("primary"){
     exec """
        samtools view -F 256 $input -o $output 
     """
   }
}

samtoolsSortSam = {
   produce("${input.prefix}.sorted.bam"){
      exec """
         samtools sort --threads 8 --output-fmt BAM -n $input -o $output 
      """
   }
}

filterUnmmaped = {
  filter("filtered"){
    exec """
      edwBamFilter $input $output  -unmapped
    """
  }
}


SortSam = {
	produce("${input.prefix}.sorted.bam"){
		exec """
			$Picard SortSam I=${input} O=${output} SO=queryname
		"""
	}
}

samtools_MergeAlignedNonAligned = {
   produce("merged.bam"){
     //from("smart_polyA.bam",".primary.bam"){
     from("smart_polyA.bam","sorted.bam"){
	     exec """
      		samtools merge   
	     """
     }
   }
}

MergeAlignedNonAligned = {
	def Refrence="/nfs2/nadhir/genomes/mm9/fasta/mm9_rmsk_pseudo.fa"
	produce("merged.bam"){
		from("smart_polyA.bam",".filtered.sorted.bam"){
			exec """
				$Picard MergeBamAlignment REFERENCE_SEQUENCE=${Refrence} 
				UNMAPPED_BAM=${input1}
				ALIGNED_BAM=${input2}
				OUTPUT=${output} 
				INCLUDE_SECONDARY_ALIGNMENTS=false
			"""
		}	
	}
}




MergeAlignedNonAligned_STAR = {
        def Refrence="/nfs2/nadhir/genomes/mm9/RepEnrich/rms_Pseudo_out/rsm_pseudo.fa"
        produce("merged.bam"){
                //from("smart_polyA.bam",".primary.bam"){
                from("smart_polyA.bam",".sorted.bam"){
                        exec """
                                $Picard MergeBamAlignment REFERENCE_SEQUENCE=${Refrence}
                                UNMAPPED_BAM=${input1}
                                ALIGNED_BAM=${input2}
                                OUTPUT=${output}
                                INCLUDE_SECONDARY_ALIGNMENTS=false
                        """
                }
        }
}


TagWithExons = {
	def AnnoFile="/nfs2/nadhir/genomes/mm9/RepEnrich/rms_Pseudo_out/repeats.refFlat"
	produce("${input.prefix}_exon_tagged.bam"){
		exec """
			$DropSeqTools/TagReadWithGeneExon 
					I=${input} 
					O=${output}
					ANNOTATIONS_FILE=${AnnoFile}
					TAG=GE
		"""
	}
}


DetectBeadSynthErros = {
	def PrimerSeq="AAGCAGTGGTATCAACGCAGAGTAC"
	produce("${input.prefix}_clean.bam"){
	exec """
		$DropSeqTools/DetectBeadSynthesisErrors I=${input}
			O=${output}
			OUTPUT_STATS=my.synthesis_stats.txt 
			SUMMARY=my.synthesis_stats.summary.txt 
			NUM_BARCODES=2000 
			PRIMER_SEQUENCE=${PrimerSeq}
	"""
       }
}

GenerateExprTable = {
        var nbCoreBarCodes : 2000
	produce("${input.prefix}.dge.txt.gz"){
		exec """
			$DropSeqTools/DigitalExpression I=${input} 
				O=${output}
				SUMMARY=out_gene_exon_tagged.dge.summary.txt 
				NUM_CORE_BARCODES=$nbCoreBarCodes
		"""
	}
}

CellNumberMap = {
	produce("cell_number_map.pdf")
	R{"""
		data=read.table("$input",header=F,stringsAsFactors=F)
		x=cumsum(data$V1)
		x=x/max(x)
		pdf("$output")
		plot(1:length(x), x, type='l', 
			 col="blue", xlab="cell barcodes sorted by number of reads [descending]",
			 ylab="cumulative fraction of reads", xlim=c(1,500));
		dev.off()
	"""}
}

BamTagHistogram = {
	produce("out_cell_readcounts.txt.gz"){
		from("_clean.bam"){
			exec """
			$DropSeqTools/BAMTagHistogram I=${input} O=${output} TAG=XC
		"""
		}		
	}
}


run {
	
	   //Using Bowtie for mapping
           //"*.fastq.gz" * [ FastqToSam + ExtractCellBarCode + ExtractMolecularBarCode +
	   //				  FilterBam + TrimAdaptor + TrimPolyA + SamToFastq + 
	   //				  bowtie1Map + SortSam + filterUnmmaped + SortSam  +  MergeAlignedNonAligned + TagWithExons +
	   //				  DetectBeadSynthErros + GenerateExprTable + 
	   //	    			  CellNumberMap + BamTagHistogram ]

           //Using STAR for mapping
           "*.fastq.gz" * [ FastqToSam + ExtractCellBarCode + ExtractMolecularBarCode +
                                          FilterBam + TrimAdaptor + TrimPolyA + SamToFastq +
                                          //STAR_align  + SortSam +  getPrimary + MergeAlignedNonAligned_STAR + TagWithExons +
                                          STAR_align  + SortSam +  MergeAlignedNonAligned_STAR + TagWithExons +
                                          DetectBeadSynthErros + GenerateExprTable +
                                          BamTagHistogram ]
}
	

