title: "RNA-Seq analysis pipeline "

// The mm10 genome
GENOME="/home/nadhir/genomes/mm10/STAR_SynDux_tdTomato_INDEX"
genome_refFlat="/home/nadhir/genomes/mm10/STAR_SynDux_tdTomato_INDEX/GRCm38.85_SynDux_TdTomato.refFlat"
gtf="/home/nadhir/genomes/mm10/STAR_SynDux_tdTomato_INDEX/Mus_musculus.GRCm38.85_fixed_SynDux_TdTomato.gtf"

// Psuedo-genome index
REPEATS_GENOME="/n/groups/zhanglab/nadhir/Pseudogenome_STARINDEX/"
REPEATS_GTF="/n/groups/zhanglab/nadhir/Pseudogenome_STARINDEX/repeats.gtf"


// Location of the Drop-seq tools
DropSeqTools="/home/nadhir/tools/Drop-seq_tools-1.13"

// Location of Trimmomatic
Trimmomatic_root="/n/app/trimmomatic/0.36/bin"
Trimmomatic_LOC="$Trimmomatic_root/trimmomatic-0.36.jar"



// Directory that contains the RNA-Seq data
RootDir="/n/groups/zhanglab/nadhir/Dux_Bulk/fastq"


// Define you samples here
// Format:
// SampleName : ["fastq_location"]

def branches = [
    //D0
    D0WT_Rep1: ["${RootDir}/D0_WT/wt0_1_combined.fq.gz"],
    D0WT_Rep2: ["${RootDir}/D0_WT/wt0_2_combined.fq.gz"],
    
    //D1 Pos
    D1_2CPos_Rep1: ["${RootDir}/D1_2C_Pos/2c1_combined.fq.gz"],
    D1_2CPos_Rep2: ["${RootDir}/D1_2C_Pos/2c2_combined.fq.gz"],
    
    //D1 Neg
    D1_2CNeg_Rep1: ["${RootDir}/D1_2C_Neg/wt1_combined.fq.gz"],
    D1_2CNeg_Rep2: ["${RootDir}/D1_2C_Neg/wt2_combined.fq.gz"],

    //D2 Pos
    D2_2CPos_Rep1: ["${RootDir}/D2_2C_Pos/2c2c1_combined.fq.gz"],
    D2_2CPos_Rep2: ["${RootDir}/D2_2C_Pos/2c2c2_combined.fq.gz"],
 
    //D2 Neg
    D2_2CNeg_Rep1: ["${RootDir}/D2_2C_Neg/2cwt1_combined.fq.gz"],
    D2_2CNeg_Rep2: ["${RootDir}/D2_2C_Neg/2cwt2_combined.fq.gz"]
]




fastqc ={
	doc title: "QC using fastQC",
      desc: "Using fastQC 0.11.4 to generate quality control reports",      
      author: "djek.nad@gmail.com"

    var SOUT : "./FastQC"
    var threds: 8
    def basename = new File(input).getName().prefix;

	exec """
      mkdir -p ${SOUT}/${basename};	  
	  fastqc --threads ${threds} --extract  --outdir "${SOUT}/${basename}" "$input"
	"""
	forward input
}





Trimmomatic ={
      doc title: "Lunching Trimmomatic",
      desc: "Using Trimmomatic to remove leading and tailing low quality base-paires. Specify PE or SE for pair-end or Single-End sequences",
      author: "djek.nad@gmail.com"

    var SOUT : "./Trimmomatic"
    var threds: 10
    var TYPE: "PE"
    var adapterPE: "/home/nadhir/tools/trimmomatic/all-PE.fa"
    var adapterSE: "/home/nadhir/tools/trimmomatic/all-SE.fa"
    var MINLEN: 30
    def basename = new File(input1).getName().prefix
    def extension= (input1 =~ /[^\.]*$/)[0]


  if ( inputs.size() == 1){
      def basename1 = branch.name
      produce( "${SOUT}/${basename1}_trimmed.${extension}"){
        exec """
            mkdir -p ${basename1}_log; 
            mkdir -p ${SOUT};

            java -jar -Xmx30G $Trimmomatic_LOC SE
            -phred33
            $input
            $output
            ILLUMINACLIP:${adapterSE}:2:30:10:5:true
            LEADING:3
            TRAILING:3
            SLIDINGWINDOW:4:15
            CROP:100
            MINLEN:${MINLEN};
            if [ -f ${basename1}_trimmed.${extension} ]; then 
              mv ${basename1}_trimmed.${extension} ${SOUT}/${basename1}_trimmed.${extension};
            fi;
        ""","trimmomatic"
    }
  }
  else{
    
    def basename1 = new File(input1).getName().prefix
    def basename2 = new File(input2).getName().prefix
      produce( "${SOUT}/${basename1}_trimmed.${extension}","${SOUT}/${basename2}_trimmed.${extension}",
            "${SOUT}/${basename1}_unpaired.${extension}","${SOUT}/${basename2}_unpaired.${extension}"){
            exec """
                      mkdir -p ${SOUT};
                      java -jar -Xmx20G  $Trimmomatic_LOC PE
                        -phred33
                        $input1 $input2
                        $output1 $output3
                        $output2 $output4
                        ILLUMINACLIP:${adapterPE}:2:30:10:5:true
                        LEADING:3
                        TRAILING:3
                        SLIDINGWINDOW:4:15
                        MINLEN:${MINLEN};

                        if [ -f ${basename1}_trimmed.${extension} ]; then 
                          mv ${basename1}_trimmed.${extension} ${SOUT}/${basename1}_trimmed.${extension};
                          mv ${basename1}_unpaired.${extension} ${SOUT}/${basename1}_unpaired.${extension};
                        fi;

                        if [ -f ${basename2}_trimmed.${extension} ]; then 
                          mv ${basename2}_trimmed.${extension} ${SOUT}/${basename2}_trimmed.${extension};
                          mv ${basename2}_unpaired.${extension} ${SOUT}/${basename2}_unpaired.${extension};
                        fi;                       
            ""","trimmomatic"
    }
    forward "${SOUT}/${basename1}_trimmed.${extension}", "${SOUT}/${basename2}_trimmed.${extension}"
  }
}





STAR_mapping_old = {
    doc title: "Mapping and expression quantification",
    desc: "In this stage we use STAR to map sequencing data and do gene expression quantification",
      author: "djek.nad@gmail.com"
 
    produce("./${branch.name}/${branch.name}Aligned.sortedByCoord.out.bam"){
    exec """
     mkcdir -p ${branch.name} &&
     STAR --genomeDir ${GENOME}   
     --outSAMtype BAM SortedByCoordinate
     --sjdbGTFfile ${gtf}
     --readFilesIn $inputs 
     --runThreadN 15 --outFileNamePrefix "${branch.name}/${branch.name}"
     --outSAMunmapped Within --outFilterType BySJout 
     --outSAMattributes NH HI AS NM MD    
     --outFilterMultimapNmax 20
     --outFilterMismatchNmax 999
     --readFilesCommand zcat 
   """
   }
   forward "${branch.name}/${branch.name}Aligned.sortedByCoord.out.bam"
}





TagWithExons = {
      doc title: "Tag genomic reagions",
      desc: "Use dropseqtool to mark reads that map overlap with exons",
      author: "djek.nad@gmail.com"

        def AnnoFile="/home/nadhir/genomes/mm10/STAR_SynDux_tdTomato_INDEX/GRCm38.85_SynDux_TdTomato.refFlat"
        produce("${branch.name}/${branch.name}Aligned.sortedByCoord.out_exon_tagged.bam"){
                exec """
                        $DropSeqTools/TagReadWithGeneExon
                                        I=${input}
                                        O=${branch.name}/${branch.name}Aligned.sortedByCoord.out_exon_tagged.bam
                                        ANNOTATIONS_FILE=${AnnoFile}
                                        TAG=GE
                ""","TagWithExons"
        }
}


get_nonGenomicReads = {
  doc title: "Extract non-genomic reads names",
      desc: "Use a cutom python script to get the read names of the non-genomic reads",
      author: "djek.nad@gmail.com"

 produce("./NonGenomicReads/${branch.name}_IntergenicReads.bed"){
 exec """   
   mkdir -p NonGenomicReads &&
   python extract_notMappedToGenes.py ${input}  > ./NonGenomicReads/${branch.name}_IntergenicReads.bed
 """
 }
}



extractNonGenomicReads  = {
  doc title: "Extract the fastq sequence of the non-genomic reads",
      desc: "Use a cutom python script to generate the fastq files corresponding to the read names",
      author: "djek.nad@gmail.com"

  from(glob("./Trimmomatic/${branch.name}*gz",  "./NonGenomicReads/${branch.name}*.bed")){
  produce("./NonGenomic_fastq/${branch.name}_trimmed.nonGenomic.fq"){
  exec """
    mkdir -p NonGenomic_fastq &&
    zcat $input1 > ${input1.prefix}.fq &&
    seqtk subseq  ${input1.prefix}.fq ${input2} >  ./NonGenomic_fastq/${branch.name}_trimmed.nonGenomic.fq
  """
  }
 }
 forward "./NonGenomic_fastq/${branch.name}_trimmed.nonGenomic.fq"
}





indexBam = {
  produce("${input}.bai"){
  exec """
        samtools index $input.bam
  """
  }
  forward input
}




  
STAR_mapping = {
    doc title: "Mapping and expression quantification",
    desc: "In this stage we use STAR to map sequencing data and do gene expression quantification",
    author: "djek.nad@gmail.com"

    def ref_genome = GENOME
    def ref_gtf = gtf
    def outdir_pref= ""
    def outdir= "${branch.name}${outdir_pref}"
    produce("./${outdir}/${branch.name}${outdir_pref}Aligned.sortedByCoord.out.bam"){
    exec """
     mkdir -p ${outdir} &&
     STAR --genomeDir ${ref_genome}
     --outSAMtype BAM SortedByCoordinate
     --sjdbGTFfile ${ref_gtf}
     --readFilesIn $inputs
     --runThreadN 15 --outFileNamePrefix "${outdir}/${branch.name}"
     --outSAMunmapped Within --outFilterType BySJout
     --outSAMattributes NH HI AS NM MD
     --outFilterMultimapNmax 20
     --outFilterMismatchNmax 999
    --readFilesCommand zcat
     --quantMode TranscriptomeSAM GeneCounts ;
   ""","STAR_mapping"
   }
   forward "./${outdir}/${branch.name}${outdir_pref}Aligned.sortedByCoord.out.bam"
}





STAR_mapping_repeats = {
    doc title: "Mapping and expression quantification",
    desc: "In this stage we use STAR to map sequencing data and do gene expression quantification",
    author: "djek.nad@gmail.com"

    def ref_genome = REPEATS_GENOME
    def ref_gtf = REPEATS_GTF
    def outdir_pref = "_repeats_final"
    def outdir= "${branch.name}${outdir_pref}"
    produce("./${outdir}/${branch.name}Aligned.sortedByCoord.out.bam"){
    exec """
     mkdir -p ${outdir} &&
     STAR --genomeDir ${ref_genome}
     --outSAMtype BAM SortedByCoordinate
     --sjdbGTFfile ${ref_gtf}
     --readFilesIn $inputs
     --runThreadN 15 --outFileNamePrefix "${outdir}/${branch.name}"
     --outSAMunmapped Within --outFilterType BySJout
     --outSAMattributes NH HI AS NM MD
     --outFilterMultimapNmax 20
     --outFilterMismatchNmax 999
     --quantMode TranscriptomeSAM GeneCounts ;
   ""","STAR_mapping_repeats"
   }
   forward "./${outdir}/${branch.name}Aligned.sortedByCoord.out.bam"
}



// Entry point of the pipeline
run {    
   branches * [ fastqc +
                Trimmomatic +  
                STAR_mapping +               
                TagWithExons + 
                indexBam + 
                get_nonGenomicReads + 
                extractNonGenomicReads + 
                STAR_mapping_repeats                
             ] 
}


