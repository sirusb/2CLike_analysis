
# Download the position of repeats from RepEnrich (it uses mm9, but it is ok as we are interested at the sequence)

#Go to: https://github.com/nskvir/RepEnrich
#Download : mm9_repeatmasker_clean.txt

# Genereate pseudo-genome
python creatermskPseudo.py mm9_repeatmasker_clean.txt ../fasta/mm9_withRandom.fa rms_Pseudo_out

cd  rms_Pseudo_out

# Create gtf
Rscript RefFlatToGff.R

# Create dict

samtools faidx rsm_pseudo.fa
java -jar /nfs2/nadhir/tools/picard.jar CreateSequenceDictionary R=rsm_pseudo.fa O=rsm_pseudo.dict

# Create STAR index
mkdir Pseudogenome_STARINDEX
cd Pseudogenome_STARINDEX
STAR   --runMode genomeGenerate   --runThreadN 16   --genomeDir .   --genomeFastaFiles ../rsm_pseudo.fa      --outFileNamePrefix out


