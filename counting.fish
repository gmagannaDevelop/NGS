#!/usr/bin/env fish

## Université Paris-Saclay
## Atelier NGS
## Student : Gustavo Magaña López

# Locations

set --path reads Data/Reads
set --path trimming Data/Trimming
set --path genome Data/Genome
set --path index Data/Index
set --path mapping Data/Mapping
set --path counts Data/Counts
set --path figures_reads Figures/Reads
set --path figures_trimming Figures/Trimming

set threads 16

# Step 1: Quality control + Reads cleaning

mkdir -p $figures_reads
echo "Creation of the fastqc files on raw reads."
fastqc -t $threads -o $figures_reads -f fastq $reads/*.fastq -q

mkdir -p $trimming

for read1_file in $reads/*.R1.fastq
    set paired_file_with_path (string split .R1.fastq $read1_file)[1]
    set paired_file_without_path (string split "/" $paired_file_with_path)[-1]
    echo Trimming (string replace ".sampled" "" $paired_file_without_path)...
    set _R1 "$paired_file_with_path.R1.fastq"
    set _R2 "$paired_file_with_path.R2.fastq"
    set _baseout $trimming/$paired_file_without_path.fastq
    trimmomatic PE -threads $threads $_R1 $_R2 -baseout $_baseout  LEADING:20 TRAILING:20 MINLEN:50 -quiet
    echo "Done."
end

# Removal of the bases from the extremity with a quality lower than 20. 
# If the final read is smaller than 50, it is discarded. 
# file with U => discard. file with P => no discard.

mkdir -p $figures_trimming
echo "Creation of the fastqc files on trimmed reads."
fastqc -t $threads -o $figures_trimming -f fastq $trimming/*.fastq -q

# Step 2: Mapping

#### Index (STAR) ####

mkdir $index -p
STAR --runMode genomeGenerate --runThreadN $threads \
        --genomeSAindexNbases 12 \
        --genomeDir $index \
        --genomeFastaFiles $genome/chr18.fa \
        --sjdbGTFfile $genome/gencode.v24lift37.basic.annotation.gtf

#### Mapping (STAR) ####

mkdir $mapping -p
for read1_file in $trimming/*1P.fastq
    set paired_file_with_path (string split _1P.fastq $read1_file)[1]
    set paired_file_without_path (string split "/" $paired_file_with_path)[-1]
    echo -n "Creating BAM file with $paired_file_without_path..."; 
    set _1P "$paired_file_with_path"_1P.fastq
    set _2P "$paired_file_with_path"_2P.fastq
    STAR --runThreadN $threads --outFilterMultimapNmax 1 \
    --genomeDir $index \
    --outSAMattributes All --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix $mapping/$paired_file_without_path\_ \
    --readFilesIn $_1P $_2P;
    echo Done
end

# Step 3: Index BAM (samtools index)

for bam_file in $mapping/*.bam
        set bam_file_without_path (string split "/" $bam_file)[-1]
        echo "Indexing BAM file with $bam_file_without_path";
        samtools index -@ $threads $bam_file;
end

# Debug BAM file indexation: 
# samtools stats $mapping/BAM_file | less


# Step 4: Counting (featureCounts)


mkdir -p $counts
featureCounts -p -T $threads -t gene -g gene_id -s 0 -a $genome/*.gtf -o $counts/counts.txt $mapping/*.bam

# Create a file with pairs between ENCODE and HUGO identifiers
perl -ne 'print "$1 $2\n" if /gene_id \"(.*?)\".*gene_name \"(.*?)\"/' \
        $genome/*.gtf | sort | uniq > $counts/encode-to-hugo.tab

sort $counts/counts.txt > $counts/sort_counts.txt

# Remove the lines which does not contain counting data before merging encode-to-hugo.tab and sort_counts.txt
sed -i '/^[#|Geneid]/d' $counts/sort_counts.txt

# Creation of hugo-counts.txt file, 
# containing for each HUGO code in Chromosome 18, 
# the numbers of reads per gene and per observation.
set hugo_lines (cat $counts/encode-to-hugo.tab | wc -l)
set count_lines (cat $counts/sort_counts.txt | wc -l)
if test $hugo_lines -eq $count_lines
    echo -n "Creating hugo-counts.txt..."
    join $counts/encode-to-hugo.tab $counts/sort_counts.txt | grep "chr18" > $counts/paired_counts.txt
    awk '{print $2 " " $8 " " $9 " " $10 " " $11 " " $12 " " $13}' $counts/paired_counts.txt > $counts/hugo-counts.txt
    echo "Done."
else
    echo
    echo ERROR
    echo 'Cannot create hugo-counts.txt because files have different numbers of lines'
    exit 1
end

