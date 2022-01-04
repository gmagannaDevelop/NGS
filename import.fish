#!/usr/bin/env fish

# Directory parameters

set --path reads Data/Reads
set --path genome Data/Genome

# Url parameters

set reads_url 'http://rssf.i2bc.paris-saclay.fr/X-fer/AtelierNGS/TPrnaseq.tar.gz'
set genome_url 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz'
set annotation_url 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz'

# Step 1: Download reads

mkdir -p $reads
echo "Downloading reads..."
wget $reads_url -P $reads
tar -zxvf $reads/TPrnaseq.tar.gz -C $reads

echo "--------Number of sequences per file--------"
for file in $reads/*.fastq
    echo "$file contains :";
    echo -n  (cat $file | wc -l) 'lines and ';
    echo (grep '^+$' $file | wc -l) entries;
end

# Step 2: download reference genome

mkdir $genome -p
echo -n "Downloading genome..."
wget $genome_url -P $genome -q
echo Done
echo -n "Downloading annotations..."
wget $annotation_url -P $genome -q
echo Done
echo "Decompress files..."
gunzip -f $genome/*.gz
echo Done