#! /bin/bash

## This script is not intended to be run, it is primarily for reference of how the toy data was prepared ##
## This is a path on my local machine, you would need your own VEP installation to run this ##
VEP_DIR="/wdata/bcbio-dna/genomes/Hsapiens/hg19/vep"

THREADS=6
     
#### Download Chromosome 22  ####
VCF="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
VCF_IDX="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"
REF="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
REF_IDX="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz.fai"

wget $VCF $VCF_IDX $REF $REF_IDX

#### Set up file names ####
NORM_VCF=$(basename -s .vcf.gz $VCF).norm.vcf.gz
VEP_VCF=$(basename -s .vcf.gz $NORM_VCF).vep.vcf.gz

#### Normalize ####
bcftools norm -c ws -f $(basename $REF) -m - $(basename $VCF)  \
  | sed -e 's/Number=A/Number=1/g' \
  | head -n 1255 \
  | bgzip -c >  $NORM_VCF 
  
tabix $NORM_VCF 

#### Annotate with VEP ####
unset PERL5LIB && \
    vep --vcf \
    -o stdout \
    -i $(basename $NORM_VCF) \
    --fork $THREADS \
    --species homo_sapiens \
    --no_stats \
    --cache \
    --offline \
    --dir $VEP_DIR \
    --symbol \
    --numbers \
    --biotype \
    --total_length \
    --canonical \
    --gene_phenotype \
    --ccds \
    --uniprot \
    --domains \
    --regulatory \
    --protein \
    --tsl \
    --appris \
    --af \
    --max_af \
    --af_1kg \
    --af_esp \
    --pubmed \
    --variant_class \
    --allele_number \
    --fasta $(basename $REF) \
    --sift b \
    --polyphen b \
    --hgvs \
    --shift_hgvs 1 \
    --merged | \
    sed '/^#/! s/;;/;/g' | \
    bgzip -c > $VEP_VCF
    
tabix $VEP_VCF

mkdir anno-vcf
mv $(basename $VEP_VCF)* anno-vcf/

mkdir norm-vcf
mv $(basename $NORM_VCF)* norm-vcf/

mkdir raw-vcf
mv $(basename $VCF)* raw-vcf/
mv $(basename $REF)* raw-vcf/





