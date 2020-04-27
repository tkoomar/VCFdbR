# VCFdbR


This processing pipeline converts a VCF into an SQLite representation, using R. This readme only covers the practical matters of executing the pipeline. Please [see the wiki](https://github.com/tkoomar/VCFdbR/wiki)  or [FAQ](https://github.com/tkoomar/VCFdbR/wiki/FAQ) for more detailed discussions and tutorials covering how to use a VCFdb once it is created. 

[How do I convert a VCF to a VCFdb?](#how-do-i-convert-a-vcf-to-a-vcfdb)

[What do I need to run VCFdbR?](#what-do-i-need-to-run-vcfdbr)

[So what exactly does VCFdbR do?](#so-what-exactly-does-vcfdbr-do)

## A few important notes

VCFdbR was designed with Linux-based high performance computing clusters in mind, but it should work on any distro. It has not been tested on Windows or Mac operating systems, and changes to specifically support them are not planned. 

The main benefit of a SQLite database in this context is quick searching of huge datasets with minimal RAM overhead. This entire concept is based on the adage that _"storage is cheaper than time"_. No matter the options you choose, plan for the resulting database to take up at least 10 times as much disk space as the gzipped VCF used to create the database in the first place. 

## How do I convert a VCF to a VCFdb?

The only input required for VCFdbR is a [tabix-indexed](https://www.biostars.org/p/59492/) and properly formatted VCF, with *no multiallelic sites*. From the command line, you can call the master pipeline script `VCFdb.R` with `Rscript`. 

```{shell}
$ Rscript VCFdb.R --prefix [character] --vcf [character] --mode ['file'|'table']
```

* The string passed to `--prefix` will be the name of the output database and associated files
* The string passed to `--vcf` should be the path to the VCF from which the database will be created
* The `--mode` argument must be passed either the string `'file'` or `'table'`
    * `'table'` mode produces a "Table-GT" database where the genotypes are stored in a table within the SQLite database. 
      * This is ideal for smaller cohorts, as it allows for querying and filtering on genotypes and qualities without reading them into memory. It also keeps all the data together in a single file. 
      * *Unless you are dealing with thousands of whole-genomes or tens of thousands of exomes, `'table'` is likely the optimal mode.*
    * `'file'` mode produces a "File-GT" database where the genotypes are stored as individual files in a directory. The SQLite database will only contain the various variant annotations present in the VCF. 
        * This is ideal when working with very large cohorts that may result in final databases that are so large as to violate filesystem or SQLite database size limits (~4TB). However, the resulting database is significnatly less portable and care should be taken when choosing where to store the database. 
        * If your hard disk is fast, this has the added benefit of allowing genotypes to be written in parallel, via the `--threads` argument
        * This mode puts a column in the database pointing to the locaiton of the genotype files on the filesystem, this is why _**it is not reccomended to move the File-GT genotype folder after creation.**_ Consider where you build the database carefully. 
            * It is possible, but will require manually alterirng the paths in the `geno` column of the `variant_impact` table. Plus, moving the thousands or millions of individual genotype files will be exceptionally slow. 

## What do I need to run VCFdbR?

First, you need to [clone this repository](https://help.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository). 

Several R packages are required to be available. If you want to ensure you have all of the require packages, and that they are up to date, run the following code (in R) install these packages:
```{r}
> install.packages('tidyverse', 'dbplyr', 'magrittr', 'progress', 'DBI', 'RSQLite')
> if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
> BiocManager::install("VariantAnnotation")
```

If you are unfamiliar with executing R scripts from the command line via `Rscript`, I would suggest [reading up on it](https://support.rstudio.com/hc/en-us/articles/218012917-How-to-run-R-scripts-from-the-command-line).

Finally, you need a ([bgzipped and tabix indexed](https://davetang.org/muse/2013/02/22/using-tabix/)) VCF to convert! 

**This VCF needs to have all multialleleic sites split.** You can do that with `bcftools`:

```
$ bcftools norm -m - -O z -o [OUTPUT VCF] [YOUR VCF]
```

Also, if your VCF has been annotated with VEP, you might need to do a little munging in order for the CSQ column to be parsed correctly. This can be done as part of the VEP run, or afterwords, using `sed`:

```
$ zcat [YOUR VCF] | sed '/^#/\! s/;;/;/g' | bgzip -c > [OUTPUT VCF]
$ tabix [OUTPUT VCF]
```

## So what exactly does VCFdbR do?

VCFdbR takes data in the difficult-to-parse and search [Variant Call File](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format and converts it to a SQLite database that can be indexed, allowing for rapid searching. This makes exploratory analyis significnatly faster and removes the need to read the entire VCF into memory at once. This design was inspired by the fantastic [GEMINI](https://gemini.readthedocs.io/en/latest/), but has some critical changes to allow for processing very large cohorts and remvoes a lot of overhead by not providing its own bespoke interface to the dattabase. 

For example, a variant specified like this in a VCF:
```{text}
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096 HG00097 HG00099
22      16120773        rs577167963     G       A       100     PASS    AC=1;AF=0.000199681;AN=5008;NS=2504;DP=20455;EAS_AF=0;AMR_AF=0;AFR_AF=0;EUR_AF=0.001;SAS_AF=0;AA=G|||;VT=SNP;CSQ=A|upstream_gene_variant|MODIFIER|LA16c-60H5.7|ENSG00000215270|Transcript|ENST00000398242|processed_pseudogene||||||||||rs577167963|1|1947|1||SNV|Clone_based_vega_gene||YES|||||||||Ensembl|G|G|||||||0.0002|0|0|0|0.001|0|||0.001|MODIFIER|NBEAP3|ENSG00000223875|Transcript|ENST00000420638|unprocessed_pseudogene||1/3|ENST00000420638.1:n.233-1860C>T|YES|||||||||Ensembl|G|G|||||||0.0002|0|0|0|0.001|0|||0.001|EUR||||||||      GT      0|0     0|0      0|0
```

Would have its information divided among three tables in the resulting SQLite database, all of which would be linked by an indexed variant number.

Information from the core variant information and `INFO` column would be located on the `variant_info` table, where each varaint exists exaclty once:

|**variant_id**|chr|start|end|ref|alt|qual|filter|ac|af|ns|an|eas_af|eur_af|afr_af|amr_af|sas_af|dp|vt|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|**843**|22|16120773|16120773|G|A|100|PASS|1|10.000199681|2504|5008|0|0.001|0|0|0|20455|SNP|

Information specifically about the variant's consequence (generated by VEP) would be located on the `variant_impact` table, where it might have multiple entries corresponding to multiple effects:

|**variant_id**|consequence|impact|symbol|gene|biotype|exon|intron|existing_variation|...
|---|---|---|---|---|---|---|---|---|---
|**843**|upstream_gene_variant|MODIFIER|LA16c-60H5.7|ENSG00000215270|processed_pseudogene|||rs577167963|...
|**843**|intron_variant|MODIFIER|NBEAP3|ENSG00000223875|unprocessed_pseudogene||1/3|rs577167963|...
|**843**|non_coding_transcript_variant|MODIFIER|NBEAP3|ENSG00000223875|unprocessed_pseudogene||1/3|rs577167963|...

Finally, the genotypes themselves would be located on the `varint_geno` table, where they are in "long" format. This has the advantage of allowign variants to be quickly filered by SQL based on genotype or quality score (not present in this example):

|**variant_id**| sample| gt| gt_raw|
|---|---|---|---|
|**843**| HG00096 | 0  |  0\|0|
|**843**| HG00097 | 0\|0|
|**843**| HG00099|0\|0|
