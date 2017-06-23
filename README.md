# vcf_concat_mixed
Concatanate tabix indexed and sorted VCFs for the same samples the same contig may appear in multiiple VCFs 

## What?
In situations where you end up with multiple sorted VCFs but where the same chromosome may appear in more than one VCF (an example would be if you've lifted variants from one assembly to another but processed each chromosome separately) most simple concatanation tools will not deal with these variants properly. More complex tools (e.g. Picard or GATK) may be able to deal with these but can be overkill for a relatively simple task.

This tool assumes all your VCFs contain the same samples (it will use the header for the first VCF file given only in the output) and concatanates variants in order, handling chromosomes that appear multiple times in an appropriate manner, leaving you with a VCF that can be compressed and indexed with tabix.

## Usage

    python vcf_concat_mixed.py input1.vcf.gz input2.vcf.gz [input3.vcf.gz ...inputN.vcf.gz] 
    
Because this writes to STDOUT you will probably want to pipe the output to bgzip to create a compressed VCF

    python vcf_concat_mixed.py *.vcf.gz | bgzip -c > output.vcf.gz

You will need to install pysam if you haven't already.

## Author
David A. Parry, University of Edinburgh
