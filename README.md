A small script to check if 

Prerequisites
-------------

You need one of trascription factor finders:
* tfscan from EMBOSS
* fimo from MEME suite

fimo requires a database -- it can be provided with `--db` option.

Usage
-----

You need a VCF input file and a FASTA file with reference genome.

The run the script like this:
    ruby snp_tf_change.rb -f ref.fa -v snp.vcf -m fimo --db HOCOMOCOv9_AD_MEME.txt > output.txt

Output would be like:
    1       100        rs1       C       G       -       gain    HEN1_si 9e-05    

Which means:
    allele G at this position creates a new TF binging site for 
    HEN1 on the reverse strand (which matches the sequence with pvalue 9e-05)
