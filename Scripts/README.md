# Code for inital data handling of >6,000 whole cancer genomes

The data for this work is located at the GenomeDK cluster in the HMF project folder. The code is designed to merge and curate >12,000 separate .vcf files (>6,000 germline and >6,000 somatic) into a single file containing variants in genes-of-interest with both SnpEff and CADD v1.6 annotation.

## Filtering .VCF files
I have made symbolic links for all .vcf files, so that they mey all be accessed from a single folder. By running the 'Merge_VCF_files.sh' these are combined into a single, large, .vcf file in a subfolder called 'merged'

I have used VCFtools or TABIX to include only variants in the genes-of-interest. Either program gives the same results; VCFtools works on the raw vcf files whereas tabix demands a tabix-index, which takes time to generate. In brief, it is worth the while to generate the tabix-index if wanting to subset the data multiple times. 

## Assembling and curating variants
SnpEff and CADD are run on the cluster. To run these tools, you need to create the correct environment and install them (see guides at the respective webpages)

I have used R-scripts to do final assembly of the data.
