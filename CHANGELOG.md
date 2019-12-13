# vcf2aln Change Log  
Michael G. Campana & Jacob A. West-Roberts, 2017-2019  
Smithsonian Conservation Biology Institute  
Contact: campanam@si.edu  

### Version 0.11.1  
Added method fix_name to resolve reserved characters in locus names  

### Version 0.11.0  
Read/write gzipped files  

### Version 0.10.0  
Probabilistic pseudohaplotype calling  
Minimum samples called as a percentage option  
Minimum alignment length to retain parameter  

### Version 0.9.0  
Sample-specific filters no longer filter whole line  
Completely filtered sites are excluded by the skip option  
Renamed filters to better match standard tags  

### Version 0.8.0  
--annotfilter option controls FILTER value filtration  
--split_regions is now functional  

### Version 0.7.0  
vcf2aln can read streamed uncompressed VCF  
Removed extraneous debugging output from get_GT_tags  

### Version 0.6.0  
get_GT_tags gets tag information from VCF rather than VCF 4.2 standard tags  
Can read GATK HaplotypeCaller PGT phasing information  
GT/PGT information do not need to be in first slot of output   

### Version 0.5.0  
Speed increase using write-cycle controls  
Indexing bug fix  

### Version 0.4.2  
Haploid VCF bug fix  

### Version 0.4.1  
Onehap bug fix  

### Version 0.4.0  
Onehap concatenation bug fix  

### Version 0.3.0  
Now gets all type fields in VCF  
Onehap flag and bug fixes  
GLE & ambiguity code handling  
Cleaned up help screen  

### Version 0.2.0  
New filters  
Ability to identify VCF tags  
Separation of methods  

### Version 0.1.0  
Preliminary script to generate FASTA alignment from multi-sample VCF  

