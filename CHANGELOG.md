# vcf2aln Change Log  
Michael G. Campana & Jacob A. West-Roberts, 2017-2025  
Smithsonian's National Zoo and Conservation Biology Institute  
Contact: campanam@si.edu  

### Version 0.13.4  
Fixed a bug in the region splits that caused a single bp alignment at the changing of chromosomes  

### Version 0.13.3  
Added -V option to remove invariant sites after MSA conversion  

### Version 0.13.2  
Added --inferref option to infer reference when base is missing  

### Version 0.13.1  
Fixed bug in reference-inclusion if missing sites are not skipped  

### Version 0.13.0  
Added option to output reference sequence as part of the alignment  

### Version 0.12.1  
Fixed glitch where the region ID was added to the sequence length value  
Fixed glitch outputting a zero-length partition at beginning of table  

### Version 0.12.0  
Added option to output contig partition table for concatenated alignments  

### Version 0.11.5  
Handling for diploid missing data calls in haploid vcfs  
encdoed-typo correction  

### Version 0.11.4  
Fixed zlib requirement bug  

### Version 0.11.3  
Now compatible with GATK deletions coded by *  
Improvement to array summation  

### Version 0.11.2  
Fixed bug when splitting regions by length that would cause misnumbering and overwriting of previous regions  
Fixed bug for haploid data that did not properly count haploid missing data in the minimum calls filter  

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

