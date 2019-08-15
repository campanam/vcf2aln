# vcf2aln
Script to convert multi-sample VCFs to FASTA alignments without assuming the reference sequence when data are missing. Users can apply a variety of data filters, produce phased/unphased, concatenated/split alignments, etc. VCF data can be read either from previously generated files or from piped uncompressed VCF streams.    

## Authors
Michael G. Campana & Jacob A. West-Roberts, 2017-2019  

## License  
The software is made available under the Smithsonian Institution [terms of use](https://wwww.si.edu/termsofuse).  

## Citation  
Campana, M.G. & West-Roberts, J.A. 2019. vcf2aln. (https://github.com/campanam/vcf2aln)  

## Installation  
In the terminal:  
`git clone https://github.com/campanam/vcf2aln`  
`cd vcf2aln`  
`chmod +x vcf2aln.rb`  

Optionally, vcf2aln.rb can be placed within the user’s $PATH so that it can be executed from any location. Depending on your operating system, you may need to change the shebang line in the script (first line starting with #!) to specify the path of your Ruby executable.  

## Input  
vcf2aln requires an all-sites VCF (e.g. such as one produced using EMIT_ALL_SITES in the [Genome Analysis Toolkit](https://software.broadinstitute.org/gatk/)).  

## Execution
Execute the script using `ruby vcf2aln.rb` (or `vcf2aln.rb` if the script is in your $PATH). This will display the help screen. Basic usage is as follows:  
`ruby vcf2aln.rb -i <input_vcf> -o <out_prefix>` 

## Available options  
### I/O options:  
`-i, --input [FILE]`: Input VCF file.  
`--pipe`: Read data from an uncompressed VCF stream rather than a file.  
`-o, --outprefix [VALUE]`: Output FASTA alignment prefix.  
`-c, --concatenate`: Concatenate markers into single alignment (e.g. concatenate multiple separate chromosomes/contigs).  
`-s, --skip`: Skip missing sites in VCF.  
`-O, --onehap`: Print only one haplotype for diploid data. If phasing information is missing, it will generate a pseudohaplotype by randomly assigning one of the alleles.  
`-a, --alts`: Print alternate (pseudo)haplotypes in same file.  
`-b, --ambig`: Print SNP sites as ambiguity codes.  
`-N, --hap_flag`: Data are haploid.  
`-g, --split_regions [VALUE]`: Split alignment into subregional alignments for phylogenetic analysis. *DO NOT USE: UNDER DEVELOPMENT*  

### Filtration options:  
`-m, --mincalls [VALUE]`: Minimum number of individuals called to include site (Default = 0).  
`-x, --maxmissing [VALUE]`: Maximum percent missing data to include sequence (Default = 100.0).  
`-q, --qual_filter [VALUE]`: Minimum accepted value for QUAL (per site) (Default = 0.0).  
`-y, --site_depth [VALUE]`: Minimum desired total depth for each site (Default = No filter).  
`-d, --sampledepth [VALUE]`: Minimum allowed sample depth for each site (Default = No filter).  
`-l, --likelihood [VALUE]`: Minimum allowed genotype log-likelihood (At least one option must satisfy this value).  
`-p, --phred [VALUE]`: Minimum accepted phred-scaled genotype likelihood (Default = No filter).  
`-P, --posterior [VALUE]`: Minimum accepted phred-scaled genotype posterior probability (Default = No filter).  
`-C, --conditional [VALUE]`: Minimum conditional genotype quality (phred-encoded) (Default = No filter).  
`-H, --haplotype_quality [VALUE]`: Minimum allowed haplotype quality (phred-encoded) (Default = No filter).  
`-r, --sample_mq [VALUE]`: Minimum allowed per-sample RMS mapping quality (Default = No filter).  
`-R, --site_mq [VALUE]`: Minimum allowed per-site mapping quality (MQ in INFO) (Default = No filter).  
`-F, --mq0f [VALUE]`: Maximum allowed value for MQ0F. Must be between 0 and 1. (Default = No filter).  
`-S, --mqsb [VALUE]`: Minimum allowed value for MQSB. (Default = No filter).  
`-A, --adepth [VALUE]`: Minimum allowed allele depth. (Default = No filter).  

### General information:
`-t, --typefields`: Display VCF genotype field information, then quit the program.  
`-W, --writecycles`: Number of variants to store in memory before writing to disk. (Default = 1000000).  
`-v, --version`: Print program version.  
`-h, --help`: Show help.  
