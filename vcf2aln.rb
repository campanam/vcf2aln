#!/usr/bin/env ruby

#-----------------------------------------------------------------------------------------------
# vcf2aln
VCF2ALNVER = "0.12.1"
# Michael G. Campana, Jacob A. West-Roberts, 2017-2023
# Smithsonian's National Zoo and Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

require 'optparse'
require 'ostruct'
require 'zlib'

class Locus
	attr_accessor :name, :seqs, :alts, :length
	#-------------------------------------------------------------------------------------------
	def initialize(name = "", seqs = [], alts = [], length = 0)
		@name = name
		@seqs = seqs
		@alts = alts
		@length = length
		@missinghap1 = []
		@missinghap2 = []
		unless @name == ""
			for sample in $samples
				if $options.onehap
					File.open(sample + "_#{$options.outprefix}#{@name}.tmp.fa", 'w') do |write|
						write.puts ">" + sample
					end
				else
					File.open(sample + "_#{$options.outprefix}#{@name}.hap1.tmp.fa", 'w') do |write|
						write.puts ">" + sample + "_hap1"
					end
					File.open(sample + "_#{$options.outprefix}#{@name}.hap2.tmp.fa", 'w') do |write|
						write.puts ">" + sample + "_hap2"
					end
				end
				@missinghap1.push(0)
				@missinghap2.push(0)
			end
		end
	end
	#-------------------------------------------------------------------------------------------
	def write_seqs
		@length += @seqs[0].length
		for i in 0...$samples.size
			if $options.onehap
				File.open($samples[i] + "_#{$options.outprefix}#{@name}.tmp.fa", 'a') do |write|
					write << @seqs[i]
				end
			else
				File.open($samples[i] + "_#{$options.outprefix}#{@name}.hap1.tmp.fa", 'a') do |write|
					write << @seqs[i]
				end
				File.open($samples[i] + "_#{$options.outprefix}#{@name}.hap2.tmp.fa", 'a') do |write|
					write << @alts[i]
				end
			end
			if $options.maxmissing < 100.0
				@missinghap1[i] += @seqs[i].scan("?").length
				unless $options.hap_flag
					@missinghap2[i] += @seqs[i].scan("?").length
				end
			end
			@seqs[i] = ""
			@alts[i] = "" unless $options.hap_flag
		end
	end
	#-------------------------------------------------------------------------------------------
	def print_locus
		@length += @seqs[0].length
		if @length >= $options.minlength
			for i in 0...$samples.size
				if $options.onehap
					File.open($samples[i] + "_#{$options.outprefix}#{@name}.tmp.fa", 'a') do |write|
						write.puts @seqs[i]
					end
				else
					File.open($samples[i] + "_#{$options.outprefix}#{@name}.hap1.tmp.fa", 'a') do |write|
						write.puts @seqs[i]
					end
					File.open($samples[i] + "_#{$options.outprefix}#{@name}.hap2.tmp.fa", 'a') do |write|
						write.puts @alts[i]
					end
				end
				if $options.maxmissing < 100.0
					@missinghap1[i] += @seqs[i].scan("?").length
					@missinghap2[i] += @seqs[i].scan("?").length unless $options.hap_flag
					if $options.onehap
						system("rm #{$samples[i] + '_' + $options.outprefix + @name}.tmp.fa") if @missinghap1[i].to_f/@length.to_f * 100.0 > $options.maxmissing
					else
						system("rm #{$samples[i] + '_' + $options.outprefix + @name}.hap1.tmp.fa") if @missinghap1[i].to_f/@length.to_f * 100.0 > $options.maxmissing
						system("rm #{$samples[i] + '_' + $options.outprefix + @name}.hap2.tmp.fa") if @missinghap2[i].to_f/@length.to_f * 100.0 > $options.maxmissing
					end
				end
			end
			if $options.alts
				system("cat *#{$options.outprefix + @name}*.tmp.fa > #{$options.outprefix + @name + '.fa'}")
				system("gzip #{$options.outprefix + @name + '.fa'}") if $options.gzip
			else
				if $options.onehap
					system("cat *#{$options.outprefix + @name}.tmp.fa > #{$options.outprefix + @name + '.fa'}")
					system("gzip #{$options.outprefix + @name + '.fa'}") if $options.gzip
				else
					system("cat *#{$options.outprefix + @name}.hap1.tmp.fa > #{$options.outprefix + @name + '.hap1.fa'}")
					system("cat *#{$options.outprefix + @name}.hap2.tmp.fa > #{$options.outprefix + @name + '.hap2.fa'}")
					system("gzip #{$options.outprefix + @name + '.hap1.fa'}") if $options.gzip
					system("gzip #{$options.outprefix + @name + '.hap2.fa'}") if $options.gzip
				end
			end
		end
		system("rm *#{$options.outprefix + @name}*.tmp.fa")
	end
	#-------------------------------------------------------------------------------------------
	def write_partitions
		@partstart ||= 1 # Start bp of the partition
		@partition ||= 1 # Partition ID
		File.open("#{$options.outprefix}#{@name}.partitions", 'a') do |write|
			write.puts 'DNA, part' + @partition.to_s + ' = ' + @partstart.to_s + "-" + @length.to_s
		end
		@partstart = @length + 1
		@partition += 1
	end
end
#-----------------------------------------------------------------------------------------------
class Parser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		args.infile = "" # Primary input file
		args.pipe = false # Read data from a pipe rather than an input file
		args.outprefix = "" # Output prefix
		args.hap_flag = false
		args.concat = false # Concatenate markers into single alignment
		args.partition = false # Output partition file for concatenated alignment
		args.skip = false # Skip missing sites in vcf
		args.mincalls = 0 # Minimum number of sample calls to include site
		args.minpercent = 0.0 # Minimum percent of sample calls to include site
		args.maxmissing = 100.0 # Maximum percent missing data to include sequence
		args.minlength = 1 # Minimum length of alignment to retain
		args.onehap = false # Print only one haplotype
		args.probps = false # Probabilistic pseudohaplotype
		args.alts = false # Print alternate haplotypes in same file
		args.ambig = false # Print SNPs as ambiguity codes
		args.qual_filter = 0.0 # Minimum quality for site (QUAL column)
		args.annot_filter = [] # Annotations to exclude
		args.site_depth = nil #Minimum site coverage depth
		args.type_fields = false #Don't display VCF genotype fields on default
		args.sample_depth = 0 #Don't filter VCF calls based on depth by default
		args.min_ll = nil #Don't filter likelihood on default
		args.phred_likelihood = nil #Don't filter phred-scaled likelihoods (PL field) on default
		args.posterior = nil #Don't filter posterior likelihood (GP field) on default
		args.conditional = nil #Don't filter conditional likelihood (GQ field) on default
		args.hap_qual = nil #Don't filter based on haplotype quality on default
		args.sample_mq = nil #Don't filter based on per-sample mapping quality
		args.site_mq = nil #Don't filter based on per-site mapping quality
		args.mq0f = nil #Don't filter MQ0F or MQSB by default
		args.mqsb = nil # ^^
		args.adepth = nil #Don't filter allele depth by default
		args.split_regions = 0 #Don't activate region-split subroutine by default
		args.write_cycle = 1000000 # Number of cycles before force of write-out
		args.gzip = false # Gzip output
		opt_parser = OptionParser.new do |opts|
			opts.banner = "Command-line usage: ruby vcf2aln.rb [options]"
			opts.separator ""
			opts.separator "I/O options:"
			opts.on("-i","--input [FILE]", String, "Input VCF file") do |vcf|
				args.infile = vcf
			end
			opts.on("--pipe", "Read data from pipe") do
				args.pipe = true
			end
			opts.on("-o", "--outprefix [VALUE]", String, "Output alignment prefix") do |pref|
				args.outprefix = pref + "_" if pref != nil
			end
			opts.on("-z", "--gzip", "Gzip output alignments") do 
				args.gzip = true
			end
			opts.on("-c", "--concatenate", "Concatenate markers into single alignment") do
				args.concat = true
			end
			opts.on("--partition", "Partition concatenated alignment by contig") do
				args.partition = true
			end
			opts.on("-s", "--skip", "Skip missing sites in VCF") do
				args.skip = true
			end
			opts.on("-O", "--onehap", "Print only one haplotype for diploid data (conflicts with -a)") do
				args.onehap = true
				args.alts = false
			end
			opts.on("--probpseudohap", "Probabilistic pseudohaplotype based on allelic depth (implies -O, conflicts with -a, -b)") do
				args.probps = true
				args.onehap = true
				args.alts = false
				args.ambig = false
			end
			opts.on("-a", "--alts", "Print alternate haplotypes in same file (conflicts with -O, --probpseudohap)") do
				args.alts = true unless (args.probps or args.onehap)
			end
			opts.on("-b", "--ambig", "Print SNP sites as ambiguity codes (conflicts with --probpseudohap)") do
				args.ambig = true unless args.probps
			end
			opts.on("-N", "--hap_flag", "Flag for haploid data") do 
				args.hap_flag = true
				args.onehap = true # Condense file writing
			end
			opts.on("-g", "--split_regions [VALUE]", Integer, "Split alignment into subregional alignments for phylogenetic analysis") do |regions|
				args.split_regions = regions if regions != nil
			end
			opts.separator ""
			opts.separator "Filtration options:"
			opts.on("-m","--mincalls [VALUE]", Integer, "Minimum number of samples called to include site (Default = 0)") do |msnps|
				args.mincalls = msnps if msnps != nil
			end
			opts.on("-M", "--minpercent [VALUE]", Float, "Minimum percentage of samples called to include site (Default = 0.0)") do |mpc|
				args.minpercent = mpc if mpc != nil
			end
			opts.on("-x","--maxmissing [VALUE]", Float, "Maximum sample percent missing data to include sequence (Default = 100.0)") do |missing|
				args.maxmissing = missing if missing != nil
			end
			opts.on("-L", "--minlength [VALUE]", Integer, "Minimum alignment length to retain (Default = 1)") do |len|
				args.minlength = len if len != nil
			end
			opts.on("--annotfilter [VALUE]", String, "Comma-separated list of FILTER annotations to exclude") do |annot|
				args.annot_filter = annot.split(",") if annot != nil
			end
			opts.on("-q", "--qual_filter [VALUE]", Integer, "Minimum accepted value for QUAL (per site) (Default = 0.0)") do |qual|
				args.qual_filter = qual if qual != nil
			end
			opts.on("-y", "--site_depth [VALUE]", Integer, "Minimum desired total depth for each site (Default = No filter)") do |site|
				args.site_depth = site if site != nil
			end
			opts.on("-d", "--sampledepth [VALUE]", Integer, "Minimum allowed sample depth for each site (Default = No filter)") do |depth|
				args.sample_depth = depth if depth != nil
			end
			opts.on("-l", "--gl [VALUE]", Float, "Minimum allowed genotype log-likelihood (tag GL). (Default = No filter)") do |likelihood|
				args.min_ll = likelihood if likelihood != nil
			end
			opts.on("-p", "--pl [VALUE]", Integer, "Minimum accepted phred-scaled genotype likelihood (tag PL). (Default = No filter)") do |phreddy|
				args.phred_likelihood = phreddy if phreddy != nil
			end
			opts.on("-G", "--gp [VALUE]", Float, "Minimum accepted phred-scaled genotype posterior probability (tag GP). (Default = No filter)") do |post|
				args.posterior = post if post != nil
			end
			opts.on("-C", "--gq [VALUE]", Float, "Minimum conditional phred-encoded genotype quality (tag GQ). (Default = No filter)") do |condi|
				args.conditional = condi if condi != nil
			end
			opts.on("-H", "--hq [VALUE]", Integer, "Minimum allowed phred-encoded haplotype quality (tag HQ). (Default = No filter)") do |haplo|
				args.hap_qual = haplo if haplo != nil
			end
			opts.on("-r", "--sample_mq [VALUE]", Integer, "Minimum allowed per-sample RMS mapping quality (tag MQ). (Default = No filter)") do |map_quality|
				args.sample_mq = map_quality if map_quality != nil
			end
			opts.on("-R", "--site_mq [VALUE]", Integer, "Minimum allowed per-site mapping quality (MQ in INFO). (Default = No filter)") do |mq|
				args.site_mq = mq if mq != nil
			end
			opts.on("-F", "--mq0f [VALUE]", Float, "Maximum allowed value for MQ0F. Must be between 0 and 1. (Default = No filter)") do |mqf|
				args.mq0f = mqf if mqf != nil
			end
			opts.on("-S", "--mqsb [VALUE]", Float, "Minimum allowed value for MQSB. (Default = No filter)") do |sb|
				args.mqsb = sb if sb != nil
			end
			opts.on("-A", "--ad [VALUE]", Integer, "Minimum allowed allele depth (tag AD). (Default = No filter)") do |ad|
				args.adepth = ad if ad != nil
			end
			opts.separator ""
			opts.separator "General information:"
			opts.on("-t", "--typefields", "Display VCF genotype field information, then quit the program.") do
				args.type_fields = true
				args.pipe = false # Prevent crash from reading from pipe
			end
			opts.on("-W", "--writecycles", Integer, "Number of variants to store in memory before writing to disk. (Default = 1000000)") do |wrt|
				args.write_cycle = wrt if wrt != nil
			end
			opts.on("-v", "--version", "Print program version.") do
				abort("vcf2aln v." + VCF2ALNVER + "\n")
			end
			opts.on_tail("-h","--help", "Show help") do
				puts opts
				exit
			end
		end
		opt_parser.parse!(options)
		return args
	end
end
#-----------------------------------------------------------------------------------------------
def fix_name(name) # Modified from BaitsTools 1.6.5 (Campana 2019) resolve_unix_path
	for reschar in ["\\", ";", "&", "(", ")","*","?","[","]","~",">","<","!","\"","\'", "$", " ", "|"]
		name.gsub!(reschar, "_") # Use odd syntax because escape character causes issues with backslash in sub
	end
	return name
end
#-----------------------------------------------------------------------------------------------
def quality_filter(line_arr, gt_index, pgt_index)
	site_qual = line_arr[5]
	filter = line_arr[6]
	site_info_fields = line_arr[7]
	info_arr = site_info_fields.split(";")
	sample_info_fields = line_arr[8].split(":")
	samples = line_arr[9..-1]
	samples.map!{ |element| element.split(":")}
	$options.hap_flag ? missingdata = "." : missingdata = "./." # Some VCFs are haploid
	genotypes = []
	if (gt_index.nil? && pgt_index.nil?)
		puts "** Missing genotype (GT or PGT) tags **"
		puts "** Treating line: #{line_arr.join("\t")} as missing data **"
		samples.each{|list| genotypes.push(missingdata)}
		new_samples = []
		samples.each{|list|
			list.insert(0, missingdata)
			new_list = list.join(":")
			new_samples.push(new_list)
		}
		line_arr[9..-1] = new_samples		
	else
		samples.each{|list| genotypes.push(list[gt_index])}
	end
	# Leave if nothing to replace
	if genotypes.all? {|x| x == missingdata}
		if $options.skip
			return "all_filtered"
		else
			return line_arr
		end
	end
	found = false # Filter ALL values from a site
	if sample_info_fields.include?("GLE")
		puts "** GLE not supported as of version #{VCF2ALNVER} **"
		puts "** Treating line: #{line_arr.join("\t")} as missing data **"
		found = true
	end
	found = true if site_qual.to_i < $options.qual_filter || ($options.annot_filter.include?(filter))
	unless ($options.mq0f.nil? && $options.mqsb.nil? && $options.site_mq.nil? && $options.site_depth.nil?) || found #Don't execute this loop if you're not trying to filter for any site-specific quality scores
		#This loop will determine if any of your site-specific quality scores don't satisfy the specified conditions. (INFO column)
		info_arr.each do |item|
			if found
				next
			else
				item_info = item.split("=")
			  	id = item_info[0]
			  	if id == "DP" && !$options.site_depth.nil?  #Scan for unacceptable site depth
			    	dp = item_info[1]
			    	found = true if dp.to_i < $options.site_depth
			  	elsif id == "MQ0F" && !$options.mq0f.nil?  #Scan for unacceptable mq0f
			   		mqf = item_info[1]
			    	found = true if mqf.to_f > $options.mq0f
			  	elsif id == "MQSB" && !$options.mqsb.nil? #Scan for unaccepable mqsb
			    	mqsb = item_info[1]
			    	found = true if mqsb.to_f < $options.mqsb
			  	elsif id == "MQ" && !$options.site_mq.nil? #Scan for unacceptable overall map quality
			    	mq = item_info[1]
			    	found = true if mq.to_i < $options.site_mq
			  	end
			end
		end
	end
	if found
		new_samples = []
		samples.each{|list|
			list[gt_index] = missingdata
			list[pgt_index] = missingdata unless pgt_index.nil?
			new_list = list.join(":")
			new_samples.push(new_list)
		}
		line_arr[9..-1] = new_samples
		if $options.skip
			return "all_filtered"
		else
			return line_arr
		end
	else
		unless ($options.sample_mq.nil? && $options.sample_depth.nil? && $options.min_ll.nil? && $options.phred_likelihood.nil? && $options.posterior.nil? && $options.conditional.nil? && $options.hap_qual.nil? && $options.adepth.nil?)
			index_hash = {}
			new_samples = []
			new_genotypes = []
			#Create index hash; info fields may be ordered in any way, so we need a generalizable method to access each field at its proper position in the sample info
			sample_info_fields.each{|v| index_hash[v] = sample_info_fields.index(v)}
			#This bit can be sped up by removing options from the index hash if the corresponding value in $options does not exist. 
			for sample in samples
				for field in sample_info_fields
					break if sample[gt_index] == missingdata
					case field
					when "GT", "PGT", "PQ"
						next
					when "DP"
						if $options.sample_depth && sample[index_hash["DP"]].to_i < $options.sample_depth
						 	sample[gt_index] = missingdata
							sample[pgt_index] = missingdata unless pgt_index.nil?
					 	end
					when "GQ"
						if $options.conditional && sample[index_hash["GQ"]].to_i < $options.conditional
							sample[gt_index] = missingdata
							sample[pgt_index] = missingdata unless pgt_index.nil?
						end
					when "MQ"
						if $options.sample_mq && sample[index_hash["MQ"]].to_i < $options.sample_mq
							sample[gt_index] = missingdata
							sample[pgt_index] = missingdata unless pgt_index.nil?
						end
					when "GL"
			  			if $options.min_ll
					 		gl_array = sample[index_hash["GL"]].split(",").map { |x| x.to_i }
					 		if gl_array.max < $options.min_ll
								sample[gt_index] = missingdata
								sample[pgt_index] = missingdata unless pgt_index.nil?
							end
						end
					when "PL"
						if $options.phred_likelihood 
							pl_array = sample[index_hash["PL"]].split(",").map { |x| x.to_i }
							if pl_array.max < $options.phred_likelihood
								sample[gt_index] = missingdata
								sample[pgt_index] = missingdata unless pgt_index.nil?
							end
						end
					when "GP"
						if $options.posterior
							gp_array = sample[index_hash["GP"]].split(",").map { |x| x.to_f }
							if gp_array.max < $options.posterior
								sample[gt_index] = missingdata
								sample[pgt_index] = missingdata unless pgt_index.nil?
							end
						end
					when "HQ"
						if $options.hap_qual
							hap_array = sample[index_hash["HQ"]].split(",").map { |x| x.to_i }
							if $options.hap_flag
								if hap_array[0] < $options.hap_qual
									sample[gt_index] = missingdata
									sample[pgt_index] = missingdata unless pgt_index.nil?
								end
							else
								if hap_array[0] < $options.hap_qual && hap_array[1] < $options.hap_qual # If both under, set as missing data
									sample[gt_index] = missingdata
									sample[pgt_index] = missingdata unles pgt_index.nil?
								elsif hap_array[0] < $options.hap_qual # If only first under, set as second haplotype
									sample[gt_index][0] = sample[gt_index][2]
									sample[pgt_index][0] = sample[ggt_index][2] unless pgt_index.nil?
								elsif hap_array[1] < $options.hap_qual # if only second under, set as first haplotype
									sample[gt_index][2] = sample[gt_index][0]
									sample[pgt_index][2] = sample[ggt_index][0] unless pgt_index.nil?
								end
							end
						end
					when "AD"
						if $options.adepth 
							ad_array = sample[index_hash["AD"]].split(",").map { |x| x.to_i }
							if $options.hap_flag
								if ad_array[sample[gt_index].to_i] < $options.adepth # If selected allele low, discard
									sample[gt_index] = missingdata
									sample[pgt_index] = missingdata unless pgt_index.nil?
								end
							else
								if ad_array[sample[gt_index][0].to_i] < $options.adepth && ad_array[sample[gt_index][3].to_i] < $options.adepth # If both under, set as missing data
									sample[gt_index] = missingdata
									sample[pgt_index] = missingdata unless pgt_index.nil?
								elsif ad_array[sample[gt_index][0].to_i] < $options.adepth # If only first under, set as second haplotype
									sample[gt_index][0] = sample[gt_index][2]
									sample[pgt_index][0] = sample[ggt_index][2] unless pgt_index.nil?
								elsif ad_array[sample[gt_index][2].to_i] < $options.adepth # if only second under, set as first haplotype
									sample[gt_index][2] = sample[gt_index][0]
									sample[pgt_index][2] = sample[ggt_index][0] unless pgt_index.nil?
								end
							end
						end
					end
				end
				new_sample = sample.join(":")
				new_genotypes.push(sample[gt_index])
				new_samples.push(new_sample)
			end
		end
		line_arr[9..-1] = new_samples
	end
	if new_genotypes.all? { |x| x== missingdata }
		return "all_filtered"
	else
		return line_arr
	end
end
#-----------------------------------------------------------------------------------------------
def build_ambig_hash
	$ambig_hash = {['A','G'] => 'R', ['C','T'] => 'Y', ['A','C'] => 'M', ['G','T'] => 'K', ['C','G'] => 'S', ['A','T'] => 'W'}
end
#-----------------------------------------------------------------------------------------------
def get_ambiguity_code(var1, var2)
	if var1 == var2
		return var1
	else
		return $ambig_hash[[var1,var2].sort]
	end
end
#-----------------------------------------------------------------------------------------------
def get_GT_fields(vcf_file)
	@fields = {}
	File.open(vcf_file, 'r') do |vcf2aln|
		while line = vcf2aln.gets
			if line[0..12] == "##FORMAT=<ID="
				fields_arr = line.split("Description")
				field = fields_arr[0].split(",Number")[0][13..-1]
				description = fields_arr[1][2...-3]
				@fields[field] = description
			end
		end
	end
	puts "Genotype field information categories:"
	for field in @fields.keys
		puts field + ": " + @fields[field]
		puts "** GLE containing vcf files not supported as of version #{VCF2ALNVER} **" if field == "GLE"
	end
	exit	
end
#-----------------------------------------------------------------------------------------------
def vcf_to_alignment(line, index, previous_index, previous_endex, previous_name, current_locus, prev_pos, write_cycle, region, regionval)
	if line[0..1] == "#C"
		$samples = line[0..-2].split("\t")[9..-1] # Get sample names
		minpc = ($samples.size.to_f * $options.minpercent / 100.0).ceil
		$options.mincalls = minpc if minpc > $options.mincalls
	elsif line[0].chr != "#"
		line_arr = line[0...-1].split("\t") # Exclude final line break
		codes = line_arr[8].split(":") #GET PHASING LOCATION
		pgt_index = nil
		gt_index = nil
		ad_index = nil
		if codes.include?("PGT")
			vars_index = codes.index("PGT")
			pgt_index = vars_index
			gt_index = codes.index("GT")
		elsif codes.include?("GT")
			vars_index = codes.index("GT")
			gt_index = vars_index
		else
			vars_index = 0
		end
		if $options.adepth
			if codes.include?("AD")
				ad_index = codes.index("AD")
			else
				puts "** Site missing AD tag **"
				puts "** Treating haplotypes as equally likely for line: #{line} **"
			end
		end
		line_arr = quality_filter(line_arr, gt_index, pgt_index) # Quality filter has to be called to check for GLE
		return index, previous_index, previous_endex, previous_name, current_locus, prev_pos, write_cycle, region, regionval if line_arr == "all_filtered"
		if $options.split_regions > 0 && regionval >= $options.split_regions # Change region name at break
			region += 1
			regionval = 0
		end
		regionval += 1
		write_cycle += 1
		$options.concat ? name = "concat_aln" : name = fix_name(line_arr[0].dup)
		name << "_region" + region.to_s if $options.split_regions > 0
		if name != current_locus.name
			current_locus.print_locus unless current_locus.name == ""
			seqs = []
			alts = []
			$samples.size.times do
				seqs.push("")
				alts.push("")
			end
			current_locus = Locus.new(name, seqs, alts)
			index = -1 # Set internal start index
			previous_index = -1 # Index of previous ending base
			previous_endex = 0 # End index of indels
			regionval = 1
		end
		if line_arr[0] != previous_name # Reset indexes for concatenated alignments
			if current_locus.name != ""
				current_locus.write_seqs
				current_locus.write_partitions if $options.partition
			end
			index = -1
			write_cycle = 1
			previous_index = -1
			previous_endex = 0
			previous_name = line_arr[0]
			region = 1
		end
		if $options.hap_flag
			missingness = 0
			for genotype in line_arr[9..-1]
				missingness += 1 if genotype.split(":")[gt_index][0] == "."
			end
		else
			missingness = line.scan("./.").length + line.scan(".|.").length
		end
		if $samples.size - missingness >= $options.mincalls
			variants = [line_arr[3], line_arr[4].split(",")].flatten!
			variants[variants.index("*")] = "-" if variants.include?("*")
			lengths =[] # Get sequence lengths for indels
			for var in variants
				lengths.push(var.length)
			end
			for i in 0...variants.size
				var = variants[i]
				if var.length < lengths.max
					(lengths.max-var.length).times {var += "-"} # Increase length with gaps. Assumes all indels properly aligned (may need local realignment for multiallelic sites)
					var[0] = variants[0][0] if var[0].chr == "." # Replace with character from reference if needed
				end
				variants[i] = var # Variable scope handling
			end
			current_base = line_arr[1].to_i # Set starting base position	
			if current_base > previous_endex
				for i in 0...$samples.size
					unless $options.skip # Adjust for missing bases
						current_locus.seqs[i] << "?" * (current_base - 1 - previous_endex)
						current_locus.alts[i] << "?" * (current_base - 1 - previous_endex)
					end
				end
				if write_cycle >= $options.write_cycle
					current_locus.write_seqs
					write_cycle = 0
				end
				index =  current_locus.seqs[0].size - 1
				previous_index = current_base - 1
				previous_endex = current_base - 1
			end
			for i in 0...$samples.size # Adjust locus lengths
				current_locus.seqs[i] << "?" * lengths.max * (current_base - previous_endex) if current_base > previous_endex
				current_locus.alts[i] << "?" * lengths.max * (current_base - previous_endex) if current_base > previous_endex
			end
			index += current_base - previous_index
			endex = index + lengths.max - 1 # Sequence end index
			unless $options.hap_flag #Accounting for ploidy
				for i in 9...line_arr.size
					vars = line_arr[i].split(":")[vars_index] # This code handles phasing and randomizes unphased diplotypes
					vars = line_arr[i].split(":")[gt_index] if (codes.include?("PGT") && vars[0].chr == ".") # Correction for homozygous reference phasing
					if $options.probps && !ad_index.nil? && vars != "./." && vars != ".|."
						ad_arr = line_arr[i].split(":")[ad_index].split(",").map { |x| x.to_i }
						if $options.adepth
							new_ad = []
							for ad in ad_arr
								ad = 0 if ad < $options.adepth
								new_ad.push(ad)
							end
							ad_arr = new_ad
						end
						ad_arr_sum = ad_arr.reduce(:+) # Code here for backwards compatibility since macOS only has Ruby 2.0
						randvar = rand(ad_arr_sum)
						ad_sum = 0
						for ad in ad_arr
							ad_sum += ad
							if randvar <= ad_sum - 1
								vars[0] = vars[2] = ad_arr.index(ad).to_s
								break
							end
						end
						randvar = 0
					else
						randvar = rand(2)
					end
					if !$options.ambig or endex - index > 0 or variants.include?("-") # Phase haplotypes even in ambiguity situations
						if vars[0].chr != "." # Code below uses string replacement due to multiallelic sites having same start index
							if vars[1].chr == "|" or randvar == 0
								current_locus.seqs[i-9][index..endex] = variants[vars[0].to_i]
							else
								current_locus.alts[i-9][index..endex] = variants[vars[0].to_i]
							end
						else
							if vars[1].chr == "|" or randvar == 0
								current_locus.seqs[i-9][index..endex] = "?" * (endex - index + 1)
							else
								current_locus.alts[i-9][index..endex] = "?" * (endex - index + 1)
							end
						end
						if vars[2].chr != "."
							if vars[1].chr == "|" or randvar == 0
								current_locus.alts[i-9][index..endex] = variants[vars[2].to_i]
							else
								current_locus.seqs[i-9][index..endex] = variants[vars[2].to_i]
							end
						else
							if vars[1].chr == "|" or randvar == 1
								current_locus.seqs[i-9][index..endex] = "?" * (endex - index + 1)
							else
								current_locus.alts[i-9][index..endex] = "?" * (endex - index + 1)
							end
						end
					elsif $options.ambig
						if vars[0] == "." and vars[2] == "."
							code = "?"
						elsif vars[0] == "."
							code = variants[vars[2].to_i]
						elsif vars[1] == "."
							code = variants[vars[0].to_i]
						else
							code = get_ambiguity_code(variants[vars[0].to_i].upcase, variants[vars[2].to_i].upcase)
						end
						current_locus.seqs[i-9][index..endex] = code
						current_locus.alts[i-9][index..endex] = code
					end
				end
			else
				for i in 9...line_arr.size
					vars = line_arr[i].split(":")[vars_index]
					if vars[0] == "." # Some tools (e.g. VCFtools) output diploid missing data calls even for haploid VCFs
						current_locus.seqs[i-9][index..endex] = "?" * (endex - index + 1)
					else 
						current_locus.seqs[i-9][index..endex] = variants[vars.to_i]
					end
				end
			end
			previous_index = current_base
			previous_endex = current_base + lengths.max - 1  if current_base + lengths.max - 1 > previous_endex
			prev_pos = current_base
		end
	end
	return index, previous_index, previous_endex, previous_name, current_locus, prev_pos, write_cycle, region, regionval
end
#-----------------------------------------------------------------------------------------------
def read_input
	index = -1 # Set internal start index
	previous_index = -1 # Index of previous ending base
	previous_endex = 0 # End index of indels
	previous_name = "" # Name comparator for concatenated alignments
	current_locus = Locus.new("") # Current locus
	prev_pos = 0
	write_cycle = 0
	region = 1
	regionval = 0
	if $options.pipe
		ARGF.each_line { |line| index, previous_index, previous_endex, previous_name, current_locus, prev_pos, write_cycle, region, regionval = vcf_to_alignment(line, index, previous_index, previous_endex, previous_name, current_locus, prev_pos, write_cycle, region, regionval) }
	elsif $options.infile[-3..-1] == ".gz"
		Zlib::GzipReader.open($options.infile) do |vcf2aln|
	 		while line = vcf2aln.gets
	 			index, previous_index, previous_endex, previous_name, current_locus, prev_pos, write_cycle, region, regionval = vcf_to_alignment(line, index, previous_index, previous_endex, previous_name, current_locus, prev_pos, write_cycle, region, regionval)
	 		end
		end
	else
		File.open($options.infile, 'r') do |vcf2aln|
	 		while line = vcf2aln.gets
	 			index, previous_index, previous_endex, previous_name, current_locus, prev_pos, write_cycle, region, regionval = vcf_to_alignment(line, index, previous_index, previous_endex, previous_name, current_locus, prev_pos, write_cycle, region, regionval)
	 		end
		end
	end
	current_locus.print_locus # Print final alignment
	current_locus.write_partitions if $options.partition
end
#-----------------------------------------------------------------------------------------------
ARGV[0] ||= "-h"
$options = Parser.parse(ARGV)
if $options.pipe
	ARGV.clear # Clear the input array for piping to ARGF
else
	while !FileTest.exist?($options.infile)
		print "Input file not found. Please re-enter.\n"
		$options.infile = gets.chomp
	end
end
build_ambig_hash if $options.ambig
get_GT_fields($options.infile) if $options.type_fields
read_input
