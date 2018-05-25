#!/usr/bin/env ruby

#-----------------------------------------------------------------------------------------------
# vcf2aln
VCF2ALNVER = "0.3.0"
# Michael G. Campana, Jacob A. West-Roberts, 2017-2018
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

require 'optparse'
require 'ostruct'

$categories = {"GT" => "GT: Genotype information", "DP" => "DP: Sample read depth", "FT" => "FT: Sample genotype filter indicating if this genotype call passed quality filters", "GL" => "GL: Genotype log likelihoods (base 10) showing likelihood of each possible phased genotype", "GLE" => "GLE: Genotype likelihoods for uncertain copy number/phasing (includes all possible arrangements)", "PL" => "PL: Phred-scaled genotype likelihoods rounded to the nearest integer", "GP" => "GP: Phred-scaled genotype posterior probabilities", "AD" => "AD: Allele Depth", "GQ" => "GQ: Conditional genotype quality (phred-encoded): probability of the wrong genotype call conditioned on this site being variant", "HQ" => "HQ: Haplotype qualities", "PS" => "PS: Phase set: the set of phased genotypes to which this genome belongs", "PQ" => "PQ: Phase Quality. Phred-scaled probability that alleles are ordered incorrectly in a heterozygote (Not commonly used)", "EC" => "Comma separated list of expected alternate allele counts", "MQ" => "MQ: RMS mapping quality (Integer)"}

class Locus
	attr_accessor :name, :seqs, :alts, :length
	def initialize(name = "", seqs = [], alts = [], length = 0)
		@name = name
		@seqs = seqs
		@alts = alts
		@length = length
		@missinghap1 = []
		@missinghap2 = []
		unless @name == ""
			for sample in $samples
				File.open(sample + "_#{$options.outprefix}#{@name}.hap1.tmp.fa", 'w') do |write|
					write.puts ">" + sample + "_hap1"
				end
				File.open(sample + "_#{$options.outprefix}#{@name}.hap2.tmp.fa", 'w') do |write|
					write.puts ">" + sample + "_hap2"
				end
				@missinghap1.push(0)
				@missinghap2.push(0)
			end
		end
	end
	def write_seqs
		@length += @seqs[0].length
		for i in 0...$samples.size
			File.open($samples[i] + "_#{$options.outprefix}#{@name}.hap1.tmp.fa", 'a') do |write|
				write << @seqs[i]
			end
			unless $options.hap_flag
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
	def print_locus(line_num)
		#print "print_locus called at line #{line_num}", "\n"
		@length += @seqs[0].length
		if $options.split_regions == 0
			for i in 0...$samples.size
				File.open($samples[i] + "_#{$options.outprefix}#{@name}.hap1.tmp.fa", 'a') do |write|
					write.puts @seqs[i]
				end
				unless $options.hap_flag
					File.open($samples[i] + "_#{$options.outprefix}#{@name}.hap2.tmp.fa", 'a') do |write|
						write.puts @alts[i]
					end
				end
				if $options.maxmissing < 100.0
					@missinghap1[i] += @seqs[i].scan("?").length
					@missinghap2[i] += @seqs[i].scan("?").length unless $options.hap_flag
					system("rm #{$samples[i] + '_' + $options.outprefix + @name}.hap1.tmp.fa") if @missinghap1[i].to_f/@length.to_f * 100.0 > $options.maxmissing
					system("rm #{$samples[i] + '_' + $options.outprefix + @name}.hap2.tmp.fa") if @missinghap2[i].to_f/@length.to_f * 100.0 > $options.maxmissing unless $options.hap_flag
				end
			end
			if $options.alts
				system("cat *#{$options.outprefix + @name}*.tmp.fa > #{$options.outprefix + @name + '.fa'}")
			else
				system("cat *#{$options.outprefix + @name}.hap1.tmp.fa > #{$options.outprefix + @name + '.hap1.fa'}")
				system("cat *#{$options.outprefix + @name}.hap2.tmp.fa > #{$options.outprefix + @name + '.hap2.fa'}") unless $options.hap_flag
			end
			system("rm *#{$options.outprefix + @name}*.tmp.fa")
		else
			for i in 0...$samples.size
				File.open($samples[i] + "_#{$options.outprefix}#{@name}.hap1.tmp.fa", 'a') do |write|
					write.puts @seqs[i]
				end
				unless $options.hap_flag
					File.open($samples[i] + "_#{$options.outprefix}#{@name}.hap2.tmp.fa", 'a') do |write|
						write.puts @alts[i]
					end
				end
				if $options.maxmissing < 100.0
					@missinghap1[i] += @seqs[i].scan("?").length
					@missinghap2[i] += @seqs[i].scan("?").length unless $options.hap_flag
					system("rm #{$samples[i] + '_' + $options.outprefix + @name}.hap1.tmp.fa") if @missinghap1[i].to_f/@length.to_f * 100.0 > $options.maxmissing
					system("rm #{$samples[i] + '_' + $options.outprefix + @name}.hap2.tmp.fa") if @missinghap2[i].to_f/@length.to_f * 100.0 > $options.maxmissing unless $options.hap_flag
				end
			end
			if $options.alts
				system("cat *#{$options.outprefix + @name}*.tmp.fa > #{$options.outprefix + @name + '.fa'}")
			else
				system("cat *#{$options.outprefix + @name}.hap1.tmp.fa > #{$options.outprefix + @name + '_region' + $num_regions.to_s + '.hap1.fa'}")
				system("cat *#{$options.outprefix + @name}.hap2.tmp.fa > #{$options.outprefix + @name + '_region' + $num_regions.to_s + '.hap2.fa'}") unless $options.hap_flag
			end
			system("rm *#{$options.outprefix + @name}*.tmp.fa")
			$num_regions += 1
		end
	end
end
#-----------------------------------------------------------------------------------------------
class Parser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		args.infile = "" # Primary input file
		args.outprefix = "" # Output prefix
		args.hap_flag = false
		args.concat = false # Concatenate markers into single alignment
		args.skip = false # Skip missing sites in vcf
		args.mincalls = 0 # Minimum number of calls to include site
		args.maxmissing = 100.0 # Maximum percent missing data to include sequence
		args.alts = false # Print alternate haplotypes in same file
		args.ambig = false # Print SNPs as ambiguity codes
		args.qual_filter = 0 #Minimum quality for site (QUAL column)
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
		opt_parser = OptionParser.new do |opts|
			opts.banner = "Command-line usage: ruby vcf2aln.rb [options]"
			opts.on("-i","--input [FILE]", String, "Input VCF file") do |vcf|
				args.infile = vcf
			end
			opts.on("-o", "--outprefix [VALUE]", String, "Output alignment prefix") do |pref|
				args.outprefix = pref + "_" if pref != nil
			end
			opts.on("-c", "--concatenate", "Concatenate markers into single alignment") do
				args.concat = true
			end
			opts.on("-s", "--skip", "Skip missing sites in vcf") do
				args.skip = true
			end
			opts.on("-m","--mincalls [VALUE]", Integer, "Minimum number of calls to include site (Default = 0)") do |msnps|
				args.mincalls = msnps if msnps != nil
			end
			opts.on("-x","--maxmissing [VALUE]", Float, "Maximum percent missing data to include sequence (Default = 100.0)") do |missing|
				args.maxmissing = missing if missing != nil
			end
			opts.on("-a", "--alts", "Print alternate haplotypes in same file") do
				args.alts = true
			end
			opts.on("-b", "--ambig", "Print SNP sites as ambiguity codes.") do
				args.ambig = true
			end
			opts.on("-N", "--hap_flag", "Flag for haplotype data") do 
				args.hap_flag = true
			end
			opts.on("-g", "--split_regions [VALUE]", Integer, "Split alignment into subregional alignments for phylogenetic analysis") do |regions|
				args.split_regions = regions if regions != nil
			end
			opts.on("-t", "--typefields", "Display VCF genotype field information, then quit the program.") do
				args.type_fields = true
			end
			opts.on("-q", "--qual_filter [VALUE]", Integer, "Minimum accepted value for QUAL (per site) (Default = 0)") do |qual|
				args.qual_filter = qual if qual != nil
			end
			opts.on("-y", "--site_depth [VALUE]", Integer, "Minimum desired depth for each site (DP, in INFO) (Default = No filter)") do |site|
				args.site_depth = site if site != nil
			end
			opts.on("-d", "--sampledepth [VALUE]", Integer, "Minimum allowed read depth (Default = No filter)") do |depth|
				args.sample_depth = depth if depth != nil
			end
			opts.on("-l", "--likelihood [VALUE]", Float, "Minimum allowed genotype log-likelihood (At least one option must satisfy this value)") do |likelihood|
				args.min_ll = likelihood if likelihood != nil
			end
			opts.on("-p", "--phred [VALUE]", Integer, "Minimum accepted phred-scaled genotype likelihood (Default = No filter)") do |phreddy|
				args.phred_likelihood = phreddy if phreddy != nil
			end
			opts.on("-P", "--posterior [VALUE]", Float, "Minimum accepted phred-scaled genotype posterior probability (Default = No filter)") do |post|
				args.posterior = post if post != nil
			end
			opts.on("-C", "--conditional [VALUE]", Float, "Minimum conditional genotype quality (phred-encoded) (Default = No filter)") do |condi|
				args.conditional = condi if condi != nil
			end
			opts.on("-H", "--haplotype_quality [VALUE]", Integer, "Minimum allowed haplotype quality (phred-encoded) (Default = No filter)") do |haplo|
				args.hap_qual = haplo if haplo != nil
			end
			opts.on("-r", "--sample_mq [VALUE]", Integer, "Minimum allowed per-sample RMS mapping quality (Default = No filter)") do |map_quality|
				args.sample_mq = map_quality if map_quality != nil
			end
			opts.on("-R", "--site_mq [VALUE]", Integer, "Minimum allowed per-site mapping quality (MQ in INFO) (Default = No filter)") do |mq|
				args.site_mq = mq if mq != nil
			end
			opts.on("-F", "--mq0f [VALUE]", Float, "Maximum allowed value for MQ0F. Must be between 0 and 1. (Default = No filter)") do |mqf|
				args.mq0f = mqf if mqf != nil
			end
			opts.on("-S", "--mqsb [VALUE]", Float, "Minimum allowed value for MQSB. (Default = No filter)") do |sb|
				args.mqsb = sb if sb != nil
			end
			opts.on("-A", "--adepth [VALUE]", Integer, "Minimum allowed allele depth. (Default = No filter)") do |ad|
				args.adepth = ad if ad != nil
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
def quality_filter(line_arr)
	site_qual = line_arr[5]
	filter = line_arr[6]
	site_info_fields = line_arr[7]
	samples = line_arr[9..-1]
	info_arr = site_info_fields.split(";")
	samples.map!{ |element| element.split(":")}
	genotypes = []
	samples.each{|list| genotypes.push(list[0])}
	# Leave if nothing to replace
	return line_arr if genotypes.all? {|x| x == "./."}
	found = false
	found = true if site_qual.to_i < $options.qual_filter || (filter != 'PASS' && filter != ".")
	unless ($options.mq0f.nil? && $options.mqsb.nil? && $options.sample_mq.nil? && $options.site_depth.nil?) || found #Don't execute this loop if you're not trying to filter for any site-specific quality scores
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
	unless found || ($options.sample_mq.nil? && $options.sample_depth.nil? && $options.min_ll.nil? && $options.phred_likelihood.nil? && $options.posterior.nil? && $options.conditional.nil? && $options.hap_qual.nil? && $options.adepth.nil?)
		sample_info_fields = line_arr[8]
		sample_info_fields = sample_info_fields.split(":")
		index_hash = {}
		#Create index hash; info fields may be ordered in any way, so we need a generalizable method to access each field at its proper position in the sample info
		sample_info_fields.each{|v| index_hash[v] = sample_info_fields.index(v)}
		#This bit can be sped up by removing options from the index hash if the corresponding value in $options does not exist. 
		for sample in samples do
			break if found
			for field in sample_info_fields do
				break if found
				case field
				when "GT"
					next
				when "DP"
					if $options.sample_depth && sample[index_hash["DP"]].to_i < $options.sample_depth
					 	found = true
					 	break
				 	end
			 	when "GL"
			  		if $options.min_ll
					 	gl_array = sample[index_hash["GL"]].split(",")
						for element in gl_array do
							if element.to_f < $options.min_ll
								found = true
								break
							end
						end
					end
				when "GLE"
					puts "** GLE not supported as of version #{VCF2ALNVER} **"
					puts "** Treating line: #{line_arr.join("\t")} as missing data **"
					found = true
					break
				when "PQ"
					next
				when "PL"
					if $options.phred_likelihood && sample[index_hash["PL"]].to_i < $options.phred_likelihood
						found = true
						break
					end
				when "GP"
					gp_array = sample[index_hash["GP"]].split(",")
					if $options.posterior
						for element in gp_array do
							if element.to_f < $options.posterior
								found = true
								break
							end
						end
					end
				when "GQ"
					if $options.conditional && sample[index_hash["GQ"]].to_i < $options.conditional
						found = true
						break
					end
				when "HQ"
					hap_array = sample[index_hash["HQ"]].split(",")
					if $options.hap_qual
						for element in hap_array do
							if element.to_f < $options.hap_qual
								found = true
								break
							end
						end
					end
				when "MQ"
					if $options.sample_mq && sample[index_hash["MQ"]].to_i < $options.sample_mq
						found = true
						break
					end
				when "AD"
					if $options.adepth && sample[index_hash["AD"]].to_i < $options.adepth
						found = true
						break
					end
				end
			end
		end
	end
	# Replace filtered genotypes
	if found
		genotypes.each { |genotype| genotype.replace("./.") if genotype != "./." }
		new_samples = []
		samples.each{|list|
			list.delete_at(0)
			list.unshift("./.")
			new_list = list.join(":")
			new_samples.push(new_list)
		}
		line_arr[9..-1] = new_samples
		return line_arr
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
	@fields = []
	File.open(vcf_file, 'r') do |vcf2aln|
		while line = vcf2aln.gets
			if line[0].chr != "#"
				line_arr = line.split("\t")
				gt_fields = line_arr[8].split(":")
				@fields.push(gt_fields).flatten!.uniq!
			end
		end
	end
	# Moved this section since fields are not constant across all sites
	puts "Genotype field information categories:"
	@fields.each { |item|
		if $categories.has_key?(item)
			puts $categories[item]
			puts "** GLE containing vcf files not supported as of version #{VCF2ALNVER} **" if item.to_s == "GLE"
		else
			puts "** Genotype code #{item} not specified in VCF manual version 4.2 **"
		end
	}
	exit
end
#-----------------------------------------------------------------------------------------------
def vcf_to_alignment
	index = -1 # Set internal start index
	previous_index = -1 # Index of previous ending base
	previous_endex = 0 # End index of indels
	previous_name = "" # Name comparator for concatenated alignments
	current_locus = Locus.new("") # Current locus
	$num_regions = 0
	prev_pos = 0
	line_num = 0
	File.open($options.infile, 'r') do |vcf2aln|
		while line = vcf2aln.gets
			line_num += 1
			if line[0..1] == "#C"
				$samples = line[0..-2].split("\t")[9..-1] # Get sample names
			elsif line[0].chr != "#"
				line_arr = line.split("\t")
				line_arr = quality_filter(line_arr) # Quality filter has to be called to check for GLE
				$options.concat ? name = "concat_aln" : name = line_arr[0]
				if name != current_locus.name
					current_locus.print_locus(line_num) unless current_locus.name == ""
					# puts "print_locus call from station A"
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
				end
				if line_arr[0] != previous_name # Reset indexes for concatenated alignments
					current_locus.write_seqs
					previous_index = -1
					previous_endex = 0
					previous_name = line_arr[0]
				end
				if $samples.size - line.scan("./.").length - line.scan(".|.").length >= $options.mincalls
					variants = [line_arr[3], line_arr[4].split(",")].flatten!
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
								current_locus.seqs[i] += "?" * (current_base - 1 - previous_endex)
								current_locus.alts[i] += "?" * (current_base - 1 - previous_endex)
							end
						end
						current_locus.write_seqs #Write sequence to end
						index = -1
						previous_index = current_base - 1
						previous_endex = current_base - 1
					end
					for i in 0...$samples.size # Adjust locus lengths
						current_locus.seqs[i] += "?" * lengths.max * (current_base - previous_endex) if current_base > previous_endex
						current_locus.alts[i] += "?" * lengths.max * (current_base - previous_endex) if current_base > previous_endex
					end
					#Changed from previous_index to previous_endex; I'm going to use that variable. (J)
					#(Also Jacob) DON'T DO THAT!!!!!!!!!!!!!!!!!
					index += current_base - previous_index
					endex = index + lengths.max - 1 # Sequence end index
					
					

					#if current_base.to_i >= 0
						#print "line_arr[3] ", line_arr[3], "\n"
						#print "Current base, ", current_base, "\n"
						#print "previous index ", previous_index, "\n"
						#print "previous endex ", previous_endex, "\n"
						#print "index, ", index, "\n"
						#print "endex ", endex, "\n"	
						#puts "current_locus.seqs[0].size: ", current_locus.seqs[0].size			
						#for i in 9...line_arr.size	
						#	print current_locus.seqs[i-9][index..endex], "\n"
						#end
						
						
					#end

					unless $options.hap_flag #Don't run this subroutine on haploid vcf
						for i in 9...line_arr.size
							vars = line_arr[i].split(":")[0] # This code handles phasing and randomizes unphased diplotypes
							randvar = rand(2)
							if !$options.ambig or endex - index > 0 # Phase haplotypes even in ambiguity situations
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
							vars = line_arr[i]
							if vars == "."
								current_locus.seqs[i-9][index..endex] = "?" * (endex - index + 1)
							else 
								current_locus.seqs[i-9][index..endex] = variants[vars.to_i]
							end
						end
					end
					
					if (current_base - 1 - prev_pos) >= $options.split_regions && previous_index != -1 && $options.split_regions != 0
						current_locus.print_locus(line_num)
					#	puts "print_locus call from station B"
					#	puts "current base: #{current_base}", "prev_pos: #{prev_pos}", "Difference: #{current_base - 1 - prev_pos}", "previous_index: #{previous_index}"

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
						
						
						current_locus.write_seqs
						previous_index = -1
						previous_endex = 0
						prev_pos = current_base
					else
						previous_index = current_base
						previous_endex = current_base + lengths.max - 1  if current_base + lengths.max - 1 > previous_endex
						prev_pos = current_base
					end
				end
			end
		end
		current_locus.print_locus(line_num) # Print final alignment
	#	puts "print_locus call from station C"
	end
end
#-----------------------------------------------------------------------------------------------
ARGV[0] ||= "-h"
$options = Parser.parse(ARGV)
while !FileTest.exist?($options.infile)
	print "Input file not found. Please re-enter.\n"
	$options.infile = gets.chomp
end
build_ambig_hash if $options.ambig
get_GT_fields($options.infile) if $options.type_fields
vcf_to_alignment

