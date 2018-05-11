#!/usr/bin/env ruby

#-----------------------------------------------------------------------------------------------
# vcf2aln
# Michael G. Campana, 2017
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------------

require 'optparse'
require 'ostruct'

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
				File.open(sample + "_#{@name}.hap1.tmp.fa", 'w') do |write|
					write.puts ">" + sample + "_hap1"
				end
				File.open(sample + "_#{@name}.hap2.tmp.fa", 'w') do |write|
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
			File.open($samples[i] + "_#{@name}.hap1.tmp.fa", 'a') do |write|
				write << @seqs[i]
			end
			File.open($samples[i] + "_#{@name}.hap2.tmp.fa", 'a') do |write|
				write << @alts[i]
			end
			if $options.maxmissing < 100.0
				@missinghap1[i] += @seqs[i].scan("?").length
				@missinghap2[i] += @seqs[i].scan("?").length
			end
			@seqs[i] = ""
			@alts[i] = ""
		end
	end	
	def print_locus
		@length += @seqs[0].length
		for i in 0...$samples.size
			File.open($samples[i] + "_#{@name}.hap1.tmp.fa", 'a') do |write|
				write.puts @seqs[i]
			end
			File.open($samples[i] + "_#{@name}.hap2.tmp.fa", 'a') do |write|
				write.puts @alts[i]
			end
			if $options.maxmissing < 100.0
				@missinghap1[i] += @seqs[i].scan("?").length
				@missinghap2[i] += @seqs[i].scan("?").length
				system("rm #{$samples[i] + '_' + @name}.hap1.tmp.fa") if @missinghap1[i].to_f/@length.to_f * 100.0 > $options.maxmissing
				system("rm #{$samples[i] + '_' + @name}.hap2.tmp.fa") if @missinghap2[i].to_f/@length.to_f * 100.0 > $options.maxmissing
			end
		end
		if $options.alts
			system("cat *#{@name}*.tmp.fa > #{@name + '.fa'}")
		else
			system("cat *#{@name}.hap1.tmp.fa > #{@name + '.hap1.fa'}")
			system("cat *#{@name}.hap2.tmp.fa > #{@name + '.hap2.fa'}")
		end
		system("rm *#{@name}*.tmp.fa")
	end
end
#-----------------------------------------------------------------------------------------------
class Parser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		args.infile = "" # Primary input file
		args.concat = false # Concatenate markers into single alignment
		args.skip = false # Skip missing sites in vcf
		args.mincalls = 0 # Minimum number of calls to include site
		args.maxmissing = 100.0 # Maximum percent missing data to include sequence
		args.alts = false # Print alternate haplotypes in same file
		opt_parser = OptionParser.new do |opts|
			opts.banner = "Command-line usage: ruby vcf2aln.rb [options]"
			opts.on("-i","--input [FILE]", String, "Input VCF file") do |vcf|
				args.infile = vcf
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
			opts.on_tail("-h","--help", "Show help") do
				puts opts
			end
		end
		opt_parser.parse!(options)
		return args
	end
end
#-----------------------------------------------------------------------------------------------
ARGV[0] = "-h" if ARGV[0].nil?
$options = Parser.parse(ARGV)
current_locus = Locus.new("") # Current locus
while !FileTest.exist?($options.infile)
	print "Input file not found. Please re-enter.\n"
	$options.infile = gets.chomp
end
index = -1 # Set internal start index
previous_index = -1 # Index of previous ending base
previous_endex = 0 # End index of indels
previous_name = "" # Name comparator for concatenated alignments
File.open($options.infile, 'r') do |vcf2aln|
	while line = vcf2aln.gets
		if line[0..1] == "#C"
			$samples = line[0..-2].split("\t")[9..-1] # Get sample names
		elsif line[0].chr != "#"
			line_arr = line.split("\t")
			$options.concat ? name = "concat_aln" : name = line_arr[0]
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
			end
			if line_arr[0] != previous_name # Reset indexes for concatenated alignments
				current_locus.write_seqs
				index = -1 #
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
				index += current_base - previous_index
				endex = index + lengths.max - 1 # Sequence end index
				for i in 9...line_arr.size
					vars = line_arr[i].split(":")[0] # This code handles phasing and randomizes unphased diplotypes
					randvar = rand(2)
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
				end
				previous_index = current_base
				previous_endex = current_base + lengths.max - 1  if current_base + lengths.max - 1 > previous_endex	
			end
		end
	end
	current_locus.print_locus # Print final alignment
end
