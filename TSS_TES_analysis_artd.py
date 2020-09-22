
"""Script for comparing enrichment of canonical promotor and polyA motifs around TSS and TES sites in arabidopsis Iso-Seq, Morton et al. TSS sites and Simpson et al. polyA sites"""



import re
import random
import matplotlib
import pandas as pd
import numpy
import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from plotnine import *



arabidopsis_genome = "/mnt/shared/scratch/mc42302/201903_RTD2/arabidopsis/ArabidopsisGenome.fa"

arabidopsis_TSS_original = "/mnt/shared/scratch/mc42302/201903_RTD2/arabidopsis/tpc_AllPeak/tpc_AllPeak.txt"
arabidopsis_TSS_IsoSeq = "/mnt/shared/scratch/mc42302/201903_RTD2/arabidopsis/TSS_enriched/TSS_enriched.txt"
arabidopsis_TSS_original_hc = "/mnt/shared/scratch/mc42302/201903_RTD2/arabidopsis/tpc_AllPeak/tpc_TSS_allPeaks2.txt"

arabidopsis_IsoSeq_TES = "/mnt/shared/scratch/mc42302/201903_RTD2/arabidopsis/TES_enriched/TES_enriched.txt"
gordon_TES = "/mnt/shared/scratch/mc42302/201903_RTD2/arabidopsis/TES_DRS/TES_DRS.txt"

outpath = "/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/TSS_TES_validation/"

def generate_motif_plot_arabidopsis(motif,arabidopsis_TSS_original,chrom_index,strand_index,TSS_index,ara_chromosomes,filename,colour):
	histo = [0] * 1101
	histokeys = list(range(-550,551))
	for n,line in enumerate(open(arabidopsis_TSS_original)):
		if not n:
			continue
		#print(line)
		line1 = line.rstrip("\n").split("\t")
		chromosome = line1[chrom_index]
		strand = line1[strand_index]
		#print(strand)
		TSS_position = int(line1[TSS_index])
		if chromosome not in ara_chromosomes.keys():
			chromosome = "C" + chromosome[1:] #Gordon's data has lower case c, this needs to be changed
		if strand == "+" or strand == "fwd":
			sequence = ara_chromosomes[chromosome][TSS_position - 550:TSS_position + 551].upper()
		elif strand == "-" or strand == "rev":
			sequence_reverse = Seq(ara_chromosomes[chromosome][TSS_position - 550:TSS_position + 551].upper(),IUPAC.unambiguous_dna)
			sequence = str(sequence_reverse.reverse_complement())
		else:
			print("Input error")
		#now find if motif is in region +/- 500
		motif_positions = motif_pos(motif,sequence)
		if motif_positions:
			#assert len(motif_position) == 2
			for motif_position in motif_positions:
				for i in range(motif_position[0],motif_position[1]):
					histo[i] += 1
	#####
	histo_dict = dict(zip(histokeys,histo))
	histo_dict2 = {a:z for a,z in histo_dict.items() if a >= -500 and a <= 500}
	#NOw do control run
	histo_control = [0] * 1101
	strands = ["+","-"]
	for i in range(0,len(open(arabidopsis_TSS_original).readlines())):
		#First gent random chromosome
		random.shuffle(strands)
		strand = strands[0]
		chromosomes = list(ara_chromosomes.keys())
		random.shuffle(chromosomes)
		chromosome = chromosomes[0]
		TSS_position = random.randint(550,len(ara_chromosomes[chromosome]) - 551)
		if strand == "+":
			sequence = ara_chromosomes[chromosome][TSS_position - 550:TSS_position + 551].upper()
		elif strand == "-":
			sequence_reverse = Seq(ara_chromosomes[chromosome][TSS_position - 550:TSS_position + 551].upper(),IUPAC.unambiguous_dna)
			sequence = str(sequence_reverse.reverse_complement())
		else:
			print("Input error")
		motif_positions = motif_pos(motif,sequence)
		if motif_positions:
			#assert len(motif_position) == 2
			for motif_position in motif_positions:
				for i in range(motif_position[0],motif_position[1]):
					histo_control[i] += 1
	histo_dict_control = dict(zip(histokeys,histo_control))
	histo_dict_control2 = {a:z for a,z in histo_dict_control.items() if a >= -500 and a <= 500}
	print(str(histo_dict_control2))
	line_plot(histo_dict2,histo_dict_control2,filename,"position in nt from TSS","Instances",colour,False)
	#Get frequency of motif
	length = len(open(arabidopsis_TSS_original).readlines())
	histo_dict_frequency = {a:z/length for a,z in histo_dict2.items()}
	histo_dict_control_frequency = {a:z/length for a,z in histo_dict_control2.items()}
	line_plot(histo_dict_frequency,histo_dict_control_frequency,filename + "_frequency","position in nt from TSS","Instances",colour,False)

def generate_motif_plot_bothcompare(motif,enriched_TSS,enriched_TSS2,chrom_index,strand_index,TSS_index,chrom_index2,strand_index2,TSS_index2,genome1,genome2,filename,colour,colour2,data_type):
	"""Create plot with enriched sites, not enriched sites and random sites"""
	histo = [0] * 1101
	histokeys = list(range(-550,551))
	for n,line in enumerate(open(enriched_TSS)):
		if not n:
			continue
		line1 = line.rstrip("\n").split("\t")
		chromosome = line1[chrom_index]
		strand = line1[strand_index]
		TSS_position = int(line1[TSS_index])
		if chromosome not in ara_chromosomes.keys():
			chromosome = "C" + chromosome[1:] #Gordon's data has lower case c, this needs to be changed
		if strand == "+" or strand == "fwd":
			sequence = genome1[chromosome][TSS_position - 550:TSS_position + 551].upper()
		elif strand == "-" or strand == "rev":
			sequence_reverse = Seq(genome1[chromosome][TSS_position - 550:TSS_position + 551].upper(),IUPAC.unambiguous_dna)
			sequence = str(sequence_reverse.reverse_complement())
		else:
			print("Input error")
		#now find if motif is in region +/- 500
		motif_positions = motif_pos(motif,sequence)
		if motif_positions:
			#assert len(motif_position) == 2
			for motif_position in motif_positions:
				for i in range(motif_position[0],motif_position[1]):
					histo[i] += 1
	histo_dict = dict(zip(histokeys,histo))
	histo_dict2 = {a:z for a,z in histo_dict.items() if a >= -500 and a <= 500}
	########
	histo_no_enriched = [0] * 1101
	for n,line in enumerate(open(enriched_TSS2)):#gene, chromosome,TSS,strand
		if not n:
			continue
		line1 = line.rstrip("\n").split("\t")
		chromosome = line1[chrom_index2]
		strand = line1[strand_index2]
		TSS_position = int(line1[TSS_index2])
		if chromosome not in ara_chromosomes.keys():
			chromosome = "C" + chromosome[1:] #Gordon's data has lower case c, this needs to be changed
		if strand == "+" or strand == "fwd":
			sequence = genome2[chromosome][TSS_position - 550:TSS_position + 551].upper()
		elif strand == "-" or strand == "rev":
			sequence_reverse = Seq(genome2[chromosome][TSS_position - 550:TSS_position + 551].upper(),IUPAC.unambiguous_dna)
			sequence = str(sequence_reverse.reverse_complement())
		else:
			print("Input error")
		#now find if motif is in region +/- 500
		motif_positions = motif_pos(motif,sequence)
		if motif_positions:
			#assert len(motif_position) == 2
			for motif_position in motif_positions:
				for i in range(motif_position[0],motif_position[1]):
					histo_no_enriched[i] += 1
	histo_dict_notenriched_full = dict(zip(histokeys,histo_no_enriched))
	histo_dict_notenriched = {a:z for a,z in histo_dict_notenriched_full.items() if a >= -500 and a <= 500}
	#####
	
	#NOw do control run
	histo_control = [0] * 1101
	strands = ["+","-"]
	for i in range(0,max([len(open(enriched_TSS).readlines()),len(open(enriched_TSS2).readlines())])):
		#First gent random chromosome
		random.shuffle(strands)
		strand = strands[0]
		chromosomes = list(genome1.keys())
		random.shuffle(chromosomes)
		chromosome = chromosomes[0]
		TSS_position = random.randint(550,len(genome1[chromosome]) - 551)
		if strand == "+":
			sequence = genome1[chromosome][TSS_position - 550:TSS_position + 551].upper()
		elif strand == "-":
			sequence_reverse = Seq(genome1[chromosome][TSS_position - 550:TSS_position + 551].upper(),IUPAC.unambiguous_dna)
			sequence = str(sequence_reverse.reverse_complement())
		else:
			print("Input error")
		motif_positions = motif_pos(motif,sequence)
		if motif_positions:
			#assert len(motif_position) == 2
			for motif_position in motif_positions:
				for i in range(motif_position[0],motif_position[1]):
					histo_control[i] += 1
	histo_dict_control = dict(zip(histokeys,histo_control))
	histo_dict_control2 = {a:z for a,z in histo_dict_control.items() if a >= -500 and a <= 500}

	line_plot_2(histo_dict2,histo_dict_notenriched,histo_dict_control2,filename,f"position relative to predicted {data_type} (nt)","Instances",colour,colour2,False,len(open(enriched_TSS).readlines()),len(open(enriched_TSS2).readlines()))
	line_plot_2(histo_dict2,histo_dict_notenriched,histo_dict_control2,filename + "_frequency",f"position relative to predicted {data_type} (nt)",f"Motif instances per {data_type}",colour,colour2,True,len(open(enriched_TSS).readlines()),len(open(enriched_TSS2).readlines()))

def line_plot_2(input_dictionary1,input_dictionary2,input_dictionary3,filename,x_axis_title,y_axis_title,bar_colour,bar_colour2,is_frequency,total1,total2):
	histokeys = list(input_dictionary1.keys())
	histovalues1 = list(input_dictionary1.values())
	histovalues2 = list(input_dictionary2.values())
	histovalues3 = list(input_dictionary3.values())
	if is_frequency:
		histovalues1 = [i/total1 for i in histovalues1]
		histovalues2 = [i/total2 for i in histovalues2]
		histovalues3 = [i/max([total1,total2]) for i in histovalues3]
	#histokeys,histovalues=histokeys[1:],histovalues[1:]
	pandadict = {x_axis_title:histokeys,y_axis_title:histovalues1,"second":histovalues2,"random":histovalues3}
	histodata = pd.DataFrame(pandadict)
	t = ggplot(aes(x = x_axis_title), data = histodata) + geom_line(aes(y = y_axis_title), color = bar_colour) + geom_line(aes(y = "second"), color = bar_colour2) + geom_line(aes(y = "random"), color = "grey")
	p = t + theme_bw() + labs(title="", x = x_axis_title,y = y_axis_title)
	p.save(filename)

def motif_pos(motif,sequence):
	p = re.compile(motif)
	if not p.search(sequence):
		return []
	else:
		coordinates = []
		iterator = p.finditer(sequence)
		for match in iterator:
			coordinates.append(match.span())
		return coordinates

def parse_genome(genome_input):
	genome = open(genome_input).read()
	#Split genome into chromosome dictionary with key = chr id and item chromosome sequence
	chromosomes = genome.split(">")
	del chromosomes[0]
	chromosomedict = {}
	for chromosome in chromosomes:
		splitchrom = chromosome.split("\n")
		chromid = splitchrom[0]
		chromsequence = "".join(splitchrom[1:])
		chromosomedict[chromid] = chromsequence #Genome sequence now in dictionary with chromosomes as key
	return chromosomedict


###Compare both Morton and Iso-Seq motifs

ara_chromosomes = parse_genome(arabidopsis_genome)

tatamotif = "TATA[AT]A[AT]"
generate_motif_plot_bothcompare(tatamotif,arabidopsis_TSS_IsoSeq,arabidopsis_TSS_original,1,3,2,0,1,5,ara_chromosomes,ara_chromosomes,outpath + "TATA TSS arabidopsis isoseq Morton","blue","red","TSS")

Inr = "[CT]TCA[ATCG]T[CT][CT]"
generate_motif_plot_bothcompare(Inr,arabidopsis_TSS_IsoSeq,arabidopsis_TSS_original,1,3,2,0,1,5,ara_chromosomes,ara_chromosomes,outpath + "Inr TSS arabidopsis isoseq Morton","blue","red","TSS")

ypatch4 = "[CT][CT][CT][CT][CT]C"
generate_motif_plot_bothcompare(ypatch4,arabidopsis_TSS_IsoSeq,arabidopsis_TSS_original,1,3,2,0,1,5,ara_chromosomes,ara_chromosomes,outpath + "ypatch4 TSS arabidopsis isoseq Morton","blue","red","TSS")

kozakmotif5 = "[ATCG][AC][AG][AC][ATCG]ATGGCG"
generate_motif_plot_bothcompare(kozakmotif5,arabidopsis_TSS_IsoSeq,arabidopsis_TSS_original,1,3,2,0,1,5,ara_chromosomes,ara_chromosomes,outpath + "Kozak motif TSS arabidopsis isoseq Morton","blue","red","TSS")



#Now do with Morton High confidence:
tatamotif = "TATA[AT]A[AT]"
generate_motif_plot_bothcompare(tatamotif,arabidopsis_TSS_IsoSeq,arabidopsis_TSS_original_hc,1,3,2,0,1,5,ara_chromosomes,ara_chromosomes,outpath + "TATA TSS arabidopsis isoseq Morton HC","blue","red","TSS")

Inr = "[CT]TCA[ATCG]T[CT][CT]"
generate_motif_plot_bothcompare(Inr,arabidopsis_TSS_IsoSeq,arabidopsis_TSS_original_hc,1,3,2,0,1,5,ara_chromosomes,ara_chromosomes,outpath + "Inr TSS arabidopsis isoseq Morton HC","blue","red","TSS")


ypatch4 = "[CT][CT][CT][CT][CT]C"
generate_motif_plot_bothcompare(ypatch4,arabidopsis_TSS_IsoSeq,arabidopsis_TSS_original_hc,1,3,2,0,1,5,ara_chromosomes,ara_chromosomes,outpath + "ypatch4 TSS arabidopsis isoseq Morton HC","blue","red","TSS")

kozakmotif5 = "[ATCG][AC][AG][AC][ATCG]ATGGCG"
generate_motif_plot_bothcompare(kozakmotif5,arabidopsis_TSS_IsoSeq,arabidopsis_TSS_original_hc,1,3,2,0,1,5,ara_chromosomes,ara_chromosomes,outpath + "Kozak motif TSS arabidopsis isoseq Morton HC","blue","red","TSS")






#TES, comparing Iso-Seq with Gordon's data
#Get Gordons data
CFlm = "TGTA"
generate_motif_plot_bothcompare(CFlm,arabidopsis_IsoSeq_TES,gordon_TES,1,3,2,1,2,3,ara_chromosomes,ara_chromosomes,outpath + "CFlm TES arabidopsis isoseq Simpson","blue","red","TES")

PAS = "AATAAA"
generate_motif_plot_bothcompare(PAS,arabidopsis_IsoSeq_TES,gordon_TES,1,3,2,1,2,3,ara_chromosomes,ara_chromosomes,outpath + "PAS TES arabidopsis isoseq Simpson","blue","red","TES")
