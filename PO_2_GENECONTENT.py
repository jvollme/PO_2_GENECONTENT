#!/usr/bin/python
#created 21.06.2014 by John Vollmers
#from Bio.Phylo.TreeConstruction import * #for using the directly biopython-implemented Distance(UPGMA + NJ) and parsimony methods. Not done/necessary yet
import os, sys, argparse, time, multiprocessing, random, traceback
from Bio import AlignIO, SeqIO
from subprocess import call
from Bio.Phylo.Applications import RaxmlCommandline, PhymlCommandline
myparser=argparse.ArgumentParser(description="\n==PO_2_GENECONTENT.py v1.01 by John Vollmers==\nCreates discrete binary character matrices and phylogenetic trees based on the gene content (presence/absence of orthologues) in comparison organisms.\nThis script is supposed to be part of a pipeline consisting of:\n\
\tA.)Conversion of Genbank/Embl-Files to ANNOTATED(!) Fastas\n\t   using CDS_extractor.pl by Andreas Leimbach\n\
\tB.)Calculation of orthologs and paralogs using proteinortho5\n\t   (WITH the '-single' and '-self' arguments!)\n\
\tC.)The creation of discrete binary character marices\n\t   (and optionally phylogenetic trees) based on:\n\t\t-the fasta sequences of step A\n\
\t\t-the proteinortho5-results from step B\n", formatter_class=argparse.RawTextHelpFormatter)
myparser.add_argument("-po", "--proteinortho", action = "store", dest = "po_resultfile", help = "(String) file with proteinortho5 results", required = True)
myparser.add_argument("-s", "--silent", action = "store_true", dest = "no_verbose", help = "non-verbose mode")
myparser.add_argument("-o", "--out", action = "store", dest="out_file", default = "output", help = "Basename for result-files\nDefault='output'\nWill create a phylip-file of the aligned discrete binary character data by default")
myparser.add_argument("-ofas", "--out_fasta", action = "store_true", dest = "outfasta", default = False, help = "Create a file with the aligned discrete binary character data in Fasta format\n(in addition to the default phylip file)")
myparser.add_argument("-omat", "--out_matrix", action = "store_true", dest = "outmatrix", default = False, help = "Create a file with the aligned discrete character data in tabular format\n(in addition to the default phylip file)")
myparser.add_argument("-odiff", "--out_difference", action = "store", dest = "outdiff", choices = ["simple", "none", "braycurtis", "canberra", "cityblock", "jaccard", "euclidean"], default = "simple", help = "Create difference/similarity matrixes for distance matrix based phylogeny inference.\n\t- simple : Simple normalized distance matrix\n\t distance=Nr shared OGs/total Nr OGs in SMALLER organism\n\t- braycurtis : Bray-Curtis distance\n\t- canberra : Canberra distance\n\t- cityblock : Manhatten distance\n\t- euclidean : Euclidean distance\n\t- jaccard : Jaccard-Needham dissimilarity\nDefault: \"simple\"") #only simply algorithm implemented for now. Add more complex methods and bootstrapping options
myparser.add_argument("-t", "--threads", action = "store", dest = "nthreads", type = int, default = 4, help = "Number of cpus to use\nDefault=4 (or maximum number of available cores, if less than 4 cores available)")
myparser.add_argument("-sd", "--seed", action = "store", dest = "seed_nr", type = int, default = 0, help = "Integer to provide as seed for RAxML or PhyML\n0=seed generated randomly\nDefault=random seed")
myparser.add_argument("-mt", "--make_tree", action = "store", dest = "tree_method", choices = ["raxml", "raxml_bs", "raxml_rapidbs", "none"], default = "none", help = "Generate ML phylogenetic trees using RAxML with \"new rapid hill climbing\"\nand substitution model: \"BINGAMMA\"\n\t-\"raxml\": single tree without bootstraps (using new rapid hill climbing)\n\t-\"raxml_bs\": thorough bootstrap analyses and search for best ML tree\n\t-\"raxml_rapidbs\": rapid bootstrap analyses and search for best ML tree\n\t in one run\n\t-\"none\"\nDefault=none")
myparser.add_argument("-tbp","--tree_builder_path", action="store", dest = "treebuilder_path", default = "", help = "Path to treebuilder (currently only raxml supported) if not listed in $PATH")
#myparser.add_argument("--raxml", action="store_true", dest="raxml", default=False, help="Generate ML phylogenetic tree without bootstrap values using RAxML and \"new rapid hill climbing\" and substitution model: \"BINGAMMA\"\nDefault=don't use")
#myparser.add_argument("--raxml_bs", action="store_true", dest="raxml_BS", default=False, help="Generate ML phylogenetic tree with thorough bootstrap analyses and search for best ML tree using RAxML and \"new rapid hill climbing\" and substitution model: \"BINGAMMA\"'\nDefault: don't use")
#myparser.add_argument("--raxml_rapidbs", action="store_true", dest="raxml_rapidBS", default=False, help="Generate ML phylogenetic tree with rapid bootstrap analyses and search for best ML tree in one run using RAxML and substitution model: \"BINGAMMA\"\nDefault: don't use")
#myparser.add_argument("--pars", action="store_true", dest="fpars", default=False, help="Generate Parsimony phylogenetic tree using fpars (the EMBOSS implementation of the PHYLIP tool)\nDefault: don't use")
#myparser.add_argument("--mix", action="store_true", dest="fmix", default=False, help="Generate Parsimony phylogenetic tree using fmix with wagner parsimony (the EMBOSS implementation of the PHYLIP tool)\nDefault=don't use")
#myparser.add_argument("--dollop", action="store_true", dest="fdollop", default=False, help="Generate Parsimony phyogenetic trees using fdollop and Dollo Parsimony (the EMBOSS implementation of the PHYLIP tool)\n Default=don't use")
myparser.add_argument("-bs", "--bootstraps", action = "store", dest = "bootstraps", type = int, default = 1000, help = "Number of times to resample for bootstrapping\nonly makes sense in combination with '-mt raxml_bs' or '-mt raxml_rapidbs'\nDefault=1000")
myparser.add_argument("--min_freq", action = "store", dest = "min_freq", type=int, default = 2, help = "Minimum frequency for Orthologeous Groups across all comparison-organisms for Group to be considered in Analyses\n Default=2 (exclude all singletons)")
#myparser.add_argument("--ignore_columns", action="store", dest="ignore_column_list", default=None, help="indexes (starting with 1) of Organisms to be ignored when counting frequencies of Orthologeous Groups (Tip: ignore all but one representative member of overrepresented species/strains)")
#implement bootsrapped versions of pars, mix and dollop soon
args = myparser.parse_args()

version = "v0.4"
OG_number = 0
available_cores = multiprocessing.cpu_count() #counts how many cores are available, to check if the user-argument for threads can be fulfilled
wstrings, estrings, lstrings = [], [], [] #warning, error and log messages respectively
out_file, PO_file=args.out_file + "_" + time.strftime("%Y%m%d%H%M%S"), args.po_resultfile
nthreads, seed, bootstraps = args.nthreads, args.seed_nr, args.bootstraps
tree_method, treebuilder_path = args.tree_method, args.treebuilder_path
min_freq = args.min_freq
raxml_prog = "raxmlHPC"
logfilename = out_file + "_PO_2_GENECONTENT_" + time.strftime("%Y%m%d%H%M%S") + ".log"
verbose = True
docontinue = True


def randomnumber():
	#returns a random integer to use as seed for rxml and pyml
	random.seed()
	return random.randint(1,2000)

def checkargs(args):
	global verbose, nthreads, seed, bootstraps, tree_method, raxml_prog
	global nthreads
	if args.no_verbose:
		verbose = False
	if not os.path.exists(PO_file) or not os.path.isfile(PO_file):
		datsANerror("ERROR: cannot find proteinortho-resultfile: " + PO_file)
	if available_cores<nthreads:
		datsAwarning("WARNING: less than " + str(nthreads) + " cores/threads available! Setting nthreads to " + str(available_cores))
		nthreads = available_cores
		
	if seed == 0:
		seed = randomnumber()
		datsAlogmessage("Setting seed randomly to " + str(seed))
	else:
		datsAlogmessage("Using " + str(seed) + " as seed")
	if tree_method != "none":
		if treebuilder_path == "":
			if which("raxmlHPC") == None and which("raxmlHPC-PTHREADS") == None and which("raxmlHPC-PTHREADS-SSE3") == None and which("raxmlHPC-SSE3") == None:
				datsANerror("ERROR: No raxmlHPC binaries found in any directory within $PATH! please provide a PATH to raxml binaries!")
				tree_method = "none"
			else:
				if which("raxmlHPC-PTHREADS-SSE3") != None:
					raxml_prog = which("raxmlHPC-PTHREADS-SSE3")
					datsAlogmessage("found " + raxml_prog + "! will use this tool") 
				elif which("raxmlHPC-PTHREADS")!=None:
					raxml_prog = which("raxmlHPC-PTHREADS")
					datsAlogmessage("found " + raxml_prog + "! will use this tool")
				else:
					datsAwarning("Warning: multithreading is only supported with 'PTHREADS'-versions of raxml. Not sure if your raxml binaries support this.\t\nif raxml calculations fail, recombile raxml with 'PTHREADS'-option")
					if which("raxmlHPC-SSE3") != None:
						raxml_prog = which("raxmlHPC-SSE3")
						datsAlogmessage("found " + raxml_prog + "! will use this tool")
					else:
						raxml_prog = which("raxmlHPC")
						datsAlogmessage("found " + raxml_prog + "! will use this tool")
				try:
					checkraxml_cline = RaxmlCommandline(raxml_prog, version = True)
					versiontext = checkraxml_cline()[0]
					startv = versiontext.find("version ")+len("version ")
					endv = versiontext[startv:].find(" ")
					version = versiontext[startv:startv+endv].split(".")
					if int(version[0]) < 8 or (int(version[0]) == 8 and int(version[1]) == 0 and int(version[2]) <20):
						datsAwarning("Warning: This script was devised for and tested with RAxML v8.0.20. Your version is v" + ".".join(version)+" !\n\tThis may very well still work, but if it doesn't it's YOUR fault!")
				except:
					datsAwarning("Warning: This script was devised for and tested with RAxML v8.0.20.\n\tNot sure which version of RAxML you're using, but it sure as hell isn't v7 or v8!\n\tThis may very well still work, but if it doesn't it's YOUR fault!")
		elif os.path.exists(treebuilder_path) and os.path.isfile(treebuilder_path):
			if nthreads > 1 and not "PTHREADS" in treebuilder_path:
				datsAwarning("Warning: multithreading is only supported with 'PTHREADS'-versions of raxml. Not sure if your choosen binaries support this.\t\nif raxml calculations fail, recombile raxml with 'PTHREADS'-option") 
			try:
				checkraxml_cline = RaxmlCommandline(version = True)
				versiontext = checkraxml_cline()[0]
				startv = versiontext.find("version ")+len("version ")
				endv = versiontext[startv:].find(" ")
				version = versiontext[startv:startv+endv].split(".")
				if int(version[0]) < 8 or (int(version[0]) == 8 and int(version[1]) == 0 and int(version[2]) < 20):
						datsAwarning("Warning: This script was devised for and tested with RAxML v8.0.20. Your version is v" + ".".join(version) + " !\n\tThis may very well still work, but if it doesn't it's YOUR fault!")
				raxml_prog=treebuilder_path
			except:
				datsAwarning("Warning: Correct raxML-version not found under " + treebuilder_path + "!\nWill NOT calculate ML trees!")
				tree_method = "none" 	#Add checks for dollop, mix and pars later
	#as soon as bootstrapping works with those tools, also add check for numpy v>7
	
def which(file):
	#locate scripts and programs in PATH
    for path in os.environ["PATH"].split(":"):
        if os.path.exists(path + "/" + file):
                return path + "/" + file
    return None

def datsAwarning(wstring):
	global wstrings
	wstrings.append(wstring)
	if verbose:
		print "\n" + wstring + "\n"
	
def datsANerror(estring):
	global estrings
	global docontinue
	estrings.append(estring)
	docontinue=False
	if verbose:
		print >>sys.stderr, "\n" + estring + "\n"

def datsAlogmessage(lstring):
	global lstrings
	lstrings.append(lstring)
	if verbose:
		print lstring

def read_PO_file(filename):
	open_PO_file = open(filename, 'r')
	firstline = True
	org_index = 3 #default index of beginning of organism-result-columns in proteinortho4+
	headers = []
	for line in open_PO_file:
		if firstline:
			if line.startswith("#species"):
				datsAwarning("WARNING! It seems you are using results from Proteinortho4. It's recommended to use Proteinortho5!\n\t(However, this script SHOULD still work)")
			elif line.startswith("# Species"):
				if verbose:
					print "\nProteinortho-results are based on Proteinortho5 or later. Good."
			else:
				datsAwarning("WARNING! Cannot clearly recognize Format of Proteinortho-results! This may produce erroneous results!\n\t For best results use Proteinortho5!") 
			headers = line.rstrip().split("\t")[org_index:]
			#dict_org_genecounts = dict.fromkeys(headers, 0) #for counting number of genes in each organism for calculation of simple difference/similarity matrizes
			columns = [[h.rstrip(),[]] for h in headers]
			firstline = False
		elif not line.startswith("#"):
			zeilentokens = line.rstrip().split("\t")[org_index:]
			if len(zeilentokens) != len(headers):
				datsANerror("ERROR: Different numbers of columns in header_line and Locus_tag_lines in " + filename)
				break
			for h in range(len(headers)):
				if zeilentokens[h] == "*":
					columns[h][1].append(0) #if no homolog = 0
				else:
					#dict_org_genecounts[headers[h]] += len(zeilentokens.split(",")) 
					columns[h][1].append(1) #if homolog = 1
					#print zeilentokens
	open_PO_file.close()
	#check if everything is ok before continuing:
	for c in columns:
		if len(c[1]) != len(columns[0][1]):
			datsANerror("ERROR: something went wrong while reading proteinortho-results. Not the same number of lines (=OGs) for all Organisms!")
	OG_number = len(columns[0][1])
	datsAlogmessage("\nrecognized " + str(OG_number) + " total 'Orthologeous Groups' (OGs) + singletons in the comparison organisms")
	OG_indices = range(1, OG_number + 1) #I know, I know! NOT the best way to keep track of the indices of each OG based on the original PO_file! 
	return headers, columns, OG_number, OG_indices
	
def filter_OGs_by_freq(minfreq, columns, column_indices):
	datsAlogmessage("Filtering out all OGs with a frequency below " + str(minfreq) + " between all comparison organisms")
	#filter out all OGs that do not occur in at leas minfreq comparison-organisms
	#remember: columns are organized as [[organism1=[header],[OGs]], [organism2=[header],[OGs]],...]
	line_index = 0
	#print "len(columns[0])="+str(len(columns[0]))
	while line_index < len(columns[0][1]):
		if verbose:
			sys.stdout.write("\rProcessing OG " + str(line_index) + " from "+str(len(columns[0][1])))
			sys.stdout.flush()
		OG_sum = 0
		for organism in columns:
			OG_sum += organism[1][line_index]
		if OG_sum < minfreq:
			for i in range(0, len(columns)):
				columns[i][1].pop(line_index) #remove all lines that do not add up to minfreq
			column_indices.pop(line_index) #adjust the inices accordinlgy
		else:
			#print "OG_sum = " + str(OG_sum)
			line_index+=1 #check next line, if this one is ok
	datsAlogmessage("\n" + str(len(columns[0][1])) + " OGs remaining after filtering")
	return columns, column_indices
	
def write_binary_matrix(outname, columns, column_indices):
	outfile = open(outname + ".tab", "w")
	headerstring = "OG"							#headerline
	for c in columns:							#headerline
		headerstring += ("\t" + c[0])				#headerline
	outfile.write(headerstring + "\n")			#headerline
	for index in range(0, len(column_indices)):		
		linestring = str(column_indices[index])
		for c in columns:
			linestring += "\t" + str(c[1][index])
		outfile.write(linestring + "\n")
	outfile.close()		
	return outfile.name

#def create_bootstrap_permutations(nr_bs, outname, columns, column_indices):
#	
def calculate_sim_matrix_professional(outname, columns, method):
	try:
		import scipy.spatial.distance as dist #check if biopythons distance calculation can be used for binary character data as well!
	except ImportError:
		datsANerror("SciPy (http://scipy.org/) needs to be installed to use the distance-calculation method \"" + method +"\"\nWithout this module only methods \"simple\" and \"none\" are available")
	print "Generating similarity/distance matrices using " + method  
	method_dict = {"braycurtis":dist.braycurtis, "canberra":dist.canberra, "cityblock":dist.cityblock, "jaccard":dist.jaccard, "euclidean":dist.euclidean}
	matrix_dict = {}
	for c1 in columns:
		matrix_dict[c1[0]] = {"similarity_dict":{}, "difference_dict":{}}
		for c2 in columns:
			matrix_dict[c1[0]]["difference_dict"][c2[0]] = method_dict[method](c1[1], c2[1])
			matrix_dict[c1[0]]["similarity_dict"][c2[0]] = 1- matrix_dict[c1[0]]["difference_dict"][c2[0]]
	return matrix_dict

def matrix_dict_to_matrix_obj(matrix_dict, headers): #headers-list is necessary to enforce original order of organisms in final matrix (dicts get all jumbled up)
	#this function converts my square tabular matrix format to a triangular format biopython matrix object 
	print "converting matrix_dict to biopython matrix_object" 
	from Bio.Phylo.TreeConstruction import _Matrix, _DistanceMatrix
	matrix_type = "difference_dict"
	matrix_list = []
	for indexa in range(0,len(headers)):
		indexb = 0
		org1 = headers[indexa]
		matrix_list.append([])
		while indexb <= indexa:
			org2 = headers[indexb]
			matrix_list[indexa].append(matrix_dict[org1][matrix_type][org2])
			indexb += 1
	matrix_obj = _DistanceMatrix(headers, matrix_list)
#	print matrix_obj
	return matrix_obj

def calculate_NJ_tree(outname, matrix_obj):
	from Bio import Phylo
	from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
	constructor = DistanceTreeConstructor()
	tree = constructor.nj(matrix_obj)
	out_file = open(outname + ".nj.tree", "w")
	Phylo.write(tree, out_file, "newick")
	out_file.close()


def calculate_sim_matrix_simple(outname, columns): 
	print "calculating simple similarity and difference matrices"
	#calculate a simple similarity/difference matrix based on the shared genecontent between the comparison organisms
	#for each comparison-pair, add up all shared OGs. Then divide the sum of each pair by the number of genes in the smaller genome of each pair
	#this provides a simple correction for genome size
	matrix_dict = {}
	sim_calc_dict = {}
	for c1 in columns:
#		print c1
		sim_calc_dict[c1[0]] = {"total_og_count":sum(c1[1]), "shared_og_dict":{}}
		matrix_dict[c1[0]] = {"similarity_dict":{}, "difference_dict":{}}
		#total_og_count=tota sum of OGs in respective organism, shared_og_dict=dictionary containing numbers ogs shared with each comparison organism
		#similarity_dict=for each comparison_organism: shared_ogs/smaller_total_og_count_of_comparison_pair
		for c2 in columns:
			sum_shared = 0
			for og_index in range(0, len(c1[1])):
				if c1[1][og_index] == 1 and c2[1][og_index] == 1:
					sum_shared += 1
			sim_calc_dict[c1[0]]["shared_og_dict"][c2[0]] = sum_shared
	print "===testing list==="
	for org1 in sim_calc_dict: #now calculate the actual similarity values
#		print org1
		for org2 in sim_calc_dict:
			print "test: " + org1 + " vs " + org2
#			print " shared: " +  str(sim_matrix_dict[org1]["shared_og_dict"][org2]) + " count_org1: " + str(sim_matrix_dict[org1]["total_og_count"]) + " count_org2: " + str(sim_matrix_dict[org2]["total_og_count"]) + " min: " + str(min(sim_matrix_dict[org1]["total_og_count"], sim_matrix_dict[org2]["total_og_count"]))
			matrix_dict[org1]["similarity_dict"][org2] = float(sim_calc_dict[org1]["shared_og_dict"][org2]) / min(sim_calc_dict[org1]["total_og_count"], sim_calc_dict[org2]["total_og_count"])
			print "similarity = " + str(float(sim_calc_dict[org1]["shared_og_dict"][org2]) / min(sim_calc_dict[org1]["total_og_count"], sim_calc_dict[org2]["total_og_count"]))
			matrix_dict[org1]["difference_dict"][org2] = 1 - matrix_dict[org1]["similarity_dict"][org2]
	return matrix_dict

def write_sim_dif_matrix(headers, outname, matrix_dict, matrix_method): #adapt for creating bootstraps (multiple marices in phylip format)
	#using the headers list to arrange the organisms as in the input data (dictionaries get all jumbled up)
	print "writing similarity and difference matrices"
	outfilesim = open(outname + "." + matrix_method + ".sim", "w")
	outfilediff = open(outname + "." + matrix_method + ".diff", "w")
	firstline = ""
	for org1 in headers:
		current_sim_line, current_diff_line = org1 , org1
		current_line = org1
		for org2 in headers:
			if firstline != None:
				firstline += "\t" + org2
			current_sim_line += "\t" + str(matrix_dict[org1]["similarity_dict"][org2])
			current_diff_line += "\t" + str(matrix_dict[org1]["difference_dict"][org2])
#			print current_diff_line
		if firstline != None:
			outfilesim.write(firstline + "\n")
			outfilediff.write(firstline + "\n")
			firstline = None
		outfilesim.write(current_sim_line + "\n")
		outfilediff.write(current_diff_line + "\n")
	return outfilesim.name, outfilediff.name

def write_binaryalignment_fasta(outputfilename, columns):
	outfile = open(outputfilename+".fas", "w")
	if verbose:
		print "producing alignment file '" + outfile.name + "' in fasta format"
	for c in columns:
		outfile.write(">" + c[0] + "\n")
		for line in c[1]:
			outfile.write(str(line))
		outfile.write("\n")
	outfile.close()
	return outfile.name
	
def write_binaryalignment_phylip(outputfilename, columns):
	outfile = open(outputfilename + ".relaxed.phy","w")
	if verbose:
		print "producing alignment file '" + outfile.name + "' in relaxed sequential phylip format"
	outfile.write(str(len(columns)) + " " + str(len(columns[0][1])))
	for c in columns:
		#print c[0]+c[0]
		linestring = "\n" + c[0] + " "
		#print "line"+str(c[1][0])+"line"
		for line in c[1]:
			linestring += str(line)
		#print linestring
		outfile.write(linestring)
	outfile.close()
	print "writing: " + outfile.name
	return outfile.name

def call_raxml_rapidbs(alignmentfile, outputfilename, seed, parameters): #parameters should be a dictionary (This dictionary thing was introduced, so that the script can be more easily adapted to accept custom commandline-parameters for raxml by the user)
	nr_threads = 4
	if "-T" in parameters:
		nr_threads=parameters["-T"]
	bootstraps=1000
	if "-N" in parameters:
		bootstraps = parameters["-N"]
	elif "-#" in parameters:
		bootstraps = parameters["-#"]
	
	datsAlogmessage("Calculating phylogenies: 'rapid bootstrap analyses and search for best-scoring ML Tree in one run' using raxmlHPC")
	try:
		outname = "GENECONTENT_rapidBS" + str(bootstraps) + "_minfreq"+str(min_freq) + "_" + time.strftime("%Y%m%d%H%M%S") + "_" + "final_tree"	
		raxml_cline = RaxmlCommandline(raxml_prog, sequences = alignmentfile, algorithm = "a", model = "BINGAMMA", name = outname, parsimony_seed = seed, rapid_bootstrap_seed = seed, num_replicates = bootstraps, threads = nr_threads) 
		datsAlogmessage("-->" + str(raxml_cline))
		raxml_cline()
		datsAlogmessage("-->SUCCESS")
		#the resultfiles will be: "RAxML_bipartitions.rapidBS_final_tree" and "RAxML_bipartitionsBranchLabels.rapidBS_final_tree"
		#Labels on nodes or branches, respectively
		outputfiles = ["RAxML_bipartitions."+outname, "RAxML_bipartitionsBranchLabels."+outname]
	except Exception as e:
		#note for future versions: create a custom exception "raxml_except" instead of using this construct over and over!
		datsAlogmessage("-->FAILURE")
		datsAwarning("WARNING: rapid bootstrap analyses and search for best-scoring ML tree using raxmlHPC has failed!\n\t" + str(e))	
		return None
	return outputfiles
				
def call_raxml_bs(alignmentfile, outputfilename, seed, parameters):
	nr_threads = 4
	if "-T" in parameters:
		nr_threads = parameters["-T"]
	bootstraps = 1000
	if "-N" in parameters:
		bootstraps = parameters["-N"]
	elif "-#" in parameters:
		bootstraps = parameters["-#"]
	datsAlogmessage("Calculating phylogenies: Thorough bootstrap analyses with raxml")
		
	datsAlogmessage("\tDetermining best ML tree of 20 raxmlHPC runs") 
	try:
		raxml_cline = RaxmlCommandline(raxml_prog, model = "BINGAMMA", name = "best_delme_tempfile", parsimony_seed = seed, num_replicates = 20, sequences = alignmentfile, threads = nr_threads) 
		datsAlogmessage("\t-->"+str(raxml_cline))
		raxml_cline()
		#the resultfile will be :"RAxML_bestTree.best_delme_tempfile"
		datsAlogmessage("\t-->SUCCESS")
	except Exception as e:
		datsAlogmessage("\t--FAILURE")
		datsAwarning("WARNING: thorough bootstrap analyses using raxmlHPC failed!\n\t"+str(e))
		return None
	
	datsAlogmessage("\tDoing bootstrap analyses with "+str(bootstraps)+" runs using raxmlHPC")
	try:
		raxml_cline = RaxmlCommandline(raxml_prog, model = "BINGAMMA", sequences = alignmentfile, name = "boot_delme_tempfile", parsimony_seed = seed, bootstrap_seed = seed, num_replicates = bootstraps, threads = nr_threads)
		datsAlogmessage("\t-->" + str(raxml_cline))
		raxml_cline()
		#the resultfile will be: "RAxML_bootstrap.boot_delme_tempfile"
		datsAlogmessage("\t-->SUCCESS")
	except Exception as e:
		datsAlogmessage("\t-->FAILURE")
		datsAwarning("WARNING: thorough bootstrap analyses using raxmlHPC failed!\n\t" + str(e))
		return None
		
	datsAlogmessage("\tDrawing bipartitions of bootstrap trees onto best ML tree using raxmlHPC")
	try:
		outname = "GENECONTENT_BS" + str(bootstraps) + "_minfreq" + str(min_freq) + time.strftime("%Y%m%d%H%M%S") + "_" + "final_tree"
		raxml_cline = RaxmlCommandline(raxml_prog, model = "BINGAMMA", parsimony_seed = seed, algorithm = "b", starting_tree = "RAxML_bestTree.best_delme_tempfile", bipartition_filename="RAxML_bootstrap.boot_delme_tempfile", name=outname)
		datsAlogmessage("\t-->"+str(raxml_cline))
		raxml_cline()
		#The resultfiles will be: RAxML_bipartitions.final_tree" and "RAxML_bipartitionsBranchLabels.final_tree"
		outputfiles=["RAxML_bipartitions."+outname,"RAxML_bipartitionsBranchLabels."+outname]
		datsAlogmessage("\t-->SUCCESS")	
	except Exception as e:
		datsAlogmessage("\t-->FAILURE")
		datsAwarning("WARNING: thorough bootstrap analyses using raxmlHPC failed!\n\t"+str(e))
		return None
	return outputfiles

def call_raxml_nobs(alignmentfile, outputfilename, seed, parameters):
	nr_threads=4
	if "-T" in parameters:
		nr_threads=parameters["-T"]
	try:
		outname="GENECONTENT_raxml_"+"_minfreq"+str(min_freq)+time.strftime("%Y%m%d%H%M%S")+"_"+"final_tree"
		datsAlogmessage("Calculating phylogeny: Determining best ML tree of 20 raxmlHPC runs")
		raxml_cline=RaxmlCommandline(raxml_prog, sequences = alignmentfile, model = "BINGAMMA", name = outname,  parsimony_seed = seed, num_replicates = 20, threads = nr_threads)
		datsAlogmessage("\t-->" + str(raxml_cline))
		raxml_cline()
		#the resultfile will be :"RAxML_bestTree.final_tree"
		outputfiles = ["RAxML_bestTree." + outname]
		datsAlogmessage("\tSUCCESS")
		print "deleting temporary files"
		for delfile in os.listdir("."):
			if delfile.startswith("RAxML_") and "." + outname + ".RUN." in delfile:
				#print "deleting temp-file: "+delfile
				os.remove(delfile)
	except Exception as e:
		datsAlogmessage("\t-->FAILURE")
		datsAwarning("WARNING: searching for best ML tree using raxmlHPC failed!\n\t"+str(e))
		return None
	return outputfiles

def relaxed_2_strict_phylip(input_file, output_file): #create library with uniqe 10 character names
	datsAlogmessage("converting relaxed phylip format ('"+input_file+"') to strict sequential phylip format('" + output_file + "')")
	datsAlogmessage("CAREFUL: Names will be shortened to 10 characters! This may affect branch labeling!")
	try:
		relaxed_file = open(input_file, "r")
		relaxed_align = AlignIO.read(relaxed_file, "phylip-relaxed")
		relaxed_file.close()
		strict_file = open(output_file, "w")
		AlignIO.write(relaxed_align, "phylip-sequential")
		strict_file.close()
		return outputfile
	except Exception as e:
		datsAlogmessage("-->FAILED")
		datsAwarning("Could not convert relaxed phylip format to strict phylip format!\nPHYLIP based programs (fpar, fdollop, fmix) can NOT be used!")
		return None
		
def call_fpars(alignmentfile, out_file, seed):
	#input should be a temporal "strict" phylip file, NOT the relaxed phylip file that is normally generated!
	#important: first you have to generate a temporal strict phylip file with 10-character placeholders for each name (stored as a dict).
	#use this temporal input file as input for fpars
	#then read in the outputfiles of fpars and translate the names back into the original names (->convert back to relaxed phylip file)
	#but do this outside of this method in order to enable to use the seqboot and condense tools also
	#UPDATE: it appears the use of discrete (binary) characters is NOT possible in Fseqboot.
	datsAlogmessage("calculating phylogeny: Discrete character parsimony using fpars (EMBOSS implementation of the PHYLIP package)")
	fpars_args=['-infile '+alignmentfile, '-outfile '+out_file+"_fpars.log", '-auto', "-seed "+str(seed)]
	fpars_command=["fpars"]+fpars_args
	datsAlogmessage("-->"+" ".join(fpars_args))
	try:
		call(fpars_command)
		datsAlogmessage("-->SUCCESS")
		#the output will be a logfile named "<output>.log" and a treefile named "strict.treefile" (because the inputfilename will be automatically concatenated to "strict.")
		#therefore rename the treefile:
		os.rename("strict.treefile",out_file + "_fpars.treefile")
		return [out_file + "_fpars.log", out_file + "_fpars.treefile"]
	except Exception as e:
		datsAlogmessage("-->FAILURE")
		datsAwarning("WARNING: Discrete character parsimony using fpars failed!\n\t" + str(e))
		return None

def call_fmix(alignmentfile, out_file, seed):
	#input should be a temporal "strict" phylip file, NOT the relaxed phylip file that is normally generated!
	#important: first you have to generate a temporal strict phylip file with 10-character placeholders for each name (stored as a dict).
	#use this temporal input file as input for fmix
	#then read in the outputfiles of fpars and translate the names back into the original names (->convert back to relaxed phylip file)
	#but do this outside of this method in order to enable to use the seqboot and condense tools also
	#UPDATE: it appears the use of discrete (binary) characters is NOT possible in Fseqboot.
	datsAlogmessage("calculating phylogeny: Discrete character parsimony using fmix with wagner parsimony (EMBOSS implementation of the PHYLIP package)")
	fpmix_args=['-infile ' + alignmentfile, '-outfile ' + out_file + "_fmix.log", "-method w", '-auto', "-seed " + str(seed)]
	fmix_command=["fmix"] + fmix_args
	datsAlogmessage("-->" + " ".join(fpars_args))
	try:
		call(fmix_command)
		datsAlogmessage("-->SUCCESS")
		#the output will be a logfile named "<output>.log" and a treefile named "strict.treefile" (because the inputfilename will be automatically concatenated to "strict.")
		#therefore rename the treefile:
		os.rename("strict.treefile",out_file+"_fmix.treefile")
		return [out_file+"_fmix.log", out_file+"_fmix.treefile"]
	except Exception as e:
		datsAlogmessage("-->FAILURE")
		datsAwarning("WARNING: Discrete character parsimony using fmix failed!\n\t"+str(e))
		return None
	

def call_fdollop(alignmentfile, out_file, seed):
	#input should be a temporal "strict" phylip file, NOT the relaxed phylip file that is normally generated!
	#important: first you have to generate a temporal strict phylip file with 10-character placeholders for each name (stored as a dict).
	#use this temporal input file as input for fmix
	#then read in the outputfiles of fpars and translate the names back into the original names (->convert back to relaxed phylip file)
	#but do this outside of this method in order to enable to use the seqboot and condense tools also
	#UPDATE: it appears the use of discrete (binary) characters is NOT possible in Fseqboot.
	datsAlogmessage("calculating phylogeny: Discrete character parsimony using fdollop with dollo parsimony (EMBOSS implementation of the PHYLIP package)")
	fpmix_args=['-infile '+alignmentfile, '-outfile '+out_file+"_fdollop.log", "-method d", '-auto', "-seed "+str(seed)]
	fmix_command=["fdollop"]+fmix_args
	datsAlogmessage("-->"+" ".join(fpars_args))
	try:
		call(fmix_command)
		datsAlogmessage("-->SUCCESS")
		#the output will be a logfile named "<output>.log" and a treefile named "strict.treefile" (because the inputfilename will be automatically concatenated to "strict.")
		#therefore rename the treefile:
		os.rename("strict.treefile",out_file+"_fdollop.treefile")
		return [out_file+"_fdollop.log", out_file+"_fdollop.treefile"]
	except Exception as e:
		datsAlogmessage("-->FAILURE")
		datsAwarning("WARNING: Discrete character parsimony using fdollop failed!\n\t"+str(e))
		return None
	#call function to convert indexed names in treefile to original names (could be done by calling "sed -i" on the resultfiles)
	#must convert treefile to other filename and return that filename
	#return something
	
def call_phyml():
	# implement phylip tools (neighborjoining etc) for distance based tree inference (phylip neighbor etc)
	pass

def write_logfile(logfilename):
	logfile = open(logfilename,'w')
	logfile.write("PO_2_GENECONTENT.py logfile")
	logfile.write("\n" + time.strftime("%c"))
	logfile.write("\nMessages\n:")
	for l in lstrings:
		logfile.write(l + "\n")
	logfile.write("\n-----------\nWarnings:\n")
	for w in wstrings:
		logfile.write(w + "\n")
	logfile.write("\n------------\nErrors:\n")
	for e in estrings:
		logfile.write(e + "\n")
	logfile.close()

def main():
	alignment_files = []
	tree_files = []
	if verbose:
		print "\n ==PO_2_GENECONTENT.py " + version + " by John Vollmers==\n"
	checkargs(args)
	if docontinue:
		try:
			headers, columns, OG_number , column_indices = read_PO_file(PO_file)
			if min_freq > 1:
				columns, column_indices = filter_OGs_by_freq(min_freq, columns, column_indices)
			if args.outfasta and docontinue:
				alignment_files.append(write_binaryalignment_fasta(out_file, columns))
			if args.outdiff != "none" and docontinue:
				if args.outdiff == "simple":
					matrix_dict = calculate_sim_matrix_simple(out_file, columns)
				else:
					matrix_dict = calculate_sim_matrix_professional(out_file, columns, args.outdiff)
				alignment_files.extend(write_sim_dif_matrix(headers, out_file, matrix_dict, args.outdiff))
				matrix_obj = matrix_dict_to_matrix_obj(matrix_dict, headers)
				calculate_NJ_tree(out_file + ".args.outdiff", matrix_obj)
			if args.outmatrix and docontinue:
				alignment_files.append(write_binary_matrix(out_file, columns, column_indices))
			if docontinue:
				alignment_files.append(write_binaryalignment_phylip(out_file, columns))
			if tree_method == "raxml_rapidbs" and docontinue:
				tree_files.extend(call_raxml_rapidbs(alignment_files[-1], out_file, seed, {"-N":bootstraps, "-T":nthreads}))
			elif tree_method == "raxml_bs" and docontinue:
				tree_files.extend(call_raxml_bs(alignment_files[-1], out_file, seed, {"-N":bootstraps, "-T":nthreads}))
			elif tree_method == "raxml" and docontinue:
				print "number of alignmentfiles: " + str(len(alignment_files))
				tree_files.extend(call_raxml_nobs(alignment_files[-1], out_file, seed, {"-T":nthreads}))
#			elif args.fpars or args.fmix or args.fdollop and docontinue: #remove this sometime
#				alignment_files.append(relaxed_2_strict_phylip(alignment_files[-1], alignment_files[-1].replace(".relaxed.", ".strict.")))
#				if alignment_files[-1]!=None:
#					if args.fpars:
#						treefiles.extend(call_fpars(alignment_files[-1], outfile, seed))
#					if args.fdollop:
#						treefiles.extend(call_fdollop(alignment_files[-1], outfile, seed))
#					if args.fmix:
#						treefiles.extend(call_fmix(alignment_files[-1], outfile,seed))
			if docontinue:
				if verbose:
					print "================================\FINISHED!"	
				datsAlogmessage("Created the following alignment-files:" + "\n\t-".join(alignment_files))
				datsAlogmessage("Created the following tree-files:" + "\n\t-".join(tree_files))
		except Exception as e:
			datsAlogmessage(str(e))
			for frame in traceback.extract_tb(sys.exc_info()[2]):
				fname,lineno,fn,text = frame
				print "Error in %s on line %d" % (fname, lineno)
		finally:
			datsAlogmessage("cleaning up...")
			for delfile in os.listdir("."):
				if "delme_tempfile" in delfile:
					os.remove(delfile)
					
	if len(wstrings) > 0 or len(estrings)>0:
		if len(wstrings) > 0:
			print "\nThere were Warnings!"
		if len(estrings) > 0:
			print "\nThere were Errors!"
	print "See " + logfilename + " for details\n"
	write_logfile(logfilename)

main()
