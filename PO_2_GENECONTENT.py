#!/usr/bin/python
#created 20.03.15 by John Vollmers
#from Bio.Phylo.TreeConstruction import * #for using the directly biopython-implemented Distance(UPGMA + NJ) and parsimony methods. Not done/necessary yet

import os, sys, logging, argparse, time, multiprocessing, random, traceback, Bio
from Bio import AlignIO, SeqIO
from subprocess import call
from Bio.Phylo.Applications import RaxmlCommandline, PhymlCommandline

myparser=argparse.ArgumentParser(description="\n==PO_2_GENECONTENT.py v1.01 by John Vollmers==\nCreates discrete binary character matrices and phylogenetic trees based on the gene content (presence/absence of orthologues) in comparison organisms.\nThis script is supposed to be part of a pipeline consisting of:\n\
\tA.)Conversion of Genbank/Embl-Files to ANNOTATED(!) Fastas\n\t   using CDS_extractor.pl by Andreas Leimbach\n\
\tB.)Calculation of orthologs and paralogs using proteinortho5\n\t   (WITH the '-single' and '-self' arguments!)\n\
\tC.)The creation of discrete binary character marices\n\t   (and optionally phylogenetic trees) based on:\n\t\t-the fasta sequences of step A\n\
\t\t-the proteinortho5-results from step B\n", formatter_class=argparse.RawTextHelpFormatter)
myparser.add_argument("-po", "--proteinortho", action = "store", dest = "po_file", help = "(String) file with proteinortho5 results", required = True)
myparser.add_argument("-s", "--silent", action = "store_true", dest = "no_verbose", help = "non-verbose mode")
myparser.add_argument("-o", "--out", action = "store", dest="out_file", default = "output", help = "Basename for result-files\nDefault='output'\nWill create a phylip-file of the aligned discrete binary character data by default")
myparser.add_argument("-ofas", "--out_fasta", action = "store_true", dest = "outfasta", default = False, help = "Create a file with the aligned discrete binary character data in Fasta format\n(in addition to the default phylip file)")
myparser.add_argument("-omat", "--out_matrix", action = "store_true", dest = "outmatrix", default = False, help = "Create a file with the aligned discrete character data in tabular format\n(in addition to the default phylip file)")
myparser.add_argument("-odiff", action = "store", dest = "outdiff", choices = ["simple", "none", "braycurtis", "canberra", "cityblock", "jaccard", "euclidean"], default = "simple", help = "Create difference/similarity matrixes for distance matrix based phylogeny inference.\n\t- simple : Simple normalized distance matrix\n\t  (distance=Nr shared OGs/total Nr OGs in SMALLER organism)\n\t- braycurtis : Bray-Curtis distance\n\t- canberra : Canberra distance\n\t- cityblock : Manhatten distance\n\t- euclidean : Euclidean distance\n\t- jaccard : Jaccard-Needham dissimilarity\nDefault: \"simple\"") #only simply algorithm implemented for now. Add more complex methods and bootstrapping options
myparser.add_argument("-cpu", "--cpu", action = "store", dest = "ncpus", type = int, default = 4, help = "Number of cpus to use\nDefault=4 (or maximum number of available cores, if less than 4 cores available)")
myparser.add_argument("-sd", "--seed", action = "store", dest = "seed", type = int, default = 0, help = "Integer to provide as seed for RAxML or PhyML\n0=seed generated randomly\nDefault=random seed")
myparser.add_argument("-mt", "--make_tree", action = "store", dest = "tree_method", choices = ["raxml", "raxml_bs", "raxml_rapidbs", "nj", "nj_bs", "none"], default = "none", help = "Generate ML phylogenetic trees using RAxML with \"new rapid hill climbing\"\nand substitution model: \"BINGAMMA\"\n\t-\"raxml\": single tree without bootstraps (using new rapid hill climbing)\n\t-\"raxml_bs\": thorough bootstrap analyses and search for best ML tree\n\t-\"raxml_rapidbs\": rapid bootstrap analyses and search for best ML tree\n\t in one run\n\t-\"nj\": Neighbor joining\n\t-\"none\"\nDefault=none")
myparser.add_argument("-tbp","--tree_builder_path", action="store", dest = "treebuilder_path", default = "", help = "Path to treebuilder (currently only raxml supported) if not listed in $PATH")
myparser.add_argument("-bs", "--bootstraps", action = "store", dest = "bootstraps", type = int, default = 1000, help = "Number of times to resample for bootstrapping\nonly bootstrapped trees will be produced, not the permutated data matrices")
myparser.add_argument("--min_freq", action = "store", dest = "min_freq", type=int, default = 2, help = "Minimum frequency for Orthologeous Groups across all comparison-organisms\nfor Group to be considered in Analyses\nDefault=2 (exclude all singletons)\nRecommended when working with highly fragmented genomes or dubious ORF-finding")
myparser.add_argument("--debug", action = "store_true", dest = "debug", default = False, help = "Log extra info for debugging")
#myparser.add_argument("--ignore_columns", action="store", dest="ignore_column_list", default=None, help="indexes (starting with 1) of Organisms to be ignored when counting frequencies of Orthologeous Groups (Tip: ignore all but one representative member of overrepresented species/strains)")
args = myparser.parse_args()

#set global variables (in addition to args)
version = "v0.7c"
OG_number = 0
available_cores = multiprocessing.cpu_count() #counts how many cores are available, to check if the user-argument for cpus can be fulfilled
out_file = args.out_file + "_" + time.strftime("%Y%m%d%H%M%S")
logtime = time.strftime('%H:%M %Z on %b %d, %Y')
#ncpus, seed, bootstraps = args.ncpus, args.seed_nr, args.bootstraps
#tree_method, treebuilder_path = args.tree_method, args.treebuilder_path
#min_freq = args.min_freq
raxml_prog = "raxmlHPC"
logfilename = "%s_PO_2_GENECONTENT.log" % out_file
#verbose = True
#docontinue = True

def setlogger(logfilename): #shall replace datsanerror, datsanlogmessage  datswarning

	#format logger
	logFormatter = logging.Formatter("[%(levelname)s]  %(message)s")
	mylogger = logging.getLogger()
	if args.debug:
		mylogger.setLevel(logging.DEBUG)
	else:
		mylogger.setLevel(logging.INFO)
		
	#set logger to write to logfile
	fileHandler = logging.FileHandler(logfilename)
	fileHandler.setFormatter(logFormatter)
	mylogger.addHandler(fileHandler)
	
	#set loger to write to stderr also
	consoleHandler = logging.StreamHandler()
	consoleHandler.setFormatter(logFormatter)
	if args.no_verbose:
		consoleHandler.setlevel(logging.WARNING)
	mylogger.addHandler(consoleHandler)
	
	return mylogger

def randomnumber(min, max):
	#returns a random integer to use as seed for rxml and pyml
	#mylogger.debug("randomnumber(%s, %s)" %(min, max))
	random.seed()
	return random.randint(min,max)

def checkargs():
	global args
	mylogger.debug("checkargs(%s)" %(args,))
	global raxml_prog
#	if args.no_verbose:
#		verbose = False
	if not os.path.exists(args.po_file) or not os.path.isfile(args.po_file):
		estring = "Cannot find proteinortho-resultfile: %s" % args.po_file
		mylogger.error(estring)
		raise IOError(estring)
	if available_cores < args.ncpus:
		mylogger.warning("less than %s cores/threads available! Setting number of cpus to %s" %(args.ncpus, available_cores))
		args.ncpus = available_cores
		
	if args.seed == 0 and "raxml" in args.tree_method:
		args.seed = randomnumber(1, 2000)
		mylogger.info("Setting seed randomly to %s" % args.seed)
	elif "raxml" in args.tree_method:
		mylogger.info("Using %s as seed" % args.seed)
	if "raxml" in args.tree_method:
		if args.treebuilder_path == "":
			if which("raxmlHPC") == None and which("raxmlHPC-PTHREADS") == None and which("raxmlHPC-PTHREADS-SSE3") == None and which("raxmlHPC-SSE3") == None:
				estring = "No raxmlHPC binaries found in any directory within $PATH! please provide a PATH to raxml binarieswhen using raxml for treebuilding!"
				mylogger.error(estring)
				raise OSError(estring)
			else:
				if which("raxmlHPC-PTHREADS-SSE3") != None:
					raxml_prog = which("raxmlHPC-PTHREADS-SSE3")
					mylogger.info("found %s ! will use this tool" % raxml_prog)
				elif which("raxmlHPC-PTHREADS")!=None:
					raxml_prog = which("raxmlHPC-PTHREADS")
					mylogger.info("found %s ! will use this tool" % raxml_prog)
				else:
					mylogger.warning("Warning: multithreading is only supported with 'PTHREADS'-versions of raxml. Not sure if your raxml binaries support this.\t\nif raxml calculations fail, recombile raxml with 'PTHREADS'-option")
					if which("raxmlHPC-SSE3") != None:
						raxml_prog = which("raxmlHPC-SSE3")
						mylogger.info("found %s ! will use this tool" % raxml_prog)
					else:
						raxml_prog = which("raxmlHPC")
						mylogger.info("found %s ! will use this tool" % raxml_prog)
				try:
					checkraxml_cline = RaxmlCommandline(raxml_prog, version = True)
					versiontext = checkraxml_cline()[0]
					startv = versiontext.find("version ")+len("version ")
					endv = versiontext[startv:].find(" ")
					version = versiontext[startv:startv+endv].split(".")
					if int(version[0]) < 8 or (int(version[0]) == 8 and int(version[1]) == 0 and int(version[2]) <20):
						mylogger.warning("Warning: This script was devised for and tested with RAxML v8.0.20. Your version is v" + ".".join(version)+" !\n\tThis may very well still work, but if it doesn't it's YOUR fault!")
				except:
					mylogger.warning("Warning: This script was devised for and tested with RAxML v8.0.20.\n\tNot sure which version of RAxML you're using, but it sure as hell isn't v7 or v8!\n\tThis may very well still work, but if it doesn't it's YOUR fault!")
		elif os.path.exists(args.treebuilder_path) and os.path.isfile(args.treebuilder_path):
			if args.ncpus > 1 and not "PTHREADS" in args.treebuilder_path:
				mylogger.warning("Warning: multithreading is only supported with 'PTHREADS'-versions of raxml. Not sure if your choosen binaries support this.\t\nif raxml calculations fail, recombile raxml with 'PTHREADS'-option")
			try:
				checkraxml_cline = RaxmlCommandline(version = True)
				versiontext = checkraxml_cline()[0]
				startv = versiontext.find("version ")+len("version ")
				endv = versiontext[startv:].find(" ")
				version = versiontext[startv:startv+endv].split(".")
				if int(version[0]) < 8 or (int(version[0]) == 8 and int(version[1]) == 0 and int(version[2]) < 20):
						mylogger.warning("Warning: This script was devised for and tested with RAxML v8.0.20. Your version is v" + ".".join(version) + " !\n\tThis may very well still work, but if it doesn't it's YOUR fault!")
				raxml_prog=args.treebuilder_path
			except:
				mylogger.warning("Warning: Correct raxML-version not found under %s !\nWill NOT calculate ML trees!" % args.treebuilder_path)

def which(afile):
#locate scripts and programs in PATH
	for path in os.environ["PATH"].split(":"):
		if os.path.exists(path + "/" + afile):
			return path + "/" + afile
	return None

def read_PO_file(filename):
	mylogger.debug("read_PO_file(%s)" % filename)
	
	open_PO_file = open(filename, 'r')
	firstline = True
	org_index = 3 #default index of beginning of organism-result-columns in proteinortho4+
	headers = []
	
	for line in open_PO_file:
		if firstline:
			if line.startswith("#species"):
				mylogger.warning("WARNING! It seems you are using results from Proteinortho4. It's recommended to use Proteinortho5!\n\t(However, this script SHOULD still work)")
			elif line.startswith("# Species"):
				if not args.no_verbose:
					print "\nProteinortho-results are based on Proteinortho5 or later. Good."
			else:
				mylogger.warning("WARNING! Cannot clearly recognize Format of Proteinortho-results! This may produce erroneous results!\n\t For best results use Proteinortho5!")
			headers = line.rstrip().split("\t")[org_index:]
			#dict_org_genecounts = dict.fromkeys(headers, 0) #for counting number of genes in each organism for calculation of simple difference/similarity matrizes
			columns = [[h.rstrip(),[]] for h in headers]
			firstline = False
		elif not line.startswith("#"):
			zeilentokens = line.rstrip().split("\t")[org_index:]
			if len(zeilentokens) != len(headers):
				estring = "Different numbers of columns in header_line and Locus_tag_lines in %s" % filename
				mylogger.error(estring)
				raise Runtime(estring)
				break
			for h in range(len(headers)):
				if zeilentokens[h] == "*":
					columns[h][1].append(0) #if no homolog = 0
				else:
					#dict_org_genecounts[headers[h]] += len(zeilentokens.split(","))
					columns[h][1].append(1) #if homolog = 1
					
	open_PO_file.close()
	
	#check if everything is ok before continuing:
	for c in columns:
		if len(c[1]) != len(columns[0][1]):
			estring = "ERROR: something went wrong while reading proteinortho-results. Not the same number of lines (=OGs) for all Organisms!"
			mylogger.error(estring)
			raise RuntimeException(estring)
			
	OG_number = len(columns[0][1])
	mylogger.info("\nrecognized " + str(OG_number) + " total 'Orthologeous Groups' (OGs) + singletons in the comparison organisms")
	OG_indices = range(1, OG_number + 1) #I know, I know! NOT the best way to keep track of the indices of each OG based on the original PO_file!
	
	return headers, columns, OG_number, OG_indices

def filter_OGs_by_freq(minfreq, columns, column_indices):
	mylogger.info("Filtering out all OGs with a frequency below %s between all comparison organisms" % minfreq)
	#filter out all OGs that do not occur in at leas minfreq comparison-organisms
	#remember: columns are organized as [[organism1=[header],[OGs]], [organism2=[header],[OGs]],...]
	
	line_index = 0
	#print "len(columns[0])="+str(len(columns[0]))
	
	while line_index < len(columns[0][1]):
		if not args.no_verbose:
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
			
	mylogger.info("\n" + str(len(columns[0][1])) + " OGs remaining after filtering")
	return columns, column_indices

def write_binary_matrix(outname, columns, column_indices):
	outfile = open(outname + ".tab", "w")
	
	#first just write the column headers
	headerstring = "OG"							#headerline
	for c in columns:							#headerline
		headerstring += ("\t" + c[0])			#headerline
	outfile.write(headerstring + "\n")			#headerline
	
	#then write the actual genecontent data matrix (line for line)
	for index in range(0, len(column_indices)):
		linestring = str(column_indices[index])
		for c in columns:
			linestring += "\t" + str(c[1][index])
		outfile.write(linestring + "\n")
	
	outfile.close()
	return outfile.name

def make_single_bootstrap_tree(columns, maxlen, method, headers, mp_output):
	debugtime_a = time.time()
	current_permutation = [[col[0], []] for col in columns]
	
	#permute original data-matrix
	for c in range(maxlen):
		rand_sample = randomnumber(0, maxlen-1)
		for org_index in range(0, len(columns)):
			current_permutation[org_index][1].append(columns[org_index][1][rand_sample])
	debugtime_b = time.time()
	print "created another bootstrap permutation in %s" % calc_timediff(debugtime_a, debugtime_b)
	#directly create distance matrix from permutation
	current_permutation_matrix_obj = matrix_dict_to_matrix_obj(calculate_matrix(current_permutation, method), headers)
	debugtime_c = time.time()
	print "created another bootstrap permutaion matrix in %s (total)" % calc_timediff(debugtime_a, debugtime_c)
	#create NJ tree from distance matrix
	current_permutation_tree = calculate_NJ_tree(current_permutation_matrix_obj)
	debugtime_d = time.time()
	print "created another bootstrap permutation tree in %s (total)" % calc_timediff(debugtime_a, debugtime_d)
	mp_output.put(current_permutation_tree)

def create_bootstrap_permutations(columns, method, headers):
	mylogger.info("creating %d bootstrap_permutations of binary character matrix (and corresponding bootstrap tree)\nusing %d cpus" %(args.bootstraps, args.ncpus))
	
	full_thread_mp_groups = args.bootstraps // args.ncpus
	remaining_mp_group_threads = args.bootstraps % args.ncpus
	maxlen = len(columns[0][1])
	permutation_tree_list = []
	bs_index = 0
	
	for mp_group in range(full_thread_mp_groups):
		permutation_tree_list.extend(run_multiprocess_group(columns, method, headers, args.ncpus))
		sys.stdout.write("\rcreated %d of %d bootstrap permutations" %(len(permutation_tree_list), args.bootstraps))
		sys.stdout.flush()
		
	#finish off any remaining permutation jobs (in case the numbe of bootstraps was not evenly divisible by the number of cpus)
	if remaining_mp_group_threads > 0:
		permutation_tree_list.extend(run_multiprocess_group(columns, method, headers, remaining_mp_group_threads))
		sys.stdout.write("\rcreated %d of %d bootstrap permutations" %(len(permutation_tree_list), args.bootstraps))
		sys.stdout.flush()
		
	return permutation_tree_list

def run_multiprocess_group(columns, method, headers, cpus): #must start working with classes to avoid such long argument lists
	if __name__ == '__main__': #just making sure function is only called within its intended context
		maxlen = len(columns[0][1])
		group_results = []
		
		mp_output = multiprocessing.Queue()
		processes = [multiprocessing.Process(target=make_single_bootstrap_tree, args=(columns, maxlen, method, headers, mp_output)) for x in range(cpus)]
		
		for p in processes:
			p.start()
			
		for p in processes:
			p.join()
			
		group_results = [mp_output.get() for p in processes]
		return group_results
		
	else:
		raise RuntimeError("FORBIDDEN TO CALL THIS (MULTIPROCESSING) FUNCTION FROM AN EXTERNAL MODULE\n-->ABORTING")

def calculate_bootstrapped_NJ_tree(ref_tree, permutation_trees):
	mylogger.debug("calculate_bootstrapped_NJ_tree(ref_tree, permutation_trees)")
	from Bio.Phylo.Consensus import get_support
	
	bs_tree = get_support(ref_tree, permutation_trees)
	
	return bs_tree

def calculate_matrix(columns, method):

	if method != "none" and method != "simple":
		matrix_dict = calculate_sim_matrix_professional(columns, method)
	elif method == "simple":
		matrix_dict = calculate_sim_matrix_simple(columns)
		
	return matrix_dict

def calculate_sim_matrix_professional(columns, method):

	try:
		import scipy.spatial.distance as dist #check if biopythons distance calculation can be used for binary character data as well!
	except ImportError:
		datsANerror("SciPy (http://scipy.org/) needs to be installed to use the distance-calculation method \"" + method +"\"\nWithout this module only methods \"simple\" and \"none\" are available")
		
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
#	mylogger.info("converting matrix_dict to biopython matrix_object")
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
	return matrix_obj

def calculate_NJ_tree(matrix_obj):
	from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
	
	constructor = DistanceTreeConstructor()
	tree = constructor.nj(matrix_obj)
	
	return tree

def remove_internal_tree_labels(thistree):
	mylogger.debug("remove_internal_tree_labels(tree)")
	mylogger.info("\nRemoving internal node labels that do nor refer to confidence values from final tree")
	
	for inode in thistree.get_nonterminals():
		inode.name = None
		
	return thistree

def write_nj_tree(outname, matrix_method, bootstraps, main_tree):
	mylogger.debug("write_nj_tree(%s, %s, %s, main_tree)" %(outname, matrix_method, args.bootstraps))
	
	from Bio import Phylo
	outname += ".%s.%s" %(matrix_method, args.tree_method)
	if args.tree_method == "nj_bs":
		outname += ".%s" % args.bootstraps
	outname += ".tree.newick"
	mylogger.info("writing treefile to %s" % outname)
	
	out_file = open(outname, "w")
	Phylo.write(main_tree, out_file, "newick")
	out_file.close()

	return outname

def calculate_sim_matrix_simple(columns):
	#calculate a simple similarity/difference matrix based on the shared genecontent between the comparison organisms
	#for each comparison-pair, add up all shared OGs. Then divide the sum of each pair by the number of genes in the smaller genome of each pair
	#this provides a simple correction for genome size
	matrix_dict = {}
	sim_calc_dict = {}
	
	for c1 in columns:
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
			
	for org1 in sim_calc_dict: #now calculate the actual similarity values
		for org2 in sim_calc_dict:
			matrix_dict[org1]["similarity_dict"][org2] = float(sim_calc_dict[org1]["shared_og_dict"][org2]) / min(sim_calc_dict[org1]["total_og_count"], sim_calc_dict[org2]["total_og_count"])
			matrix_dict[org1]["difference_dict"][org2] = 1 - matrix_dict[org1]["similarity_dict"][org2]
			
	return matrix_dict

def write_sim_dif_matrix(headers, outname, matrix_dict, matrix_method): #adapt for creating bootstraps (multiple marices in phylip format)
	#using the headers list to arrange the organisms as in the input data (dictionaries get all jumbled up)
	mylogger.debug("write_sim_dif_matrix(headers, %s, matrix_dict, %s)" %(outname, matrix_method))
	mylogger.info("writing similarity and difference matrices")
	outfilesim = open("%s.%s.sim" %(outname, matrix_method), "w")
	outfilediff = open("%s.%s.diff" %(outname, matrix_method), "w")
	firstline = ""
	
	for org1 in headers:
		current_sim_line, current_diff_line = org1 , org1
		current_line = org1
		for org2 in headers:
			if firstline != None:
				firstline += "\t%s" % org2
			current_sim_line += "\t%s" % matrix_dict[org1]["similarity_dict"][org2]
			current_diff_line += "\t%s" % matrix_dict[org1]["difference_dict"][org2]
#			print current_diff_line
		if firstline != None:
			outfilesim.write(firstline + "\n")
			outfilediff.write(firstline + "\n")
			firstline = None
		outfilesim.write(current_sim_line + "\n")
		outfilediff.write(current_diff_line + "\n")
		
	return outfilesim.name, outfilediff.name

def write_binaryalignment_fasta(outputfilename, columns):
	outfile = open(outputfilename + ".fas", "w")
	if not args.no_verbose:
		mylogger.info("producing alignment file '%s' in fasta format" % outfile.name)
		
	for c in columns:
		outfile.write(">%s\n" % c[0])
		for line in c[1]:
			outfile.write(str(line))
		outfile.write("\n")
		
	outfile.close()
	return outfile.name

def write_binaryalignment_phylip(outputfilename, columns):
	mylogger.debug("write_binaryalignment_phylip(%s, columns)" % outputfilename)
	outfile = open(outputfilename + ".relaxed.phy","w")
	if not args.no_verbose:
		mylogger.info("producing alignment file '%s' in relaxed sequential phylip format" % outfile.name)
		
	outfile.write(str(len(columns)) + " " + str(len(columns[0][1])))
	for c in columns:
		linestring = "\n%s " % c[0]
		for line in c[1]:
			linestring += str(line)
		#print linestring
		outfile.write(linestring)
		
	outfile.close()
	mylogger.info("writing: " + outfile.name)
	return outfile.name

def call_raxml_rapidbs(alignmentfile, outputfilename, seed, parameters): #parameters should be a dictionary (This dictionary thing was introduced, so that the script can be more easily adapted to accept custom commandline-parameters for raxml by the user)
	mylogger.debug("call_raxml_rapidbs(alignmentfile, %s, %s, %s)" % outputfilename, seed, parameters)
	
	mylogger.info("Calculating phylogenies: 'rapid bootstrap analyses and search for best-scoring ML Tree in one run' using raxmlHPC")
	#try:
	
	outname = "GENECONTENT_rapidBS%s_minfreq%s_%s_final_tree" %(args.bootstraps, minfreq, time.strftime("%Y%m%d%H%M%S"))
	raxml_cline = RaxmlCommandline(raxml_prog, sequences = alignmentfile, algorithm = "a", model = "BINGAMMA", name = outname, parsimony_seed = args.seed, rapid_bootstrap_seed = args.seed, num_replicates = args.bootstraps, threads = nr_threads)
	mylogger.info("-->" + str(raxml_cline))
	raxml_cline()
	mylogger.info("-->SUCCESS")
	#the resultfiles will be: "RAxML_bipartitions.rapidBS_final_tree" and "RAxML_bipartitionsBranchLabels.rapidBS_final_tree"
	#Labels on nodes or branches, respectively
	outputfiles = ["RAxML_bipartitions."+outname, "RAxML_bipartitionsBranchLabels."+outname]
	
	#except Exception as e:
		#note for future versions: create a custom exception "raxml_except" instead of using this construct over and over!
	#	mylogger.info("-->FAILURE")
	#	mylogger.warning("rapid bootstrap analyses and search for best-scoring ML tree using raxmlHPC has failed!\n\t" + str(e))
	#	return None

	return outputfiles

def call_raxml_bs(alignmentfile, outputfilename, seed, parameters):
	mylogger.debug("call_raxmlbs(alignmentfile, %s, %s, %s)" % outputfilename, args.seed, parameters)
	mylogger.info("Calculating phylogenies: Thorough bootstrap analyses with raxml")
	
	mylogger.info("\tDetermining best ML tree of 20 raxmlHPC runs")
	#try:
	raxml_cline = RaxmlCommandline(raxml_prog, model = "BINGAMMA", name = "best_delme_tempfile", parsimony_seed = args.seed, num_replicates = 20, sequences = alignmentfile, threads = nr_threads)
	mylogger.info("\t-->"+str(raxml_cline))
	raxml_cline()
	#the resultfile will be :"RAxML_bestTree.best_delme_tempfile"
	mylogger.info("\t-->SUCCESS")
	#except Exception as e:
	#	mylogger.info("\t--FAILURE")
	#	mylogger.warning("WARNING: thorough bootstrap analyses using raxmlHPC failed!\n\t"+str(e))
	#	return None
	
	mylogger.info("\tDoing bootstrap analyses with %s runs using raxmlHPC" % args.bootstraps)
	#try:
	raxml_cline = RaxmlCommandline(raxml_prog, model = "BINGAMMA", sequences = alignmentfile, name = "boot_delme_tempfile", parsimony_seed = args.seed, bootstrap_seed = args.seed, num_replicates = args.bootstraps, threads = nr_threads)
	mylogger.info("\t-->" + str(raxml_cline))
	raxml_cline()
	#the resultfile will be: "RAxML_bootstrap.boot_delme_tempfile"
	mylogger.info("\t-->SUCCESS")
	#except Exception as e:
	#	mylogger.info("\t-->FAILURE")
	#	mylogger.warning("WARNING: thorough bootstrap analyses using raxmlHPC failed!\n\t" + str(e))
	#	return None
	
	mylogger.info("\tDrawing bipartitions of bootstrap trees onto best ML tree using raxmlHPC")
	#try:
	outname = "GENECONTENT_BS%s_minfreq%s_%s_final_tree" %(args.bootstraps, args.min_freq, time.strftime("%Y%m%d%H%M%S"))
	raxml_cline = RaxmlCommandline(raxml_prog, model = "BINGAMMA", parsimony_seed = args.seed, algorithm = "b", starting_tree = "RAxML_bestTree.best_delme_tempfile", bipartition_filename="RAxML_bootstrap.boot_delme_tempfile", name=outname)
	mylogger.info("\t-->"+str(raxml_cline))
	raxml_cline()
	#The resultfiles will be: RAxML_bipartitions.final_tree" and "RAxML_bipartitionsBranchLabels.final_tree"
	outputfiles=["RAxML_bipartitions."+outname,"RAxML_bipartitionsBranchLabels."+outname]
	mylogger.info("\t-->SUCCESS")
	#except Exception as e:
	#	mylogger.info("\t-->FAILURE")
	#	mylogger.warning("WARNING: thorough bootstrap analyses using raxmlHPC failed!\n\t"+str(e))
	#	return None
	
	return outputfiles

def call_raxml_nobs(alignmentfile, outputfilename, seed, parameters):
	mylogger.debug("call_raxml_nobs(alignmentfile, %s, %s, %s)" % outputfilename, args.seed, parameters)
	nr_threads=4
	if "-T" in parameters:
		nr_threads=parameters["-T"]
	#try:
	outname="GENECONTENT_raxml_minfreq%s_%s_final_tree" %(args.min_freq, time.strftime("%Y%m%d%H%M%S"))
	mylogger.info("Calculating phylogeny: Determining best ML tree of 20 raxmlHPC runs")
	raxml_cline=RaxmlCommandline(raxml_prog, sequences = alignmentfile, model = "BINGAMMA", name = outname,  parsimony_seed = args.seed, num_replicates = 20, threads = nr_threads)
	mylogger.info("\t-->" + str(raxml_cline))
	raxml_cline()
	#the resultfile will be :"RAxML_bestTree.final_tree"
	outputfiles = ["RAxML_bestTree." + outname]
	mylogger.info("\tSUCCESS")
	print "deleting temporary files"
	for delfile in os.listdir("."):
		if delfile.startswith("RAxML_") and "." + outname + ".RUN." in delfile:
			#print "deleting temp-file: "+delfile
			os.remove(delfile)
	#except Exception as e:
	#	mylogger.info("\t-->FAILURE")
	#	mylogger.warning("WARNING: searching for best ML tree using raxmlHPC failed!\n\t"+str(e))
	#	return None
	return outputfiles

def calc_timediff(timea, timeb):
	timediff = timeb- timea
	m, s = divmod(timediff, 60)
	h, m = divmod(m, 60)
	return "%d:%02d:%02d" % (h, m, s)

mylogger = setlogger(logfilename)

def main():
	alignment_files = []
	tree_files = []
	mylogger.info("\n ==PO_2_GENECONTENT.py %s by John Vollmers==\n" % version)
	mylogger.info(logtime)
	
	try:
		checkargs()
		headers, columns, OG_number, column_indices = read_PO_file(args.po_file)
		
		if args.min_freq > 1:
			columns, column_indices = filter_OGs_by_freq(args.min_freq, columns, column_indices)
		if args.outfasta:
			alignment_files.append(write_binaryalignment_fasta(out_file, columns))
		
		if args.outdiff != "none":
			mylogger.info("\nGenerating similarity/distance matrices using %s" % args.outdiff)
			if args.outdiff == "simple":
				matrix_dict = calculate_sim_matrix_simple(columns)
			else:
				matrix_dict = calculate_sim_matrix_professional(columns, args.outdiff)
				
			alignment_files.extend(write_sim_dif_matrix(headers, out_file, matrix_dict, args.outdiff))
			matrix_obj = matrix_dict_to_matrix_obj(matrix_dict, headers)
			
			if args.tree_method in ["nj", "nj_bs"]:
				mylogger.info("inferring Neigbor Joining base tree")
				starttime = time.time()
				main_tree = calculate_NJ_tree(matrix_obj)
				endtime = time.time()
				mylogger.debug(" -->inferring Neighbor base tree finished after %s (h:m:s)" % calc_timediff(starttime, endtime))
				if args.tree_method == "nj_bs" and args.bootstraps > 1:
					mylogger.info("inferring Neigbor Joining bootstrapped tree")
					starttime = time.time()
					main_tree = calculate_bootstrapped_NJ_tree(main_tree, create_bootstrap_permutations(columns, args.outdiff, headers))
					endtime = time.time()
					mylogger.debug(" -->inferring Neighbor Joining bootstrapped tree finished after %s (h:m:s)" % calc_timediff(starttime, endtime))
				main_tree = remove_internal_tree_labels(main_tree) #remove internal noda labels which are not confidence values
				tree_files.append(write_nj_tree(out_file, args.outdiff, args.bootstraps, main_tree))
				if not args.no_verbose:
					print "\nascii representation of your NJ tree (not showing confidence values):\n"
					Bio.Phylo.draw_ascii(main_tree)
					print "plain text version of your NJ tree object (WITH confidence values) is included in the log file"
				lstrings.append("\n" + "-" * 50 + "\nplain text version of your NJ tree (with confidence values):\n" + str(main_tree) + "\n" + "-" * 50)
				
		if args.outmatrix:
			alignment_files.append(write_binary_matrix(out_file, columns, column_indices))
		alignment_files.append(write_binaryalignment_phylip(out_file, columns))
		
		if args.tree_method == "raxml_rapidbs":
			tree_files.extend(call_raxml_rapidbs(alignment_files[-1], out_file, args.seed, {"-N":args.bootstraps, "-T":args.ncpus}))
		elif args.tree_method == "raxml_bs":
			tree_files.extend(call_raxml_bs(alignment_files[-1], out_file, args.seed, {"-N":args.bootstraps, "-T":args.ncpus}))
		elif args.tree_method == "raxml":
			print "number of alignmentfiles: %s" % len(alignment_files)
			tree_files.extend(call_raxml_nobs(alignment_files[-1], out_file, args.seed, {"-T":args.ncpus}))
			
		if not args.no_verbose:
			print "\n================================\FINISHED!"
		mylogger.info("Created the following alignment-files:\n\t-" + "\n\t-".join(alignment_files))
		mylogger.info("Created the following tree-files:\n\t-" + "\n\t-".join(tree_files))
		
	except IOError as e:
		for frame in traceback.extract_tb(sys.exc_info()[2]):
			fname,lineno,fn,text = frame
			#print "Error in %s on line %d :> %text" % (fname, lineno, text)
		mylogger.error("Error in %s on line %d :> %s" % (fname, lineno, text))
		mylogger.error(str(e))
		mylogger.error("One or more Inputfiles could not be found")
		
	except OSError as e:
		for frame in traceback.extract_tb(sys.exc_info()[2]):
			fname,lineno,fn,text = frame
			#print "Error in %s on line %d :> %text" % (fname, lineno, text)
		mylogger.error("Error in %s on line %d :> %s" % (fname, lineno, text))
		mylogger.error(str(e))
		mylogger.error("One or more external dependancies for the specified workflow could not be found")
		
	except Exception as e:
		for frame in traceback.extract_tb(sys.exc_info()[2]):
			fname,lineno,fn,text = frame
			#print "Error in %s on line %d :> %text" % (fname, lineno, text)
		mylogger.error("Error in %s, %s, on line %d :> %s" % (fname, fn, lineno, text))
		mylogger.error(str(e))
		
	finally:
		mylogger.info("\ncleaning up...")
		for delfile in os.listdir("."):
			if "delme_tempfile" in delfile:
				os.remove(delfile)
		print "See " + logfilename + " for details\n"

main()
