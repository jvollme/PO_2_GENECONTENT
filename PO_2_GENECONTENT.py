#!/usr/bin/python
#created 21.06.2014 by John Vollmers
#from Bio.Phylo.TreeConstruction import * #for using the directly biopython-implemented Distance(UPGMA + NJ) and parsimony methods. Not done/necessary yet
import os, sys, argparse, time, multiprocessing, random
from Bio import AlignIO, SeqIO
from subprocess import call
from Bio.Phylo.Applications import RaxmlCommandline, PhymlCommandline
myparser=argparse.ArgumentParser(description="\n==PO_2_GENECONTENT.py v1.01 by John Vollmers==\nCreates discrete binary character matrices and phylogenetic trees based on the gene content (presence/absence of orthologues) in comparison organisms.\nThis script is supposed to be part of a pipeline consisting of:\n\tA.)Conversion of Genbank/Embl-Files to ANNOTATED(!) Fastas using CDS_extractor.pl by Andreas Leimbach\n\tB.)Calculation of orthologs and paralogs using proteinortho5 (WITH the '-single' and '-self' arguments!)\n\tC.)The creation of discrete binary character marices (and optionally phylogenetic trees) based on:\n\t\t-the fasta sequences of step A\n\t\t-the proteinortho5-results from step B\n", formatter_class=argparse.RawTextHelpFormatter)
myparser.add_argument("-po", "--proteinortho", action="store", dest="po_resultfile", help="(String) file with proteinortho5 results", required=True)
myparser.add_argument("-s", "--silent", action="store_true", dest="no_verbose", help="non-verbose mode")
myparser.add_argument("-o", "--out", action="store", dest="out_file", default="output", help="Basename for result-files\nDefault='output'\nWill create a phylip-file of the aligned discrete binary character data by default")
myparser.add_argument("-ofas", "--out_fasta", action="store_true", dest="outfasta", help="Create a file with the aligned discrete binary character data in Fasta format (in addition to the default phylip file)")
myparser.add_argument("-omat", "--out_matrix", action="store_true", dest="outmatrix", help="Create a file with the aligned discrete character data in tabular format (in addition to the default phylip file)")
myparser.add_argument("-t", "--threads", action="store", dest="nthreads", type=int, default=4, help="Number of cpus to use\nDefault=4 (or maximum number of available cores, if less than 4 cores available)")
myparser.add_argument("-sd", "--seed", action="store", dest="seed_nr", type=int, default=0, help="Integer to provide as seed for RAxML or PhyML\n0=seed generated randomly\nDefault=random seed")
myparser.add_argument("--raxml", action="store_true", dest="raxml", default=False, help="Generate ML phylogenetic tree without bootstrap values using RAxML and \"new rapid hill climbing\" and substitution model: \"BINGAMMA\"\nDefault=don't use")
myparser.add_argument("--raxml_bs", action="store_true", dest="raxml_BS", default=False, help="Generate ML phylogenetic tree with thorough bootstrap analyses and search for best ML tree using RAxML and \"new rapid hill climbing\" and substitution model: \"BINGAMMA\"'\nDefault: don't use")
myparser.add_argument("--raxml_rapidbs", action="store_true", dest="raxml_rapidBS", default=False, help="Generate ML phylogenetic tree with rapid bootstrap analyses and search for best ML tree in one run using RAxML and substitution model: \"BINGAMMA\"\nDefault: don't use")
myparser.add_argument("--pars", action="store_true", dest="fpars", default=False, help="Generate Parsimony phylogenetic tree using fpars (the EMBOSS implementation of the PHYLIP tool)\nDefault: don't use")
myparser.add_argument("--mix", action="store_true", dest="fmix", default=False, help="Generate Parsimony phylogenetic tree using fmix with wagner parsimony (the EMBOSS implementation of the PHYLIP tool)\nDefault=don't use")
myparser.add_argument("--dollop", action="store_true", dest="fdollop", default=False, help="Generate Parsimony phyogenetic trees using fdollop and Dollo Parsimony (the EMBOSS implementation of the PHYLIP tool)\n Default=don't use")
myparser.add_argument("-bs", "--bootstraps", action="store", dest="bootstraps", type=int, default=1000, help="Number of times to resample for bootstrapping\nonly makes sense in combination with '--raxml_bs' or '--raxml_rapidbs'\nDefault=1000")
#implement bootsrapped versions of pars, mix and dollop soon
args=myparser.parse_args()

version="v0.1"
OG_number=0
available_cores=multiprocessing.cpu_count() #counts how many cores are available, to check if the user-argument for threads can be fulfilled
wstrings, estrings, lstrings=[], [], [] #warning, error and log messages respectively
out_file, PO_file=args.out_file+"_"+time.strftime("%Y%m%d%H%M%S"), args.po_resultfile
nthreads, seed, bootstraps=args.nthreads, args.seed_nr, args.bootstraps

logfilename=out_file+"_PO_2_GENECONTENT_"+time.strftime("%Y%m%d%H%M%S")+".log"
verbose=True
docontinue=True


def randomnumber():
	#returns a random integer to use as seed for rxml and pyml
    random.seed()
    return random.randint(1,2000)

def checkargs(args):
	global verbose, nthreads, seed, bootstraps
	global nthreads
	if args.no_verbose:
		verbose=False
	if not os.path.exists(PO_file) or not os.path.isfile(PO_file):
		datsANerror("ERROR: cannot find proteinortho-resultfile: "+PO_file)
	if available_cores<nthreads:
		datsAwarning("WARNING: less than "+str(nthreads)+" cores/threads available! Setting nthreads to "+str(available_cores))
		nthreads=available_cores
		
	if seed == 0:
		seed=randomnumber()
		datsAlogmessage("Setting seed randomly to "+str(seed))
	else:
		datsAlogmessage("Using "+str(seed)+" as seed")
	if args.raxml or args.raxml_BS or args.raxml_rapidBS:
		if which("raxmlHPC")==None:
			datsANerror("ERROR: raxmlHPC not found in any directory within $PATH")
		else:
			try:
				checkraxml_cline=RaxmlCommandline(version=True)
				versiontext=checkraxml_cline()[0]
				startv=versiontext.find("version ")+len("version ")
				endv=versiontext[startv:].find(" ")
				version=versiontext[startv:startv+endv].split(".")
				if int(version[0])<7 or (int(version[0])==7 and int(version[1])<3) or (int(version[0])==7 and int(version[1])==3 and int(version[2])<5):
					datsAwarning("Warning: This script was devised for and tested with RAxML v7.3.5. Your version is v"+".".join(version)+" !\n\tThis may very well still work, but if it doesn't it's YOUR fault!")
			except:
				datsAwarning("Warning: This script was devised for and tested with RAxML v7.3.5.\n\tNot sure which version of RAxML you're using, but it sure as hell isn't v7.3.5!\n\tThis may very well still work, but if it doesn't it's YOUR fault!")	
	#Add checks for dollop, mix and pars later
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
		print "\n"+wstring+"\n"
	
def datsANerror(estring):
	global estrings
	global docontinue
	estrings.append(estring)
	docontinue=False
	if verbose:
		print >>sys.stderr, "\n"+estring+"\n"

def datsAlogmessage(lstring):
	global lstrings
	lstrings.append(lstring)
	if verbose:
		print lstring

def read_PO_file(filename):
	open_PO_file=open(filename, 'r')
	firstline=True
	org_index=3 #default index of beginning of organism-result-columns in proteinortho4+
	headers=[]
	for line in open_PO_file:
		if firstline:
			if line.startswith("#species"):
				datsAwarning("WARNING! It seems you are using results from Proteinortho4. It's recommended to use Proteinortho5!\n\t(However, this script SHOULD still work)")
			elif line.startswith("# Species"):
				if verbose:
					print "\nProteinortho-results are based on Proteinortho5 or later. Good."
			else:
				datsAwarning("WARNING! Cannot clearly recognize Format of Proteinortho-results! This may produce erroneous results!\n\t For best results use Proteinortho5!") 
			headers=line.rstrip().split("\t")[org_index:]
			columns=[[h.rstrip(),[]] for h in headers]
			firstline=False
		elif not line.startswith("#"):
			zeilentokens=line.rstrip().split("\t")[org_index:]
			if len(zeilentokens)!=len(headers):
				datsANerror("ERROR: Different numbers of columns in header_line and Locus_tag_lines in "+filename)
				break
			for h in range(len(headers)):
				if zeilentokens[h]=="*":
					columns[h][1].append(0) #if no homolog =0
				else:
					columns[h][1].append(1) #if homolog =1
					#print zeilentokens
	open_PO_file.close()
	#check if everything is ok before continuing:
	for c in columns:
		if len(c[1])!=len(columns[0][1]):
			datsANerror("ERROR: something went wrong while reading proteinortho-results. Not the same number of lines (=OGs) for all Organisms!")
	datsAlogmessage("\nrecognized "+str(len(columns[0]))+" total 'Orthologeous Groups' (OGs) + singletons in the comparison organisms")
	OG_number=len(columns[0])
	return headers, columns, OG_number
	
def write_binary_matrix(outname, columns):
	outfile=open(outname+".tab", "w")
	headerstring="OG"
	for c in columns:
		headerstring+=("\t"+c[0])
	outfile.write(headerstring+"\n")
	for og_nr in range(len(columns[0][1])):
		linestring=str(og_nr)
		for c in columns:
			linestring+="\t"+str(c[1][og_nr])
		outfile.write(linestring+"\n")
	outfile.close()		
	return outfile.name

def write_binaryalignment_fasta(outputfilename, columns):
	outfile=open(outputfilename+".fas", "w")
	if verbose:
		print "producing alignment file '"+outfile.name+"' in fasta format"
	for c in columns:
		outfile.write(">"+c[0]+"\n")
		for line in c[1]:
			outfile.write(str(line))
		outfile.write("\n")
	outfile.close()
	return outfile.name
	
def write_binaryalignment_phylip(outputfilename, columns):
	outfile=open(outputfilename+".relaxed.phy","w")
	if verbose:
		print "producing alignment file '"+outfile.name+"' in relaxed sequential phylip format"
	outfile.write(str(len(columns))+" "+str(len(columns[0][1])))
	for c in columns:
		#print c[0]+c[0]
		linestring="\n"+c[0]+" "
		#print "line"+str(c[1][0])+"line"
		for line in c[1]:
			linestring+=str(line)
		#print linestring
		outfile.write(linestring)
	outfile.close()
	print "writing: "+outfile.name
	return outfile.name

def call_raxml_rapidbs(alignmentfile, outputfilename, seed, parameters): #parameters should be a dictionary (This dictionary thing was introduced, so that the script can be more easily adapted to accept custom commandline-parameters for raxml by the user)
	nr_threads=4
	if "-T" in parameters:
		nr_threads=parameters["-T"]
	bootstraps=1000
	if "-N" in parameters:
		bootstraps=parameters["-N"]
	elif "-#" in parameters:
		bootstraps=parameters["-#"]
	
	datsAlogmessage("Calculating phylogenies: 'rapid bootstrap analyses and search for best-scoring ML Tree in one run' using raxmlHPC")
	try:
		outname="rapidBS_"+time.strftime("%Y%m%d%H%M%S")+"_"+"final_tree"	
		raxml_cline=RaxmlCommandline(sequences=alignmentfile, algorithm="a", model="BINGAMMA", name=outname, parsimony_seed=seed, rapid_bootstrap_seed=seed, num_replicates=bootstraps, threads=nr_threads) 
		datsAlogmessage("-->"+str(raxml_cline))
		raxml_cline()
		datsAlogmessage("-->SUCCESS")
		#the resultfiles will be: "RAxML_bipartitions.rapidBS_final_tree" and "RAxML_bipartitionsBranchLabels.rapidBS_final_tree"
		#Labels on nodes or branches, respectively
		outputfiles=["RAxML_bipartitions."+outname, "RAxML_bipartitionsBranchLabels."+outname]
	except Exception as e:
		#note for future versions: create a custom exception "raxml_except" instead of using this construct over and over!
		datsAlogmessage("-->FAILURE")
		datsAwarning("WARNING: rapid bootstrap analyses and search for best-scoring ML tree using raxmlHPC has failed!\n\t"+str(e))	
		return None
	return outputfiles
				
def call_raxml_bs(alignmentfile, outputfilename, seed, parameters):
	nr_threads=4
	if "-T" in parameters:
		nr_threads=parameters["-T"]
	bootstraps=1000
	if "-N" in parameters:
		bootstraps=parameters["-N"]
	elif "-#" in parameters:
		bootstraps=parameters["-#"]
	datsAlogmessage("Calculating phylogenies: Thorough bootstrap analyses with raxml")
		
	datsAlogmessage("\tDetermining best ML tree of 20 raxmlHPC runs") 
	try:
		raxml_cline=RaxmlCommandline(model="BINGAMMA", name="best_delme_tempfile", parsimony_seed=seed, num_replicates=20, sequences=alignmentfile, threads=nr_threads) 
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
		raxml_cline=RaxmlCommandline(model="BINGAMMA", sequences=alignmentfile, name="boot_delme_tempfile", parsimony_seed=seed, bootstrap_seed=seed, num_replicates=bootstraps, threads=nr_threads)
		datsAlogmessage("\t-->"+str(raxml_cline))
		raxml_cline()
		#the resultfile will be: "RAxML_bootstrap.boot_delme_tempfile"
		datsAlogmessage("\t-->SUCCESS")
	except Exception as e:
		datsAlogmessage("\t-->FAILURE")
		datsAwarning("WARNING: thorough bootstrap analyses using raxmlHPC failed!\n\t"+str(e))
		return None
		
	datsAlogmessage("\tDrawing bipartitions of bootstrap trees onto best ML tree using raxmlHPC")
	try:
		outname="BS"+str(bootstraps)+"_"+time.strftime("%Y%m%d%H%M%S")+"_"+"final_tree"
		raxml_cline=RaxmlCommandline(model="BINGAMMA", parsimony_seed=seed, algorithm="b", starting_tree="RAxML_bestTree.best_delme_tempfile", bipartition_filename="RAxML_bootstrap.boot_delme_tempfile", name=outname)
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
		outname="raxml_"+time.strftime("%Y%m%d%H%M%S")+"_"+"final_tree"
		datsAlogmessage("Calculating phylogeny: Determining best ML tree of 20 raxmlHPC runs")
		raxml_cline=RaxmlCommandline(sequences=alignmentfile, model="BINGAMMA", name=outname,  parsimony_seed=seed, num_replicates=20, threads=nr_threads)
		datsAlogmessage("\t-->"+str(raxml_cline))
		raxml_cline()
		#the resultfile will be :"RAxML_bestTree.final_tree"
		outputfiles=["RAxML_bestTree."+outname]
		datsAlogmessage("\tSUCCESS")
		print "deleting temporary files"
		for delfile in os.listdir("."):
			if delfile.startswith("RAxML_") and "."+outname+".RUN." in delfile:
				#print "deleting temp-file: "+delfile
				os.remove(delfile)
	except Exception as e:
		datsAlogmessage("\t-->FAILURE")
		datsAwarning("WARNING: searching for best ML tree using raxmlHPC failed!\n\t"+str(e))
		return None
	return outputfiles

def relaxed_2_strict_phylip(input_file, output_file):
	datsAlogmessage("converting relaxed phylip format ('"+input_file+"') to strict sequential phylip format('"+output_file+"')")
	datsAlogmessage("CAREFUL: Names will be shortened to 10 characters! This may affect branch labeling!")
	try:
		relaxed_file=open(input_file, "r")
		relaxed_align=AlignIO.read(relaxed_file, "phylip-relaxed")
		relaxed_file.close()
		strict_file=open(output_file, "w")
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
		os.rename("strict.treefile",out_file+"_fpars.treefile")
		return [out_file+"_fpars.log", out_file+"_fpars.treefile"]
	except Exception as e:
		datsAlogmessage("-->FAILURE")
		datsAwarning("WARNING: Discrete character parsimony using fpars failed!\n\t"+str(e))
		return None

def call_fmix(alignmentfile, out_file, seed):
	#input should be a temporal "strict" phylip file, NOT the relaxed phylip file that is normally generated!
	#important: first you have to generate a temporal strict phylip file with 10-character placeholders for each name (stored as a dict).
	#use this temporal input file as input for fmix
	#then read in the outputfiles of fpars and translate the names back into the original names (->convert back to relaxed phylip file)
	#but do this outside of this method in order to enable to use the seqboot and condense tools also
	#UPDATE: it appears the use of discrete (binary) characters is NOT possible in Fseqboot.
	datsAlogmessage("calculating phylogeny: Discrete character parsimony using fmix with wagner parsimony (EMBOSS implementation of the PHYLIP package)")
	fpmix_args=['-infile '+alignmentfile, '-outfile '+out_file+"_fmix.log", "-method w", '-auto', "-seed "+str(seed)]
	fmix_command=["fmix"]+fmix_args
	datsAlogmessage("-->"+" ".join(fpars_args))
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
	# NOPE!
	pass

def write_logfile(logfilename):
	logfile=open(logfilename,'w')
	logfile.write("PO_2_GENECONTENT.py logfile")
	logfile.write("\n"+ time.strftime("%c"))
	logfile.write("\nMessages\n:")
	for l in lstrings:
		logfile.write(l+"\n")
	logfile.write("\n-----------\nWarnings:\n")
	for w in wstrings:
		logfile.write(w+"\n")
	logfile.write("\n------------\nErrors:\n")
	for e in estrings:
		logfile.write(e+"\n")
	logfile.close()

def main():
	print "MOINMOIN"
	alignment_files=[]
	tree_files=[]
	if verbose:
		print "\n ==PO_2_GENECONTENT.py "+version+" by John Vollmers==\n"
	checkargs(args)
	if docontinue:
		try:
			headers, columns, OG_number = read_PO_file(PO_file)
			if args.outfasta and docontinue:
				alignment_files.append(write_binaryalignment_fasta(out_file, columns))
			if args.outmatrix and docontinue:
				alignment_files.append(write_binary_matrix(out_file, columns))
			if docontinue:
				alignment_files.append(write_binaryalignment_phylip(out_file, columns))
			if args.raxml_rapidBS and docontinue:
				tree_files.extend(call_raxml_rapidbs(alignment_files[-1], out_file, seed, {"-N":bootstraps, "-T":nthreads}))
			if args.raxml_BS and docontinue:
				tree_files.extend(call_raxml_bs(alignment_files[-1], out_file, seed, {"-N":bootstraps, "-T":nthreads}))
			if args.raxml and docontinue:
				print "number of alignmentfiles: "+str(len(alignment_files))
				tree_files.extend(call_raxml_nobs(alignment_files[-1], out_file, seed, {"-T":nthreads}))
			if args.fpars or args.fmix or args.fdollop and docontinue:
				alignment_files.append(relaxed_2_strict_phylip(alignment_files[-1], alignment_files[-1].replace(".relaxed.", ".strict.")))
				if alignment_files[-1]!=None:
					if args.fpars:
						treefiles.extend(call_fpars(alignment_files[-1], outfile, seed))
					if args.fdollop:
						treefiles.extend(call_fdollop(alignment_files[-1], outfile, seed))
					if args.fmix:
						treefiles.extend(call_fmix(alignment_files[-1], outfile,seed))
			if docontinue:
				if verbose:
					print "================================\FINISHED!"	
				datsAlogmessage("Created the following alignment-files:"+"\n\t-".join(alignment_files))
				datsAlogmessage("Created the following tree-files:"+"\n\t-".join(tree_files))					
		except Exception as e:
			datsAlogmessage(str(e))
		finally:
			datsAlogmessage("cleaning up...")
			for delfile in os.listdir("."):
				if "delme_tempfile" in delfile:
					os.remove(delfile)
					
	if len(wstrings)>0 or len(estrings)>0:
		if len(wstrings)>0:
			print "\nThere were Warnings!"
		if len(estrings)>0:
			print "\nThere were Errors!"
	print "See "+logfilename+" for details\n"
	write_logfile(logfilename)

print "moin"
main()
