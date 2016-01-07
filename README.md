##PO_2_GENECONTENT.py

Automizes the creation of discrete binary character matrices and phylogenetic trees based on the gene content (presence/absence of orthologues) in comparison organisms. Several Options for generating Distance matrices are available and can be automatically clustered (with or without bootstraps) using the neighbor-joining algorithm. If RAxML is installed, maximum likelyhood trees based (with or without bootstraps) on these character matrices can also be created automatically. Any trees that are generated will be output in the newick format.

This script only need a proteinortho result file as sole input. 
However, This script is supposed to be part of a pipeline consisting of:

1. Conversion of Genbank/Embl-Files to **annotated**(!) Fastas using [CDS_extractor.pl][] by Andreas Leimbach

2. Calculation of orthologs and paralogs using proteinortho5 (**with** the *'-single'* and *'-self'* arguments!)

2. The creation of discrete binary character marices (and optionally phylogenetic trees) based on:
    * the [proteinortho5][]-results from step 2.)

####usage: 
 
     PO_2_GENECONTENT.py [-h] -po PO_RESULTFILE [-s] [-o OUT_FILE] [-ofas]
                         [-omat] [-t NTHREADS] [-sd SEED_NR]              
                         [-mt {raxml,raxml_bs,raxml_rapidbs,none}]  
                         [-tbp TREEBUILDER_PATH] [-bs BOOTSTRAPS]  
                         [--min_freq MIN_FREQ]  

####Dependancies:

external tools:
 - [proteinortho5][]
 - [CDS_extractor.pl][] (optional; necessary for completing pipeline with PO_2_MLSA.py and PO_2_ANNOTAION.py)
 - [RaXML][] (optional; necessary for inferring maximum likelihood-based gene content trees) 
 - [EMBOSS package] [] (optional; necessary for infering neighbor-joining-based gene content trees)

NOTE: neighbor-joining-based tree inference was recently switched from biopython implemented methods to the external tool fneighbor (of the PHYLIP package implemented in the EMBOSS suite), because this was much faster for large datasets. This is, however, no longer compatible with all distance matrixes created with SciPy functions (when not using "-odiff none" or "-odiff simple").
Will make this optional in future releases.

python modules:
 - [BioPython][] version 2.6.+
 - [SciPy][] version 0.13.+ or higher (optional; necessary for distance matrix based inferrence of gene content trees)
 
**optional arguments:**
````
  -h, --help            show this help message and exit
  -po PO_FILE, --proteinortho PO_FILE
                        (String) file with proteinortho5 results
  -s, --silent          non-verbose mode
  -o OUT_FILE, --out OUT_FILE
                        Basename for result-files
                        Default='output'
                        Will create a phylip-file of the aligned discrete binary character data by default
  -ofas, --out_fasta    Create a file with the aligned discrete binary character data in Fasta format
                        (in addition to the default phylip file)
  -omat, --out_matrix   Create a file with the aligned discrete character data in tabular format
                        (in addition to the default phylip file)
  -odiff {simple,none,braycurtis,canberra,cityblock,jaccard,euclidean}
                        Create difference/similarity matrixes for distance matrix based phylogeny inference.
                        	- simple : Simple normalized distance matrix
                        	  (distance=Nr shared OGs/total Nr OGs in SMALLER organism)
                        	- braycurtis : Bray-Curtis distance
                        	- canberra : Canberra distance
                        	- cityblock : Manhatten distance
                        	- euclidean : Euclidean distance
                        	- jaccard : Jaccard-Needham dissimilarity
                        Default: "simple"
  -cpu NCPUS, --cpu NCPUS
                        Number of cpus to use
                        Default=4 (or maximum number of available cores, if less than 4 cores available)
  -sd SEED, --seed SEED
                        Integer to provide as seed for RAxML or PhyML
                        0=seed generated randomly
                        Default=random seed
  -mt {raxml,raxml_bs,raxml_rapidbs,nj,nj_bs,none}, --make_tree {raxml,raxml_bs,raxml_rapidbs,nj,nj_bs,none}
                        Generate ML phylogenetic trees using RAxML with "new rapid hill climbing"
                        and substitution model: "BINGAMMA"
                        	-"raxml": single tree without bootstraps (using new rapid hill climbing)
                        	-"raxml_bs": thorough bootstrap analyses and search for best ML tree
                        	-"raxml_rapidbs": rapid bootstrap analyses and search for best ML tree
                        	 in one run
                        	-"nj": Neighbor joining
                        	-"none"
                        Default=none
  -tbp TREEBUILDER_PATH, --tree_builder_path TREEBUILDER_PATH
                        Path to treebuilder (currently only raxml supported) if not listed in $PATH
  -bs BOOTSTRAPS, --bootstraps BOOTSTRAPS
                        Number of times to resample for bootstrapping
                        only bootstrapped trees will be produced, not the permutated data matrices
  --min_freq MIN_FREQ   Minimum frequency for Orthologeous Groups across all comparison-organisms
                        for Group to be considered in Analyses
                        Default=2 (exclude all singletons)
                        Recommended when working with highly fragmented genomes or dubious ORF-finding
  --debug               Log extra info for debugging
````
[proteinortho5]: https://www.bioinf.uni-leipzig.de/Software/proteinortho/
[CDS_extractor.pl]: https://github.com/aleimba/bac-genomics-scripts.git
[raxml]: http://sco.h-its.org/exelixis/web/software/raxml/index.html
[EMBOSS package]: http://emboss.sourceforge.net/
[BioPython]: http://biopython.org/wiki/Main_Page
[SciPy]: http://www.scipy.org/
