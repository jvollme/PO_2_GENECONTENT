usage: PO_2_GENECONTENT.py [-h] -po PO_RESULTFILE [-s] [-o OUT_FILE] [-ofas]
                           [-omat] [-t NTHREADS] [-sd SEED_NR]
                           [-mt {raxml,raxml_bs,raxml_rapidbs,none}]
                           [-tbp TREEBUILDER_PATH] [-bs BOOTSTRAPS]
                           [--min_freq MIN_FREQ]

==PO_2_GENECONTENT.py v1.01 by John Vollmers==
Creates discrete binary character matrices and phylogenetic trees based on the gene content (presence/absence of orthologues) in comparison organisms.
This script is supposed to be part of a pipeline consisting of:
	A.)Conversion of Genbank/Embl-Files to ANNOTATED(!) Fastas using CDS_extractor.pl by Andreas Leimbach
	B.)Calculation of orthologs and paralogs using proteinortho5 (WITH the '-single' and '-self' arguments!)
	C.)The creation of discrete binary character marices (and optionally phylogenetic trees) based on:
		-the fasta sequences of step A
		-the proteinortho5-results from step B

optional arguments:
  -h, --help            show this help message and exit
  -po PO_RESULTFILE, --proteinortho PO_RESULTFILE
                        (String) file with proteinortho5 results
  -s, --silent          non-verbose mode
  -o OUT_FILE, --out OUT_FILE
                        Basename for result-files
                        Default='output'
                        Will create a phylip-file of the aligned discrete binary character data by default
  -ofas, --out_fasta    Create a file with the aligned discrete binary character data in Fasta format (in addition to the default phylip file)
  -omat, --out_matrix   Create a file with the aligned discrete character data in tabular format (in addition to the default phylip file)
  -t NTHREADS, --threads NTHREADS
                        Number of cpus to use
                        Default=4 (or maximum number of available cores, if less than 4 cores available)
  -sd SEED_NR, --seed SEED_NR
                        Integer to provide as seed for RAxML or PhyML
                        0=seed generated randomly
                        Default=random seed
  -mt {raxml,raxml_bs,raxml_rapidbs,none}, --make_tree {raxml,raxml_bs,raxml_rapidbs,none}
                        Generate ML phylogenetic trees using RAxML with "new rapid hill climbing" and substitution model: "BINGAMMA"
                        	choices:	"raxml": single tree without bootstraps (using new rapid hill climbing)
                        		raxml_bs: thorough bootstrap analyses and search for best ML tree
                        		raxml_rapidbs: rapid bootstrap analyses and search for best ML tree in one run
                        		none
                        Default=none
  -tbp TREEBUILDER_PATH, --tree_builder_path TREEBUILDER_PATH
                        Path to treebuilder (currently only raxml supported) if not listed in $PATH
  -bs BOOTSTRAPS, --bootstraps BOOTSTRAPS
                        Number of times to resample for bootstrapping
                        only makes sense in combination with '-mt raxml_bs' or '-mt raxml_rapidbs'
                        Default=1000
  --min_freq MIN_FREQ   Minimum frequency for Orthologeous Groups across all comparison-organisms for Group to be considered in Analyses
                         Default=2 (exclude all singletons)
