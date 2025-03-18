#### POLARYZER
Accesses the Y-ARS to identify the polarization of alleles (i.e. , ancestral versus derived alleles) of variants reported in vcf files.

### 1) Installation
Operating system: Linux

Type in the shell:
```
git clone git@github.com:ZehraKoksal/Y-ARS.git
cd Y-ARS/Polaryzer/
python polarYzer.py -h
```


### 2) Requirements and Commands
Polaryzer is written using python and requires a small set of packages that can be installed in a new virtual environment. 

In a (linux) command line you can run the following commands:
```
python3 -m venv Y_is_fun
source Y_is_fun/bin/activate
pip install numpy
pip install cyvcf2
pip install pandas
pip install pysam
```
Now the environment is prepared to run polaryzer!

### 3) Run polaryzer
Polaryzer can be applied to single sample vcf files and multisample vcf files. The first results in an annotated vcf file with the polarization (ancestral/derived) in the ID column. For the latter, a tab-separated .csv file is generated with samples in columns and loci in rows.

#### a) Single sample vcf input file
When running polaryzer on single sample vcf files, the user needs to specify the exact name of the Y chromosome used in the CHR column of the vcf file using parameter **-chromosome**
```
python polarYzer.py -chromosome NC_060948.1 -reference T2T -input_single_vcf vcf_T2T_test/ -output_loci_dict
```

#### b) Multisample vcf file
```
python polarYzer.py -chromosome NC_060948.1 -reference T2T -multi_sample_vcf multisample_vcf_t2t.vcf
```

SNPtotree generates the phylogenetic tree in two file formats: in a tab-separated csv file and a traditional phyloxml file.

The csv output file is to be read from left to right. Downstream variants are located in the cells below and to the right. In this example variants b, d, e and f are downstream of variant a, and variant c is downstream of variant b.
Variants in parallel branches (sister clades) within a clade are presented in a column: Variants g and h are parallel to each other. Not separable variants based on the available data are presented in one cell divided by a comma, like variants i and j.

<img src="/Images/output_phyltree.png" alt="Input file style" width="700"/>

The phyloxml output file contains annotated branches and nodes and can be viewed in phylogenetic tree visualization tools that support phyloxml format, e.g. the interactive Tree Of Life (iTOL). This tree contains the same information as the csv output file.

<img src="/Images/output_phyloxml.png" alt="Input file style" width="900"/>
<img src="/Images/output2_phyloxml.png" alt="Input file style" width="900"/>

#### c) Additional Output File: Certainty Value File

In a separate csv file, statistical support values for each variant present in the phylogenetic tree are given. The support or certainty values are the fraction of variants in the tree that support that variant's position based on their informative (=upstream/basal/rootward, downstream/terminal, parallel) pairwise relationships out of all remaining variants in the tree.

<img src="/Images/certainty_values_example.png" alt="Input file style" width="400"/>

In this example, the tree location of variant i is supported by informative pairwise relationships of 8 of the remaining 9 variants (8/9 = 0.88888). The ninth variant has no informative pairwise relationship to variant i. Variants with contradictory relationships are not present in the final tree, since they have been filtered out during the tree generation process.

#### d) Optional Output Files:

**metadata_individuals**


In t

<img src="/Images/output_phyltree_metadata.png" alt="Input file style" width="700"/>




### 4) Additional Information and Contact
More information on the software are available in [our publication:]()

For reporting bugs, comments or questions, you are welcome to contact zehra.koksal@liu.se

### 5) Referencing

Please cite: 


