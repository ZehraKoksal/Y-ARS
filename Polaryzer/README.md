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

### 3) Run polaryzer

#### a) Input file
The user is required to provide the path to the input file in the _.csv_ format. The input file contains the ancestral **A** or derived **D** allelic state or missing information **X** for each **_polymorphic variant_** in a tab-separated format. The rows present variants, the columns individuals.
The header row should present the individuals' labels and the first (index) column the variant names.

<img src="/Images/inputfile_snptotree.png" alt="Input file style" width="700"/>

The allelic states "ancestral" and "derived" of the most used model organisms are reported in public repositories. For novel SNPs or for not well investigated organisms without already reported relevant SNP information, the ancestral and derived allelic states have to be identified by the user. The ancestral allele is found in a common ancestor of the group of analysed individuals. Thus, it is helpful to conduct sequence alignments to a common ancestor rather than an arbitrarily selected reference genome, e.g. GRCh38 for humans.


#### b) Main Output Files: Phylogenetic Tree

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


In the metadata output file, the individuals presented in each row correspond to the respective row of the phylogenetic tree (tree layer). The variants in each tree layer were observed in the sequences of the respective metadata output row. In this example, variant a was found in all individuals 1 to 10, and variant b was only found in individuals 2, 3 and 4. Variants that could not be separated into different branches were represented in one tree layer (like variants i and j). In this case, the sequences corresponding to this tree layer (individuals 11 and 12) were each found in at least one of the variants (i and j).

<img src="/Images/output_phyltree_metadata.png" alt="Input file style" width="700"/>


**contradictory_variants**

For certain variants - including those resulting from sequencing errors, recurrent mutations or backmutations - variants with contradictory pairwise hierarchical relationships are ignored for the tree generation, but saved as "contradictory variants". 

**ambiguous_variants**

Based on the pairwise relationships, the final hierarchical order of the variants is established.
For some variants, an explicit position in the tree could not be determined. These variants have ambiguous positions in the tree.


### 4) Additional Information and Contact
More information on the software are available in [our publication:](https://www.mdpi.com/2073-4425/14/10/1837)

For reporting bugs, comments or questions, you are welcome to contact zehra.koksal@sund.ku.dk.

### 5) Referencing

Please cite: 
Köksal, Z.; Børsting, C.; Gusmão, L.; Pereira, V. SNPtotree—Resolving the Phylogeny of SNPs on Non-Recombining DNA. Genes 2023, 14, 1837. https://doi.org/10.3390/genes14101837

