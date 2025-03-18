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
Make sure to use a newer version of python to avoid problems with the package pyarrow. We recommend python 3.11.11.

In a (linux) command line you can run the following commands:
```
python3.11 -m venv Y_is_fun
source Y_is_fun/bin/activate
pip install numpy
pip install cyvcf2
pip install pandas
pip install pysam
pip install pyarrow
```
The following versions of the packages are used:
numpy v2.2.4, cyvcf2 v0.31.1, pandas v2.2.3, pysam v0.23.0, pyarrow v19.0.1

Now the environment is prepared to run polaryzer!

### 3) Run polaryzer
Polaryzer can be applied to single sample vcf files and multisample vcf files. The first results in an annotated vcf file with the polarization (ancestral/derived) in the ID column. For the latter, a tab-separated .csv file is generated with samples in columns and loci in rows.

#### a) Single sample vcf input file
When running polaryzer on single sample vcf files, the user needs to specify the exact name of the Y chromosome used in the CHR column of the vcf file using parameter **-chromosome**, defining the reference sequence for alignment/in SNP array among GRCh37, GRCh38 or T2T following **-reference**. Define the path to the folder containing all .vcf files following parameter **-input_single_vcf**.

```
python polarYzer.py -chromosome NC_060948.1 -reference T2T -input_single_vcf vcf_T2T_test/
```
Optionally, the parameter **-output_loci_dict** can be added to obtain a tab-separated csv file containing all loci and the ancestral and derived allele information over all vcf files.
```
python polarYzer.py -chromosome NC_060948.1 -reference T2T -input_single_vcf vcf_T2T_test/ -output_loci_dict
```

The single sample vcf file mode is run in parallel mode to reduce computing time.

#### b) Multi sample vcf input file
When running polaryzer on **one** multiple sample vcf file, the user needs to specify the exact name of the Y chromosome used in the CHR column of the vcf file using parameter **-chromosome**, defining the reference sequence for alignment/in SNP array among GRCh37, GRCh38 or T2T following **-reference**. Define the path to the vcf file following parameter **-multiple_sample_vcf**.
```
python polarYzer.py -chromosome NC_060948.1 -reference T2T -multi_sample_vcf multisample_vcf_t2t.vcf
```

The resulting output file is a tab-separated .csv file with loci being different rows, and the samples different columns. 




#### c) Example files

Different test vcf files from the 1000 Genomes Project for inputs in the **single sample vcf file mode** of the three reference genomes are available in the [subfolder example_vcfs](https://github.com/ZehraKoksal/Y-ARS/tree/main/Polaryzer/example_vcfs)



### 4) Additional Information and Contact
More information on the software are available in [our publication:]()

For reporting bugs, comments or questions, you are welcome to contact zehra.koksal@liu.se

### 5) Referencing

Please cite: 


