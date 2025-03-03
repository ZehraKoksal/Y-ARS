from cyvcf2 import VCF, Writer
import argparse
import os
import glob
import pandas as pd
import pysam

parser = argparse.ArgumentParser()
parser.add_argument('-chromosome', type=str, required=True, help='Label of Y chromosome in the vcf file. Can be "chrY", "Y" or e.g. "NC_000024.9" for the GRCh37')
parser.add_argument('-reference', type=str, choices=['T2T', 'GRCh37','GRCh38'], required=True,
                    help='Specify the reference genome. Choices: "T2T", "GRCh37", "GRCh38"')
parser.add_argument('-input_vcf', type=str, required=True,
                    help='Folder containing .vcf files')

def revert_delta_encoding(delta_encoded_list):
    original_positions = [delta_encoded_list[0]]  # Start with the first value
    # Accumulate deltas to get original positions
    for delta in delta_encoded_list[1:]:
        original_positions.append(original_positions[-1] + delta)
    return original_positions
        
args = parser.parse_args()

chromosome = args.chromosome
input_folder = args.input_vcf

#Load the dictionaries
if args.reference == "GRCh37":
    grch37 = pd.read_parquet('grch37_dictionary.parquet', engine='pyarrow')
    df_exploded = grch37.explode(1)
elif args.reference == "GRCh38":
    grch38 = pd.read_parquet('grch38_dictionary.parquet', engine='pyarrow')
    df_exploded = grch38.explode(1)
elif args.reference == "T2T":
    t2t = pd.read_parquet('t2t_dictionary_delta_enc.parquet', engine='pyarrow')
    t2t[1] = t2t[1].apply(revert_delta_encoding)
    df_exploded = t2t.explode(1)

# Check if the folder exists
if not os.path.isdir(input_folder):
    raise ValueError(f"The folder '{input_folder}' does not exist.")
# Use glob to find all .vcf files in the folder
vcf_files = glob.glob(os.path.join(input_folder, '*.vcf'))
# Check if any .vcf files were found
if not vcf_files:
    print(f"No .vcf files found in the directory '{input_folder}'.")
else:
    print(f"Found the following .vcf files in '{input_folder}':")
    for vcf_file in vcf_files[:1]:
        print(vcf_file)
        vcf_in = pysam.VariantFile(vcf_file, 'r')
        modified_filename = vcf_file.replace(".vcf", "")
        vcf_out = pysam.VariantFile(f"{modified_filename}_polarized.vcf", 'w', header=vcf_in.header)
        filtered_variants = []
        for record in vcf_in:
            if record.chrom == chromosome:
                position = record.pos
                ancestral = df_exploded[df_exploded[1] == position]
                if not ancestral.empty:
                    ancestral = ancestral[2].values[0]  # Extract the value of the second column
                    print(ancestral)
                    print(record.alts)
                    if record.alts != ".":
                        if ancestral == record.alts[0]:
                            record.id = "Ancestral"
                        else:
                            record.id = "Derived"
                    elif record.alts == ".":
                        if ancestral == record.ref:
                            record.id = "Ancestral"
                        else:
                            record.id = "Derived"
                    vcf_out.write(record)   
                else:
                    vcf_out.write(record) 
        vcf_in.close()
        vcf_out.close()
        print(f"Finished processing {vcf_file}.")

print("Finished polarizing")
