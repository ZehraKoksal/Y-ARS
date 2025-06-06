from cyvcf2 import VCF, Writer
import argparse
import os
import glob
import pandas as pd
import pysam
import concurrent.futures
from multiprocessing import Manager

parser = argparse.ArgumentParser()
parser.add_argument('-chromosome', type=str, required=True, help='Label of Y chromosome in the vcf file. Can be "chrY", "Y" or e.g. "NC_000024.9" for the GRCh37')
parser.add_argument('-reference', type=str, choices=['T2T', 'GRCh37','GRCh38'], required=True,
                    help='Specify the reference genome. Choices: "T2T", "GRCh37", "GRCh38"')
parser.add_argument('-input_single_vcf', type=str, required=False,
                    help='Folder containing .vcf files')
parser.add_argument('-multi_sample_vcf', type=str, required=False,
                    help='Path to multisample vcf file in .vcf file.')
parser.add_argument('-output_loci_dict', action='store_true', required=False,
                    help='If specified, will print all loci and ancestral/derived alleles if available in Y-ARS for single sample vcfs.')
parser.add_argument('-output', type=str, required=False,
                    help='Define folder to save output files. If not defined, uses directory of polaryzer script.')


def revert_delta_encoding(delta_encoded_list):
    original_positions = [delta_encoded_list[0]]  # Start with the first value
    # Accumulate deltas to get original positions
    for delta in delta_encoded_list[1:]:
        original_positions.append(original_positions[-1] + delta)
    return original_positions
        
args = parser.parse_args()

chromosome = args.chromosome
input_folder = args.input_single_vcf
output_dict = []
output_dict_anc = []
output_dict_der = []

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


if args.output != None:
    if not os.path.exists(args.output):
        print(f"output directory {args.output} does not exist and will be created.")
        # Create the directory if it doesn't exist
        os.makedirs(args.output)
        print(f"Directory '{args.output}' was created.")
    else:
        print(f"Directory '{args.output}' exists for storing polaryzer output.")



df_exploded.columns=["anc","pos"]
print(df_exploded)


def process_vcf_file(vcf_file, chromosome, position_to_ancestral, output_dict, output_dict_anc, output_dict_der, args):
    vcf_in = pysam.VariantFile(vcf_file, 'r')
    modified_filename = vcf_file.replace(".vcf", "")
    out_vcf_path = f"{modified_filename}_polarized.vcf"
    if args.output != None:
        modified_filename = modified_filename.rsplit('/', 1)[-1]
        out_vcf_path = f"{args.output}/{modified_filename}_polarized.vcf"
    print(out_vcf_path)
    vcf_out = pysam.VariantFile(out_vcf_path, 'w', header=vcf_in.header)
    # Write records in batches to avoid memory overload
    for record in vcf_in:
        if record.chrom == chromosome:
            position = record.pos
            if position in position_to_ancestral:
                ancestral = position_to_ancestral[position]
                output_dict.append(position)
                output_dict_anc.append(ancestral)
                
                # Iterate over samples and update the record
                for sample in record.samples:
                    genotype = record.samples[sample]['GT']
                    print(genotype)
                    if genotype is not None:
                        alleles = [record.ref] + list(record.alts)
                        allele = [alleles[g] if g is not None else '.' for g in genotype]
                        if allele[0] == ".":
                            record.id = f"{record.id};Missing data"
                        elif allele[0] != ancestral:
                            record.id = f"{record.id};Derived"
                            output_dict_der.append(allele[0])
                        elif allele[0] == ancestral:
                            record.id = f"{record.id};Ancestral"
                
                vcf_out.write(record)
            else:
                vcf_out.write(record)
                print(f"Position {position} is not in Reference Sequence and will therefore not be annotated with polarization. Please check if the position is correct. If you are using GRCh37 or GRCh38-aligned variants, consider lifting over to T2T first and then repeating polaryzer.")

    vcf_in.close()
    vcf_out.close()
    print(f"Finished processing {vcf_file}")


def main(args, input_folder, chromosome, df_exploded):
    if args.input_single_vcf:
        if args.input_single_vcf.endswith(".vcf"):
            # print("Single vcf file recognized.")
            vcf_files = []
            vcf_files.append(args.input_single_vcf)
            print(vcf_files)
        else:
            # Check if the folder exists
            if not os.path.isdir(input_folder):
                raise ValueError(f"The folder '{input_folder}' does not exist.")

            # Use glob to find all .vcf files in the folder
            vcf_files = glob.glob(os.path.join(input_folder, '*.vcf'))
            print(vcf_files)
            # Check if any .vcf files were found
            if not vcf_files:
                print(f"No .vcf files found in the directory '{input_folder}'.")
            else:
                print(f"Found the following .vcf files in '{input_folder}':")
       
                
        for vcf_file in vcf_files:
            print(vcf_file)
        
        # Prepare the ancestral dictionary for fast lookup
        position_to_ancestral = dict(zip(df_exploded['pos'], df_exploded['anc']))
        
        # Use Manager to handle shared data across processes
        with Manager() as manager:
            output_dict = manager.list()
            output_dict_anc = manager.list()
            output_dict_der = manager.list()

            # Use ProcessPoolExecutor to process files in parallel
            with concurrent.futures.ProcessPoolExecutor() as executor:
                futures = [
                    executor.submit(process_vcf_file, vcf_file, chromosome, position_to_ancestral, output_dict, output_dict_anc, output_dict_der, args)
                    for vcf_file in vcf_files
                ]
                # Wait for all processes to complete
                concurrent.futures.wait(futures)

            # Process the loci dictionary if requested
            if args.output_loci_dict:
                print(f"Successfully generated output file with Y-ARS recognized loci and ancestral/derived alleles (-output_loci_dict).")
                output_loci_dict_df = pd.DataFrame({
                    f"loci_{args.reference}": list(output_dict),
                    "ancestral": list(output_dict_anc),
                    "derived": list(output_dict_der)
                }).drop_duplicates().sort_values(by=f"loci_{args.reference}")

                if args.output != None:
                    output_loci_dict_df.to_csv(f'{args.output}/output_loci_dict.csv', sep="\t", index=False)
                else:
                    output_loci_dict_df.to_csv('output_loci_dict.csv', sep="\t", index=False)


#Single VCF file approach
if args.input_single_vcf:
    # Argument parsing
    if __name__ == "__main__":
        main(args, args.input_single_vcf, args.chromosome, df_exploded)   

elif args.multi_sample_vcf:
    vcf_in = pysam.VariantFile(args.multi_sample_vcf, 'r')
    output = pd.DataFrame()
    sample_names = list(vcf_in.header.samples)
    header_columns = ["chr","pos"] + sample_names + ["anc","der"]
    output = pd.DataFrame(columns = header_columns)
    for record in vcf_in:
        if record.chrom == chromosome:
            position = record.pos
            pos_entry = []
            pos_derived_allele = []
            pos_entry.append(chromosome)
            pos_entry.append(position)
            ancestral = df_exploded[df_exploded["pos"] == position]
            if not ancestral.empty:
                output_dict.append(position)
                ancestral = ancestral["anc"].values[0]  # Extract the value of the second column / ancestral allele
                for sample in record.samples:
                    # Get the genotype (GT) field for the current sample
                    genotype = record.samples[sample]['GT']
                    if genotype is not None:
                        alleles = [record.ref] + list(record.alts)  # List of alleles (ref + alts)
                        allele = [alleles[g] if g is not None else '.' for g in genotype]
                        if allele[0] == ".":
                            pos_entry.append("Missing data")
                        elif allele[0] != ancestral:
                            pos_entry.append("Derived")
                            pos_derived_allele.append(allele[0])
                        elif allele[0] == ancestral:
                            pos_entry.append("Ancestral")
                    else:
                        pos_entry.append(".")
            else:
                print(f"{position} is not on Y chromosome reference. Check if the position is correct.")
                for sample in record.samples:
                    pos_entry.append("-")
                ancestral = "-"
                #for derived allele
                pos_derived_allele.append("-")
            pos_derived_allele = set(pos_derived_allele)
            pos_derived_allele = ''.join(str(item) for item in pos_derived_allele)
            pos_entry.append(ancestral)
            pos_entry.append(pos_derived_allele)
            output.loc[len(output)] = pos_entry
    print(output)
    if args.output != None:
        output.to_csv(f'{args.output}/polaryzer_output.csv', sep="\t", index=False)
    else:
        output.to_csv('polaryzer_output.csv', sep="\t", index=False)
    vcf_in.close()
    print(f"Finished processing {args.multi_sample_vcf}")

else:
    print("No input vcf file was provided! Polaryzer aborted. Please provide path to vcf folder with several singel vcf files (-input_single_vcf) or path to a multi sample vcf file (-multi_sample_vcf).")
print("Finished polarizing")
