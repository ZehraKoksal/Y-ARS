from cyvcf2 import VCF, Writer
import argparse
import os
import glob
import pandas as pd
import pysam
import concurrent.futures

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

def process_vcf_file(vcf_file, chromosome, df_exploded, output_dict, output_dict_anc, output_dict_der, args):
    vcf_in = pysam.VariantFile(vcf_file, 'r')
    modified_filename = vcf_file.replace(".vcf", "")
    vcf_out = pysam.VariantFile(f"{modified_filename}_polarized.vcf", 'w', header=vcf_in.header)

    for record in vcf_in:
        if record.chrom == chromosome:
            position = record.pos
            ancestral = df_exploded[df_exploded[1] == position]
            if not ancestral.empty:
                output_dict.append(position)
                ancestral = ancestral[2].values[0]  # Extract the value of the second column
                output_dict_anc.append(ancestral)
                # For the output
                if args.output_loci_dict:
                    if ancestral != record.alts:
                        output_dict_der.append(record.alts[0])
                    if ancestral != record.ref:
                        output_dict_der.append(record.ref)
                if record.alts != ".":
                    if ancestral == record.alts[0]:
                        record.id = f"{record.id};Ancestral"
                    else:
                        record.id = f"{record.id};Derived"
                elif record.alts == ".":
                    if ancestral == record.ref:
                        record.id = f"{record.id};Ancestral"
                    else:
                        record.id = f"{record.id};Derived"
                vcf_out.write(record)
            else:
                vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()
    print(f"Finished processing {vcf_file}")

def main(args, input_folder, chromosome, df_exploded):
    if args.input_single_vcf:
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
            for vcf_file in vcf_files:
                print(vcf_file)

            # Prepare the output lists
            output_dict = []
            output_dict_anc = []
            output_dict_der = []

            # Use ThreadPoolExecutor to process files in parallel
            with concurrent.futures.ThreadPoolExecutor() as executor:
                # Pass arguments to process_vcf_file function and start processing in parallel
                futures = [
                    executor.submit(process_vcf_file, vcf_file, chromosome, df_exploded, output_dict, output_dict_anc, output_dict_der, args)
                    for vcf_file in vcf_files
                ]
                # Wait for all threads to complete
                concurrent.futures.wait(futures)

            # Process the loci dictionary if requested
            if args.output_loci_dict:
                print(f"Successfully generated output file with Y-ARS recognized loci and ancestral/derived alleles (-output_loci_dict).")
                output_loci_dict_df = pd.DataFrame()
                output_loci_dict_df[f"loci_{args.reference}"] = output_dict
                output_loci_dict_df[f"ancestral"] = output_dict_anc
                output_loci_dict_df[f"derived"] = output_dict_der
                output_loci_dict_df = output_loci_dict_df.drop_duplicates()
                output_loci_dict_df = output_loci_dict_df.sort_values(by=f"loci_{args.reference}")
                output_loci_dict_df.to_csv("output_loci_dict.csv", sep="\t", index=False)

    
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
            ancestral = df_exploded[df_exploded[1] == position]
            if not ancestral.empty:
                output_dict.append(position)
                ancestral = ancestral[2].values[0]  # Extract the value of the second column
                for sample in record.samples:
                    # Get the genotype (GT) field for the current sample
                    genotype = record.samples[sample]['GT']
                    if genotype is not None:
                        alleles = [record.ref] + list(record.alts)  # List of alleles (ref + alts)
                        allele = [alleles[g] if g is not None else '.' for g in genotype]
                        if allele[0] != ancestral:
                            pos_entry.append("Derived")
                            pos_derived_allele.append(allele[0])
                        elif allele[0] == ancestral:
                            pos_entry.append("Ancestral")
                    else:
                        pos_entry.append(".")
            else:
                pos_entry.append("-") #for ancestral
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
    output.to_csv("polaryzer_output.csv", sep="\t", index=False)
    vcf_in.close()
    print(f"Finished processing {args.multi_sample_vcf}")

else:
    print("No input vcf file was provided! Polaryzer aborted. Please provide path to vcf folder with several singel vcf files (-input_single_vcf) or path to a multi sample vcf file (-multi_sample_vcf).")
print("Finished polarizing")