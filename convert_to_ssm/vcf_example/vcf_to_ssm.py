##############################################################################
# vcf_to_ssm.py
#
# Subclass of SSM_Base_Converter meant to be an example of how to
# translate some file format (e.g. vcf) to an equivalent .ssm and .params.json
# that can be used by Pairtree.
#
# Please note that this is not meant to be a complete and final solution. 
#
##############################################################################

import os, sys
import pandas as pd
import argparse 
import re
from dataclasses import dataclass

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from ssm_base_converter import SSM_Base_Converter
from columns import SSM_Columns

@dataclass
class INFO_Fields:
    """Dataclass used to define VCF INFO field tags"""
    VCF_ID: str = "ID" # id field
    FORMAT: str = "FORMAT" # format field

@dataclass
class FORMAT_Fields:
    """Dataclass used to define VCF FORMAT field tags"""
    GT: str = "GT" # genotype
    RD: str = "RD" # reference depth 
    AD: str = "AD" # alternate allele depth 
    CN: str = "CN" # copy number for each allele
    REF_GT: str = "0" # reference genotype number (we assume it's 0)
    ALT_GT: str = "1" # variant genotype number (we assume it's 1)

# simple lambda function for producing array strings 
no_brackets_array_string = lambda arr: ", ".join(map(str, arr))

#########################################################################################
# How to run using example.vcf file:
#
#   python3 vcf_to_ssm.py example.vcf example --write-params
#
#########################################################################################

class VCF_To_SSM(SSM_Base_Converter):
    """Basic class to demonstrate how we could translate a VCF file to an equivalent .ssm and .params.json.
    We make some large assumptions (e.g. there is only 1 ALT for each variant, the FORMAT VCF columns have the necessary
    fields specified, etc.)."""

    def __init__(self, in_file):
        """Constructor"""
        super().__init__(in_file)

    def read_in_file(self, in_file):
        """Override of SSM_Base_Converter.read_in_file() which reads a VCF file into a Pandas DataFrame"""
        # initialize the variables used to construct our Pandas DataFrame in case there's nothing in the file
        header, variants = [], [] 

        # opens the vcf file, and pulls out the columns and variant data
        with open(in_file, "r") as f:
            lines = f.readlines()

            for i, line in enumerate(lines):
                if line.strip().startswith("#CHROM"):
                    break

            data = lines[i:]
            header = data[0].strip().split("\t")
            variants = [d.strip().split("\t") for d in data[1:]]

        # THIS IS REQUIRED! Sets self.in_df to a Pandas DataFrame
        self.in_df = pd.DataFrame(variants, columns=header)

        # keep a separate dataframe to access only the sample data
        sample_data = self.in_df.iloc[:, self.in_df.columns.get_loc(INFO_Fields.FORMAT):].apply(lambda x: x.str.split(":"))

        # separate out formats for each variant from variant data
        self.formats = sample_data[INFO_Fields.FORMAT]

        # keep sample data and sort the samples by name
        self.sample_data = sample_data.drop(INFO_Fields.FORMAT, axis=1)
        self.sample_data = self.sample_data[sorted(self.sample_data.columns)]

    def format_out_df(self):
        """
        Override to move the newly added columns to self.in_df 
        to self.out_df which we'll use to write out to an .ssm
        """

        # initialize out_df
        self.out_df = pd.DataFrame()     

        # set name
        self.out_df[SSM_Columns.NAME] = self.in_df[SSM_Columns.NAME]

        # set var_reads
        self.out_df[SSM_Columns.VAR_READS] = self.in_df[SSM_Columns.VAR_READS]

        # set total_reads
        self.out_df[SSM_Columns.TOTAL_READS] = self.in_df[SSM_Columns.TOTAL_READS]

        # set var_read_prob
        self.out_df[SSM_Columns.VAR_READ_PROB] = self.in_df[SSM_Columns.VAR_READ_PROB]

        # set ID
        # Provides each variant with a unique id (r's\d+').
        self.out_df[SSM_Columns.ID] = ["s" + str(number) for number in list(range(0, len(self.out_df)))]

    def p_names(self):
        """Processes name by extracting the VCF ID directly"""
        self.in_df[SSM_Columns.NAME] = self.in_df[INFO_Fields.VCF_ID].apply(str)

    def p_reads(self):
        """Processes read count data by indexing the variant read count and 
        reference read count using the format defined for that row"""
        for idx, row in self.sample_data.iterrows():

            # compile reference reads
            if FORMAT_Fields.RD in self.formats.iloc[idx]:
                reference_reads = row.apply(lambda s: int(s[self.formats[idx].index(FORMAT_Fields.RD)]))

            # compile variant reads
            if FORMAT_Fields.AD in self.formats.iloc[idx]:
                variant_reads = row.apply(lambda s: int(s[self.formats[idx].index(FORMAT_Fields.AD)]))

            # place processed var_reads and total_reads back into our self.in_df so we can access it later
            self.in_df.loc[idx, SSM_Columns.VAR_READS] = no_brackets_array_string(variant_reads.tolist())
            self.in_df.loc[idx, SSM_Columns.TOTAL_READS] = no_brackets_array_string((reference_reads + variant_reads).tolist())

    def p_var_read_prob(self):
        """Computes variant read probability by indexing the variant and reference copy 
        numbers using the format defined for that row"""
        for idx, row in self.sample_data.iterrows():
            if FORMAT_Fields.CN in self.formats.iloc[idx]:

                # assume genotypes are separated by either '|' or '/'
                genotypes = row.apply(lambda s: re.split('[\|\/]', s[self.formats[idx].index(FORMAT_Fields.GT)]))

                # assume copy numbers are separated by ','
                copy_numbers = row.apply(lambda s: s[self.formats[idx].index(FORMAT_Fields.CN)].split(","))

                # extract variant read probabilities using genotype and copy number
                var_read_probs = []
                for gts, cns in zip(genotypes, copy_numbers):
                    variant_cn = int(cns[gts.index(FORMAT_Fields.ALT_GT)])
                    ref_cn = int(cns[gts.index(FORMAT_Fields.REF_GT)])
                    var_read_probs.append(variant_cn / (variant_cn + ref_cn))

                # place our var_read_prob back into our self.in_df so we can access it later
                self.in_df.loc[idx, SSM_Columns.VAR_READ_PROB] = no_brackets_array_string(var_read_probs)

    def sample_names(self):
        """Override of base class. Used to return the sorted sample names from the VCF file."""
        return self.sample_data.columns.values


def main():
    """
    Parses arguments and uses VCF_TO_SSM class to generate a .params.json and .ssm file from a .vcf file
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("in_file", 
                        help="vcf file that will be converted to an .ssm and .params.json file")
    parser.add_argument("out_file", 
                        help="Name which will be used to name the .ssm and .params.json file that \
                             are outputted (e.g. out_file.ssm, out_file.params.json")
    parser.add_argument("--write-params", default=True, action="store_true",
                        help="Flag used to tell program to write out a .params.json file")

    args = parser.parse_args()


    print("Processing vcf file ...")

    # initializing class calls all processing functions and prepares the out_df
    vcf_to_ssm = VCF_To_SSM(args.in_file)

    print("Writing ssm file ...")

    # write SSM file
    vcf_to_ssm.write_ssm_file(args.out_file + ".ssm")

    # write the params file
    if args.write_params:
        print("Writing params.json file ...")
        vcf_to_ssm.write_out_params(args.out_file + ".params.json")

    print("Completed.")

if __name__ == "__main__":
    main()