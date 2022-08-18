##############################################################################
# ssm_base_converter.py
#
# File containing the class definition for the SSM_Base_Converter.
# The SSM_Base_Converter class can be subclassed and used for
# generating an .ssm and .params.json file from another file type (e.g. VCF).
#
# The trick this class does is search all of its subclasses functions
# for a prefix "p_". The functions that are prefixed with "p_" 
# (e.g. "p_function1()") must be intended to process the pandas 
# dataframe (self.in_df) which contains the data from the starting format (e.g. vcf). 
# Once all processing functions ("p_*") are called, then the .ssm and 
# .params.json files can be written.
#
# For an example of how to use this class, please see the vcf_to_ssm.py
# source file.
##############################################################################

import pandas as pd
from columns import SSM_Columns, PARAMS_Columns

class SSM_Base_Converter:
    """
    Base class that can be used to convert to an equivalent .ssm and .params.json file.
    For an example of its use, see 'vcf_to_ssm.py'.
    """

    def __init__(self, in_file):
        """Constructor"""
        # set up everything necessary to read/process/write
        self._init_constants()
        self._init_variables()

        # start processing if we have all of the information we need (I/O file names)
        if in_file:

            # read and set the in-file
            self.read_in_file(in_file)

            # run all processing functions
            self.process()

    def _init_constants(self):
        """Initializes constants that are used throughout"""
        # P_SIGNATURE is a prefix placed in front of a function which is meant to 
        # process something from the input format (e.g. vcf) to the .ssm format
        self.P_SIGNATURE = "p_.*"

    def _init_variables(self):
        """Initializes dataframes and processing functions list"""
        # initializations
        self.in_df = None
        self.out_df = None

        self.processing_functions = [
            # all functions used to translate input file to SSM file
        ]

    def read_in_file(self, in_file):
        """
        This function must be overriden in a subclass to read some file format of interest. It should 
        set the self.in_df which must be a Pandas DataFrame type.
        """
        # must set self.in_df such that isinstance(self.in_df, pd.DataFrame)
        # evaluates to True after this function is called
        raise NotImplementedError

    def format_out_df(self):
        """
        Empty base function that should be overridden in child class that properly configures 
        self.out_df with the necessary columns to be written out as a .ssm file
        """
        # self.out_df must contain all of the columns that are in the .ssm file format.
        raise NotImplementedError

    def write_ssm_file(self, ssm_fn):
        """Writes self.out_df to file"""

        assert isinstance(self.out_df, pd.DataFrame),  "self.out_df must be a Pandas DataFrame"
        assert set(SSM_Columns.COL_ORDER).issubset(self.out_df.columns), \
               "self.out_df does not contain all of the columns necessary for an .ssm file. \
                The columns for self.out_df are: " + ",".join(self.out_df.columns.values)
        
        self.out_df[list(SSM_Columns.COL_ORDER)].to_csv(ssm_fn, sep="\t", index=False)

    def write_out_params(self, params_file):
        """Writes out a .params.json file using the overrides of the 
        self.sample_names() and self.garbage_mutations() functions"""
        import json

        def convert(obj):
            """
            Overrides default function when something isn't json serializable
            """
            import numpy as np

            if isinstance(obj, np.ndarray): return obj.tolist()
            elif isinstance(obj, list): pass
            else: raise TypeError

        with open(params_file, 'w') as outfile:

            json_dict = json.dumps(
                            {
                                PARAMS_Columns.SAMPLES: self.sample_names(),
                                PARAMS_Columns.GARBAGE: self.garbage_mutations()
                            },
                            default=convert
                        )

            outfile.write(json_dict)

    def garbage_mutations(self):
        """
        This may be replaced or implemented at the subclass level to detect garbage mutations to ignore.
        """
        return []

    def sample_names(self):
        """
        This may be replaced or implemented at the subclass level to return the list of sample names.
        If it's not, we'll assume there is a column 'samples' in self.in_df which lists the sample names.
        """
        assert set([PARAMS_Columns.SAMPLE_NAMES]).issubset(self.in_df.columns), \
               "self.in_df must contain a column 'samples' which lists the name of the samples the variants were called from."

        return self.in_df[PARAMS_Columns.SAMPLE_NAMES].unique()

    def process(self):
        """
        Collects all processing functions (functions that match r'p\_.*'), calls all processing
        functions, then transforms the processed dataframe into a dataframe that can written as a .ssm
        """
        from tqdm import tqdm

        # collect all processing functions
        self._aggregate_processing_functions()

        # set up progress bar with all processing functions
        pbar = tqdm(self.processing_functions)

        # call all processing functions
        for function in pbar:

            pbar.set_description("Running %s" % function.__name__.lstrip(self.P_SIGNATURE))

            function()

            pbar.set_description("Completed.")

        pbar.set_description("All processing complete.")

        # translate processed_df into out_df for writing to file
        self.format_out_df()

    def _aggregate_processing_functions(self):
        """
        Searches through object attributes and methods for functions that match r'p\_.*'
        """
        import re

        regex = re.compile(self.P_SIGNATURE)

        # finds all attributes which match our prefix signature (should be processing functions)
        processing_functions = list(filter(regex.match, dir(self)))

        self.processing_functions = [getattr(self, f_name) for f_name in processing_functions]