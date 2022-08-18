##############################################################################
# columns.py
#
# File containing dataclasses that define the header columns used by the 
# .ssm and .params.json file formats.
##############################################################################

from dataclasses import dataclass, field

@dataclass(frozen=True)
class SSM_Columns:
    """Dataclass used to define .ssm column headers"""
    NAME: str = "name"
    ID: str = "id"
    VAR_READS: str = "var_reads"
    TOTAL_READS: str = "total_reads"
    VAR_READ_PROB: str = "var_read_prob"

    # order in which columns should appear in .ssm file
    COL_ORDER: tuple = ("id", "name", "var_reads", "total_reads","var_read_prob")


@dataclass(frozen=True) 
class PARAMS_Columns:
    """Dataclass used to define .params.json column headers"""
    SAMPLES: str = "samples"
    CLUSTERS: str = "clusters"
    GARBAGE: str = "garbage"