"""Shared input and logging utilities."""

import logging
import os

import numpy as np
import pandas as pd

from ._version import __version__
from .art import logo


def get_openaschem_logger():
    logging.basicConfig(
        level=logging.INFO, format="OpenASChem: %(asctime)s - %(levelname)s - %(message)s"
    )
    logger = logging.getLogger()

    return logger


def welcome(logger):
    logger.info(logo)
    logger.info("Welcome to OpenASChem %s!", __version__)


def check_raw_data_files(file_dir, standard_prefix, logger):
    def get_run_name(x):
        if x.endswith(".mzXML"):
            return x[:-6]
        else:
            return x[:-5]

    rawdata_files = [
        i for i in sorted(os.listdir(file_dir)) if i.endswith(".mzML") or i.endswith(".mzXML")
    ]
    if not rawdata_files:
        raise Exception("Cannot find any valid MS raw data files.")

    standard_files = [i for i in rawdata_files if i.startswith(standard_prefix)]
    chemical_pools = set(list(map(lambda x: x.split("_")[1], standard_files)))
    logger.info(
        "%d chemical pools found: %s" % (len(chemical_pools), ";".join(list(chemical_pools)))
    )

    raw_data_info = {}
    for pool in chemical_pools:
        standard_files_pool = [i for i in standard_files if pool in i]
        if len(standard_files_pool) != 3:
            logger.error("Invalid chemical pool: %s. Not enough standard runs." % pool)
            continue
        experimental_files_pool = [
            i for i in rawdata_files if ((pool in i) and (not i.startswith(standard_prefix)))
        ]
        experimental_runs_pool = [get_run_name(file) for file in experimental_files_pool]
        raw_data_info[pool] = [standard_files_pool, experimental_files_pool, experimental_runs_pool]

    return raw_data_info


def load_chemical_pool(chemical_pool_file, acquisition_mode, chemical_pool_columns):
    if chemical_pool_file.endswith(".tsv"):
        chemical_pool = pd.read_csv(chemical_pool_file, sep="\t")
    elif chemical_pool_file.endswith(".csv"):
        chemical_pool = pd.read_csv(chemical_pool_file)
    else:
        chemical_pool = pd.read_excel(chemical_pool_file)

    if acquisition_mode == "positive":
        chemical_pool[chemical_pool_columns["COMPOUND_MZ"]] += 1.007285
    else:
        chemical_pool[chemical_pool_columns["COMPOUND_MZ"]] -= 1.007285

    chemical_pool = chemical_pool.dropna(subset=[chemical_pool_columns["COMPOUND_ID"]])

    return chemical_pool


def split_list(input_list, num_sublists):
    if num_sublists <= 0:
        raise ValueError("Number of sub-lists must be greater than zero.")

    # Calculate the size of each chunk
    avg_size = len(input_list) // num_sublists
    remainder = len(input_list) % num_sublists

    # Create sublists
    sublists = [
        input_list[i * avg_size + min(i, remainder) : (i + 1) * avg_size + min(i + 1, remainder)]
        for i in range(num_sublists)
    ]

    return sublists


def jitter_array(a_array, scale=1):
    """Preserve the existing small perturbation used in statistical comparisons."""
    return a_array + np.random.uniform(-scale, scale, a_array.shape)
