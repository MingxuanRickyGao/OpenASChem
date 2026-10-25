"""Legacy file replication helper."""

import os
import shutil

import pandas as pd

from .utils import get_openaschem_logger, welcome


def open_filling(file_dir, out, n_replicates):
    logger = get_openaschem_logger()
    welcome(logger)

    logger.info("Filling missing data replicates.")

    def get_run_name(x):
        if x.endswith(".mzXML"):
            return x[:-6]
        else:
            return x[:-5]

    def get_file_type(x):
        if x.endswith(".mzXML"):
            return ".mzXML"
        else:
            return ".mzML"

    rawdata_files = [
        i for i in sorted(os.listdir(file_dir)) if i.endswith(".mzML") or i.endswith(".mzXML")
    ]
    if not rawdata_files:
        raise Exception("Cannot find any valid MS raw data files.")

    rawdata_runs = [get_run_name(file) for file in rawdata_files]
    rawdata_samples = list(set(["_".join(run.strip().split("_")[:-1]) for run in rawdata_runs]))

    filling_logs = {"original_run": [], "filling_run": []}

    for sample in rawdata_samples:
        existing_files = [file for file in rawdata_files if file.startswith(sample)]
        existing_runs = [get_run_name(file) for file in existing_files]

        for replicate_id in range(1, n_replicates + 1):
            target_run = "%s_%d" % (sample, replicate_id)
            if target_run not in existing_runs:
                source_file = os.path.join(file_dir, existing_files[0])
                dst_file = os.path.join(
                    file_dir, "%s_%d%s" % (sample, replicate_id, get_file_type(existing_files[0]))
                )
                shutil.copy(source_file, dst_file)
                filling_logs["original_run"].append(existing_runs[0])
                filling_logs["filling_run"].append(target_run)
                logger.info(
                    "%s is missing. Use %s to fill the missing file."
                    % (existing_runs[0], target_run)
                )

    filling_logs = pd.DataFrame(filling_logs)

    if filling_logs.shape[0] != 0:
        logger.info(
            "There are %d missing files in total. Be careful when you perform downstream analysis on these runs. Check this file to see which runs are filling runs: %s"
            % (filling_logs.shape[0], os.path.join(out, "file_filling_logs.tsv"))
        )
        filling_logs.to_csv(os.path.join(out, "file_filling_logs.tsv"), index=False, sep="\t")
    else:
        logger.info("File filling done! No missing file was found!")

    logger.info("Done!")
