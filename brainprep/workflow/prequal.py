# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2021 - 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
Interface for prequal.
"""

# Imports
import os
import glob
import shutil
import tempfile
from brainprep.color_utils import (
    print_result, print_subtitle, print_title, print_command)
import pandas as pd
import subprocess


def brainprep_prequal(dwi, bvec, bval, pe, readout_time, output_dir, t1=None):
    """ Define the dMRI pre-processing workflow.

    Parameters
    ----------
    dwi: str
        path to the diffusion weighted image.
    bvec:
        path to the bvec file.
    bval:
        path to the bval file.
    pe: str
        the de phase encoding direction (i, i-, j, j-, k, k-).
    readout_time: str
        readout time of the dwi image.
    output_dir: str
        path to the output directory.
    t1: str
        path to the t1 image in case of synb0 use.

    Notes
    -----
    In order to use the synb0 feature you must bind your freesurfer license as
    such: -B /path/to/freesurfer/license.txt:/APPS/freesurfer/license.txt
    """
    print_title("PreQual dtiQA pipeline")
    if pe in ["i", "j", "k"]:
        pe_axis = pe
        pe_signe = "+"
    elif pe in ["i-", "j-", "k-"]:
        pe_axis = pe[0]
        pe_signe = pe[1]
    else:
        raise Exception("Valid input for pe are (i, i-, j, j-, k, k-)")

    print_subtitle("Making dtiQA_config.csv")
    dtiQA_config = [os.path.basename(dwi).split('.')[0],
                    pe_signe,
                    readout_time]
    df_dtiQA_config = pd.DataFrame(dtiQA_config)
    print_result("dtiQA_config file content :\n")
    print_result(dtiQA_config)

    print_subtitle("Copy before launch")
    with tempfile.TemporaryDirectory() as tmpdir:
        df_dtiQA_config.T.to_csv(os.path.join(tmpdir, "dtiQA_config.csv"),
                                 sep=",", header=False, index=False)
        shutil.copy(dwi, tmpdir)
        shutil.copy(bvec, tmpdir)
        shutil.copy(bval, tmpdir)
        if t1 is not None:
            shutil.copy(t1, os.path.join(tmpdir, "t1.nii.gz"))

        print_subtitle("Launch prequal...")
        cmd = ["xvfb-run",  "-a", "--server-num=1",
               "--server-args='-screen 0 1600x1280x24 -ac'",
               "bash", "/CODE/run_dtiQA.sh", tmpdir, output_dir, pe_axis]
        print_command(" ".join(cmd))
        with subprocess.Popen(cmd,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT) as process:
            for line in process.stdout:
                print(line.decode("utf8"))


def brainprep_prequal_qc(data_regex, outdir, sub_idx=-4):
    """ Define the dMRI pre-processing quality control workflow.

    Parameters
    ----------
    datadir: str
        regex to the dmriprep 'stats.csv' files.
    outdir: str
        path to the destination folder.
    sub_idx: int, default -4
        the position of the subject identifier in the input path.
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn import preprocessing

    print_title("Loading PreQual dtiQA stats files...")
    stat_files = glob.glob(data_regex)
    stats = []
    for path in stat_files:
        df = pd.read_csv(path, index_col=0, header=None).T
        sid = path.split(os.sep)[sub_idx]
        df["participant_id"] = sid
        stats.append(df)
    df = pd.concat(stats)
    df.to_csv(os.path.join(outdir, "transformation.tsv"), sep="\t",
              index=False)

    print_title("Computing box plot by category...")
    for score_name in ("fa", "md", "rd", "ad"):
        _cols = [name for name in df.columns
                 if name.endswith(f"_{score_name}")]
        _df = pd.melt(df, id_vars=["participant_id"], value_vars=_cols,
                      var_name="metric", value_name="value")
        le = preprocessing.LabelEncoder()
        _df["label"] = le.fit_transform(_df.metric.values)
        sns.set(style="whitegrid")
        ax = sns.boxplot(x="label", y="value", data=_df)
        plt.setp(ax.get_xticklabels(), rotation=90)
        ax.tick_params(labelsize=8)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{score_name}.png"), dpi=400)
    _cols = ["eddy_avg_rotations_x", "eddy_avg_rotations_y",
             "eddy_avg_rotations_z", "eddy_avg_translations_x",
             "eddy_avg_translations_y", "eddy_avg_translations_z",
             "eddy_avg_abs_displacement", "eddy_avg_rel_displacement"]
    _df = pd.melt(df, id_vars=["participant_id"], value_vars=_cols,
                  var_name="metric", value_name="value")
    sns.set(style="whitegrid")
    ax = sns.boxplot(x="metric", y="value", data=_df)
    plt.setp(ax.get_xticklabels(), rotation=70)
    ax.tick_params(labelsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "transformation.png"), dpi=400)
