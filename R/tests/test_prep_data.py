import pytest
import shutil
import json
import os
import subprocess
import numpy as np

DATADIR = "data"
SRCDIR="/corexfs/schroedj/src/IPODerator"
TARGETDIR = "targets"
RANK_TARGETDIR = "rank_targets"
RESULTDIR = "results"
RANK_RESULTDIR = "rank_results"

#ARGS = "--cfg_file config.toml"
#FITARGS = "--cfg_file fit_config.toml"
FITRESULTDIR = "fit_results"
TESTPARSECSVDIR = "parse_csv_testing"
TESTPARSEBLOBDIR = "parse_blob_testing"

ARGS = f"--info {os.path.join(DATADIR, 'very_minimal.csv')} "\
    f"--data_direc {os.path.join(RESULTDIR, 'data')} "\
    f"--debug_outdir {RESULTDIR} "\
    f"--debug "

RANK_ARGS = f"--info {os.path.join(DATADIR, 'rank_info.csv')} "\
    f"--data_direc {os.path.join(RANK_RESULTDIR, 'data')} "\
    f"--debug_outdir {RANK_RESULTDIR} "\
    f"--debug "

FITARGS = f"--stan_file {os.path.join(SRCDIR, 'stan/test_control_sub.stan')} "\
    f"--draws_direc {os.path.join(FITRESULTDIR, 'draws')} " \
    f"--seed 42 " \
    f"--fit_file {os.path.join(FITRESULTDIR, 'fit.RData')} " \
    f"--debug_outdir {FITRESULTDIR} "\
    f"--debug "\
    f"--data_file {os.path.join(RESULTDIR, 'data.RData')}"

PARSEARGS = f"--cfg-file {os.path.join(TESTPARSECSVDIR, 'draws/config.toml')}"

SUMMARIZEARGS = f"--cfg-file {os.path.join(TESTPARSEBLOBDIR, 'draws/config.toml')}"

def fit_model():
    if os.path.isdir(FITRESULTDIR):
        rm_files = os.listdir(FITRESULTDIR)
        for f in rm_files:
            thingy = os.path.join(FITRESULTDIR, f)
            if os.path.isdir(thingy):
                shutil.rmtree(thingy)
            else:
                os.remove(thingy)
    cmd = f"Rscript {SRCDIR}/R/cmdstan_fit_enrichment_model.R {FITARGS}"
    result = subprocess.run(cmd, shell=True, capture_output=True)
    print()
    print("stderr from fit_model")
    print("=============================")
    print(result.stderr.decode())
    print()
    print("stdout from fit_model")
    print("=============================")
    print(result.stdout.decode())

def parse_csv_to_blob():
    # remove prior output if it exists
    out_fname = os.path.join(TESTPARSECSVDIR, "params.blob")
    if os.path.isfile(out_fname):
        os.remove(out_fname)

    cmd = f"{SRCDIR}/bin/parse_csv_to_blob {PARSEARGS}"
    result = subprocess.run(cmd, shell=True, capture_output=True)
    print()
    print("stdout from parse_csv_to_blob")
    print("=============================")
    print(result.stdout.decode())
    print()
    print("stderr from parse_csv_to_blob")
    print("=============================")
    print(result.stderr.decode())

def parse_enricherator_blob():
    # remove prior output if it exists
    out_fname = os.path.join(TESTPARSEBLOBDIR, "params.blob")
    if os.path.isfile(out_fname):
        os.remove(out_fname)

    cmd = f"{SRCDIR}/bin/parse_enricherator_blob {SUMMARIZEARGS}"
    result = subprocess.run(cmd, shell=True, capture_output=True)
    print()
    print("stdout from parse_enricherator_blob")
    print("=============================")
    print(result.stdout.decode())
    print()
    print("stderr from parse_enricherator_blob")
    print("=============================")
    print(result.stderr.decode())

def prep_data():
    if os.path.isdir(RESULTDIR):
        rm_files = os.listdir(RESULTDIR)
        for f in rm_files:
            thingy = os.path.join(RESULTDIR, f)
            if os.path.isdir(thingy):
                shutil.rmtree(thingy)
            else:
                os.remove(thingy)
    cmd = f"Rscript {SRCDIR}/R/prep_data.R {ARGS}"
    result = subprocess.run(cmd, shell=True, capture_output=True)
    print()
    print("stdout from prep_data")
    print("=============================")
    print(result.stdout.decode())
    print()
    print("stderr from prep_data")
    print("=============================")
    print(result.stderr.decode())

def prep_rank_data():
    if os.path.isdir(RANK_RESULTDIR):
        rm_files = os.listdir(RANK_RESULTDIR)
        for f in rm_files:
            thingy = os.path.join(RANK_RESULTDIR, f)
            if os.path.isdir(thingy):
                shutil.rmtree(thingy)
            else:
                os.remove(thingy)
    cmd = f"Rscript {SRCDIR}/R/prep_data.R {RANK_ARGS}"
    result = subprocess.run(cmd, shell=True, capture_output=True)
    print()
    print("stdout from prep_rank_data")
    print("=============================")
    print(result.stdout.decode())
    print()
    print("stderr from prep_rank_data")
    print("=============================")
    print(result.stderr.decode())

prep_data()
prep_rank_data()
fit_model()
parse_csv_to_blob()
parse_enricherator_blob()

def test_get_par_num():
    target = {
        "alpha": 4,
        "beta": 1,
        "gamma": 2,
    }
    with open(os.path.join(RESULTDIR, "nums.json"), "r") as result_file:
        result = json.load(result_file)
    assert all([v == result[k] for k,v in target.items()])

def test_get_par_idx_array():
    # read entire ground truth csv to string
    target = np.loadtxt(os.path.join(TARGETDIR, "param_idx_arr.csv"), delimiter=",")
    result = np.loadtxt(os.path.join(RESULTDIR, "param_idx_arr.csv"), delimiter=",")
    assert np.all(target == result)

def test_get_beta_sub_idx():
    # read entire ground truth csv to string
    target = np.loadtxt(os.path.join(RANK_TARGETDIR, "sub_inds1.csv"), delimiter=",")
    result = np.loadtxt(os.path.join(RANK_RESULTDIR, "sub_inds1.csv"), delimiter=",")
    assert np.all(target == result)

def test_get_expt_key():
    # read entire ground truth json to list
    with open(os.path.join(TARGETDIR, "expt_key.json"), "r") as target_file:
        target = json.load(target_file)
    with open(os.path.join(RESULTDIR, "data/expt_key.json"), "r") as result_file:
        result = json.load(result_file)
    assert np.all(target == result)

