import pytest
import shutil
import json
import os
import subprocess
import numpy as np

DATADIR = "data"
SRCDIR="/corexfs/schroedj/src/IPODerator"
TARGETDIR = "targets"
RESULTDIR = "results"

ARGS = f"--info {os.path.join(DATADIR, 'very_minimal.csv')} "\
    f"--draws_direc {os.path.join(RESULTDIR, 'draws')} "\
    f"--debug_outdir {RESULTDIR} "\
    f"--debug "\
    f"--data_file {os.path.join(RESULTDIR, 'test_data.RData')}"

FITRESULTDIR = "fit_results"
FITARGS = f"--stan_file {os.path.join(SRCDIR, 'stan/sig_noise_alloc.stan')} "\
    f"--draws_direc {os.path.join(FITRESULTDIR, 'draws')} " \
    f"--seed 42 " \
    f"--fit_file {os.path.join(FITRESULTDIR, 'fit.RData')} " \
    f"--debug_outdir {FITRESULTDIR} "\
    f"--debug "\
    f"--data_file {os.path.join(RESULTDIR, 'test_data.RData')}"

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
 

def get_par_idx_array():
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
    print("stdout from get_par_idx_array")
    print("=============================")
    print(result.stdout.decode())
    print()
    print("stderr from get_par_idx_array")
    print("=============================")
    print(result.stderr.decode())

#get_par_idx_array()
fit_model()

def test_get_par_num():
    target = {
        "alpha": 4,
        "beta": 1,
        "gamma": 2,
    }
    # read entire result design matrix to string
    with open(os.path.join(RESULTDIR, "nums.json"), "r") as result_file:
        result = json.load(result_file)
    assert all([v == result[k] for k,v in target.items()])

def test_get_par_idx_array():
    # read entire ground truth csv to string
    target = np.loadtxt(os.path.join(TARGETDIR, "param_idx_arr.csv"), delimiter=",")
    # read entire result design matrix to string
    result = np.loadtxt(os.path.join(RESULTDIR, "param_idx_arr.csv"), delimiter=",")
    assert np.all(target == result)

def test_get_gamma_idx():
    # read entire ground truth csv to string
    target = np.loadtxt(os.path.join(TARGETDIR, "gamma_idx.csv"), delimiter=",")
    # read entire result design matrix to string
    result = np.loadtxt(os.path.join(RESULTDIR, "gamma_idx.csv"), delimiter=",")
    assert np.all(target == result)

def test_get_expt_key():
    # read entire ground truth json to list
    with open(os.path.join(TARGETDIR, "expt_key.json"), "r") as target_file:
        target = json.load(target_file)
    # read entire result design matrix to string
    with open(os.path.join(RESULTDIR, "draws/expt_key.json"), "r") as result_file:
        result = json.load(result_file)
    assert np.all(target == result)

