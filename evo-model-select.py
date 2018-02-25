#!/usr/bin/env python3
version = "1.0.1" #parallel features
description = \
"""

    %prog --fasta_file FASTA --model_list MODEL_LIST

    Code for a simple model tester for ML analysis using Garli for
    likelihood calculations.

    version {}
    by Nicolás D. Franco-Sierra

    Expected inputs:
    (1) A nucleotide alignment in FASTA format
    (2) A model list (model name) to be evaluated in a plain text file.
    One model per line (eg.
                                GTR
                                JC+I+G
                                K80+G)
    Invariant sites and gamma distributions are supported.

    Output: a tab separated file table containing likelihood scores for
    each requested model

    What it does? it creates the proper garli.conf file and performs a
    search with 1 rep for each model. Then it gets the likelihood
    score for the search.
""".format(version)

#imports
import sys, os
from optparse import OptionParser
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import subprocess
import time
from datetime import timedelta
from multiprocessing.dummy import Pool as ThreadPool
import math

start_time = time.monotonic() # I want to record execution time.


exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
sys.path.insert(0,os.path.abspath(os.path.join(exe_path,"..", "..")))
exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]

parser = OptionParser(version="%prog {}".format(version),
            usage = description)
parser.add_option("-f", "--fasta_file", dest="fasta_file",
                    metavar="FILE",
                    type="string",
                    help="Name and path of input alignment in FASTA format. ")
parser.add_option("-m", "--model_list", dest="model_list",
                    metavar="FILE",
                    type="string",
                    #default="all",
                    help="Name and path of file contanining models to be \
                    tested.")
parser.add_option("-o", "--result_file", dest="result_file",
                    metavar="FILE",
                    type="string",
                    default="final.phylo_results.tab",
                    help="Name and path of where the files should end up. ")
parser.add_option("-t", "--temp_directory", dest="temp_directory",
                    metavar="PATH",
                    type="string",
                    default="tmp",
                    help="Name and path of temporary directory where \
                        calculations should take place. ")
parser.add_option("-j", "--jobs", dest="jobs",
                    default=1,
                    metavar="jobs",
                    type="int",
                    help="Specifies the number of jobs (operations) to run \
                    in parallel.")

(options, remaining_args) = parser.parse_args()

jobs = options.jobs

if not options.fasta_file:
    parser.error("\n\n\tMissing parameter --fasta_file FILE\n\n")
if not options.model_list:
    parser.error("\n\n\tMissing parameter --model_list FILE\n\n")

def run_cmd(cmd_str, work_path=None):
    """
    This helper function excecute a bash command as a separate process
    and returns and expection if the command fails.

    It returns the standart output and the standard error of the command as
    str objects.

    """
    process = subprocess.Popen(cmd_str, stdout = subprocess.PIPE, stderr =
    subprocess.PIPE, cwd=work_path)
    stdout_str, stderr_str = process.communicate()
    if process.returncode != 0:
        raise Exception("Failed to run '%s'\n%s%sNon-zero exit status %s" %
                            (cmd_str, stdout_str, stderr_str,
                            process.returncode))
    return stdout_str.decode(), stderr_str.decode()

def model_parser(model_string):
    """
    Parses a model string to the equivalent command block for Garli
    e.g. converts "GTR+I+G" to: "ratematrix = 6rate
                                 statefrequencies = estimate
                                 ratehetmodel = gamma
                                 numratecats = 4
                                 invariantsites = estimate"

    """
    #define equivalencies in dictionaries
    model_dict = {"GTR":"ratematrix = 6rate\nstatefrequencies = estimate\n",
                  "SYM":"ratematrix = 6rate\nstatefrequencies = equal\n",
                  "HKY":"ratematrix = 2rate\nstatefrequencies = estimate\n",
                  "K80":"ratematrix = 2rate\nstatefrequencies = equal\n",
                  "F81":"ratematrix = 1rate\nstatefrequencies = estimate\n",
                  "JC":"ratematrix = 1rate\nstatefrequencies = equal\n"}
    gamma_dict = {True:"ratehetmodel = gamma\nnumratecats = 4\n",
                    False:"ratehetmodel = none\nnumratecats = 1"}
    inv_dict = {True:"invariantsites = estimate\n",
                False:"invariantsites = none"}
    # parse model name
    model_info = model_string.split("+")
    model = model_info[0]
    inv_status = True if "I" in model_info else False
    gamma_status = True if "G" in model_info else False
    # compose Garli block
    composed_string = model_dict[model] + gamma_dict[gamma_status] + \
                      inv_dict[inv_status]
    return composed_string

def model_reader(model_file):
    """
    Reads a model file list and returns models with commands.
    It opens the file and retrieve a dictionary with the model
    names paired to their respective command block for Garli

    """
    models_to_eval = {}
    with open(model_file) as models:
        model_count = len(models.readlines())
    print("A total of {} DNA evolution models to be \
processed.".format(model_count))
    with open(model_file) as models:
        for model in models:
            model = model.strip()
            models_to_eval[model] = model_parser(model)

    return models_to_eval

def name_reformat(name):
    """
    This is a helper function to reformat FASTA headers from sequences
    retrieved from GenBank
    e.g.
    It converts headers of the form: gi|984656999|gb|KU593393.1| Turdus ignobilis voucher AMNH:DOT4802 glyceraldehyde-3-phosphate dehydrogenase (G3PDH) gene, partial cds
    To safe sequence names such as: 'Turdus_ignobilis_voucher_AMNH:DOT4802'

    """
    safe_name = "_".join(name.split("|")[-1].split()[0:4])
    return safe_name

def fasta2nexus(fasta_filename, location):
    """
    This helper function takes a filename corresponding to a FASTA file (must
    have the ".fasta" extension) and opens the file to transform it to NEXUS
    format using BioPython utilities.
    The function returns a a tuple containing the followinf info:
    (filename_of_the_new_NEXUS_file, number_of_taxa_in_file, number_of_chars
    _in_the_alignment)

    """
    nex_filename = fasta_filename.replace("fasta", "nexus")
    with open(fasta_filename) as fasta, open(os.path.join(location,
        nex_filename),"a") as nexus:
        sequences = SeqIO.parse(fasta, "fasta", alphabet=IUPAC.ambiguous_dna)
        new_seqs = []
        for sequence in sequences:
            new_name = name_reformat(sequence.description)
            sequence.id = new_name
            new_seqs.append(sequence)
        SeqIO.write(new_seqs, nexus, "nexus")
        ntax = len(new_seqs)
        nchars = len(new_seqs[0].seq)
    return nex_filename, ntax, nchars

def garliconf_gen(formatted_models, nexus_name, location):
    """
    This function takes dict of models (from model_reader function), NEXUS
    filename (from fasta2nexus function) and a path to write Garli conf files.

    It composes the appropiate number of Garli configurations files (one per
    model) in the specified folder (in location). The conf file searchreps
    parameter is always 1 because we only want to calculate the likelihood for
    our dataset given certain evolution model.

    This function outputs a initial likelihood dict with model names as keys
    and 0 likelihood for values.

    """
    num_models = len(formatted_models)
    base_conf = "[general]\ndatafname = {0}\nconstraintfile = none\nstreefname = stepwise\nattachmentspertaxon = 50\nofprefix = {1}\nrandseed = -1\navailablememory = 512\nlogevery = 10\nsaveevery = 100\nrefinestart = 1\noutputeachbettertopology = 0\noutputcurrentbesttopology = 0\nenforcetermconditions = 1\ngenthreshfortopoterm = 20000\nscorethreshforterm = 0.05\nsignificanttopochange = 0.01\noutputphyliptree = 0\noutputmostlyuselessfiles = 0\nwritecheckpoints = 0\nrestart = 0\noutgroup = 1\nresampleproportion = 1.0\ninferinternalstateprobs = 0\noutputsitelikelihoods = 0\noptimizeinputonly = 0\ncollapsebranches = 1\n\nsearchreps = 1\nbootstrapreps = 0\n\n[model1]\n{2}\n\n"
    master = "[master]\nnindivs = 4\nholdover = 1\nselectionintensity = 0.5\nholdoverpenalty = 0\nstopgen = 5000000\nstoptime = 5000000\n\nstartoptprec = 0.5\nminoptprec = 0.01\nnumberofprecreductions = 10\ntreerejectionthreshold = 50.0\ntopoweight = 1.0\nmodweight = 0.05\nbrlenweight = 0.2\nrandnniweight = 0.1\nrandsprweight = 0.3\nlimsprweight =  0.6\nintervallength = 100\nintervalstostore = 5\nlimsprrange = 6\nmeanbrlenmuts = 5\ngammashapebrlen = 1000\ngammashapemodel = 1000\nuniqueswapbias = 0.1\ndistanceswapbias = 1.0"
    likelihood_dict = {}
    for counter, (model, model_block) in zip(range(1, num_models + 1),
        formatted_models.items()):
        conf_name = model + ".conf"
        with open(os.path.join(location,conf_name), "a") as conf_file:
            conf_string = base_conf.format(nexus_name, model, model_block) + \
            master
            conf_file.write(conf_string)
        print("Generated Garli config file for {0} model ({1}/{2})"
                .format(model, counter, num_models))
        likelihood_dict[conf_name] = 0
    return likelihood_dict


def get_likelihood(garli_stdout):
    """
    This helper function takes a stdout from a Garli run and properly extracts
    the likelihood value of the first (and only) search replicate in the run.
    It outputs the likelihood value as a float.

    """
    likelihood = float(garli_stdout.partition("Results:\nReplicate 1 : ")[2]
                    .partition("\n\nParameter estimates:")[0])
    return likelihood

# def compute_likelihoods(model_dict):
#    """
#    This function was intended to calculate all likelihoods in Garli from
#    a model dict. I commented out this function because it was written as a
#    serial task (slower)
#    """
#     num_models = len(model_dict)
#     for counter, model in zip(range(1, num_models + 1), model_dict):
#         print("Computing tree with for {0}. ({1}/{2})".format(model, counter,
#                 num_models))
#         garli_stdout, garli_stderr = run_cmd(['Garli',model], temp_directory)
#         likelihood_score = get_likelihood(garli_stdout)
#         model_dict[model] = likelihood_score
#         print("Likelihood value of {0} for {1}".format(str(likelihood_score),
#                 model.replace(".conf","")))
#     return model_dict

def compute_likelihood(input_tuple):
    """
    This function recieves a tuple containing three arguments (model name,
    number of task, total of tasks)
    It runs Garli from the proper conf file and outputs a tuple of model name
    and computed likelihood score using the previously defined helper functions

    """
    model, counter, num_models, ntax, nchars = input_tuple
    print("Computing tree with for {0}. ({1}/{2})".format(model, counter,
            num_models))
    garli_stdout, garli_stderr = run_cmd(['Garli',model], temp_directory)
    likelihood_score = get_likelihood(garli_stdout)
    k = compute_k(model, ntax)
    aic_score = compute_aic(likelihood_score, k)
    bic_score = compute_bic(likelihood_score, k, nchars)
    print("Likelihood value of {:.2f} for {}".format(likelihood_score,
            model.replace(".conf","")))
    print("AIC value of {:.2f} for {}".format(aic_score,
            model.replace(".conf","")))
    print("BIC value of {:.2f} for {}".format(bic_score,
            model.replace(".conf","")))
    return model, likelihood_score, aic_score, bic_score

def calculateParallel(model_dict, ntax, nchars, threads=jobs):
    """
    This function excecute all Garli runs in parallel given the number of
    threads defined in --jobs option.

    This function returns a full likelihoods dict (containing likelihood scores
    ) from a init likelihood dict (output from garliconf_gen function).

    """
    print("Computing using {} threads.".format(str(threads)))
    num_models = len(model_dict)
    counter_list = range(1, num_models + 1)
    total_task_list = [num_models] * num_models
    ntax_list = [ntax] * num_models
    nchars_list = [nchars] * num_models
    pool = ThreadPool(threads)
    likelihoods = pool.map(compute_likelihood, zip(model_dict.keys(), counter_list, total_task_list, ntax_list, nchars_list))
    pool.close()
    pool.join()
    likelihood_dict = {}
    best_l_pair = ["no model", -math.inf]
    best_aic_pair = ["no model", math.inf]
    best_bic_pair = ["no model", math.inf]
    for model, likelihood, aic, bic in likelihoods:
        likelihood_dict[model] = [likelihood, aic, bic]
        if likelihood > best_l_pair[1]:
            best_l_pair[0] = model
            best_l_pair[1] = likelihood
        if bic < best_bic_pair[1]:
            best_bic_pair[0] = model
            best_bic_pair[1] = bic
        if aic < best_aic_pair[1]:
            best_aic_pair[0] = model
            best_aic_pair[1] = aic
    print("\nModel with highest likelihood score is {}!".format(best_l_pair[0]))
    print("Best model selected according AIC is {}!".format(best_aic_pair[0]))
    print("Best model selected according BIC is {}!".format(best_aic_pair[0]))
    return likelihood_dict

def write_results(final_scores, result_file):
    """
    This function takes final_scores from calculateParallel function and saved
    them in a file (defined by --result_file option).

    This file is separated by tabs and it is of the form of:

    MODEL_NAME  LIKELIHOOD_SCORE AIC_SCORE BIC_SCORE

    """
    print("Writting results to {}.".format(result_file))
    with open(result_file, "a") as output_file:
        header = "model\t-log(L)\tAIC\tBIC\n"
        output_file.write(header)
        for model, scores in final_scores.items():
            l_score, aic_score, bic_score = scores
            result_string = \
            "{}\t{:.2f}\t{:.2f}\t{:.2f}\n".format(model.replace(".conf",""),
            l_score, aic_score, bic_score)
            output_file.write((result_string))
    return "File saved!"

def compute_k(model_string, ntax):
    params_per_model = {"GTR":8, "SYM":5, "HKY":4, "K80":1, "F81":3,
                        "JC":0}
    model_info = model_string.replace(".conf","").split("+")
    model_name = model_info[0]
    inv_status = 1 if "I" in model_info else False
    gamma_status = 1 if "G" in model_info else False
    n_branches = ntax * 2 - 3
    k = params_per_model[model_name] + inv_status + gamma_status + n_branches
    return k

def compute_aic(likelihood, k):
    aic = -2 * likelihood + 2 * k
    return aic

def compute_bic(likelihood, k, nchars):
    bic = -2 * likelihood + k * math.log(nchars)
    return bic

#get argvs from flags

original_fasta = options.fasta_file
model_file  = options.model_list
temp_directory = options.temp_directory
result_file    = options.result_file

#run the program using the functions defined above

print("\nInput file: {} in format FASTA".format(original_fasta))
print("Creating temporary dir named: {}\n".format(temp_directory))

run_cmd(['mkdir', temp_directory])

print("Converting FASTA file to NEXUS file format")
nexus_filename, ntax, nchars = fasta2nexus(original_fasta, temp_directory)
print("NEXUS file written to {}".format(nexus_filename))

print("Processed a NEXUS file containing {0} taxa and {1} chars.\n"
        .format(ntax, nchars))

print("Reading model list from {}".format(model_file))
model_dictionary = model_reader(model_file)

print("Generating Garli conf files for specified models")

likelihood_init = garliconf_gen(model_dictionary, nexus_filename,
        temp_directory)

print("\n\nComputing likelihood scores for {} models...\n"
        .format(len(likelihood_init)))

likelihood_scores = calculateParallel(likelihood_init, ntax, nchars)

print("\nLikelihood & information criteria calculations completed!\n")

print(write_results(likelihood_scores, result_file))

print("Done!\n")

#calculate execution time
end_time = time.monotonic()
time_delta = str(timedelta(seconds=end_time - start_time))

print("Task completed in: {}".format(time_delta))
print("... plus 8 hours of coding! >.< ")
