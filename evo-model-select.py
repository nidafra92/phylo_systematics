#!/usr/bin/env python3
description = \
"""
    
    %prog --fasta_file FASTA --model_list MODEL_LIST

    Code for a simple model tester for ML analysis using Garli for
    likelihood calculations.

    version 1.0
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
"""

#imports
import sys, os
from optparse import OptionParser
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import subprocess
import time
from datetime import timedelta

start_time = time.monotonic() # I want to record execution time.


exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
sys.path.insert(0,os.path.abspath(os.path.join(exe_path,"..", "..")))
exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]

parser = OptionParser(version="%prog 1.0",
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
parser.add_option("--result_file", dest="result_file",
                    metavar="FILE",
                    type="string",
                    default="final.phylo_results",
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

if not options.fasta_file:
    parser.error("\n\n\tMissing parameter --fasta_file FILE\n\n")
if not options.model_list:
    parser.error("\n\n\tMissing parameter --model_list FILE\n\n")

def run_cmd(cmd_str, work_path=None):
    """
    Throw exception if run command fails
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
    Reads a model file list and returns models with commands
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
    safe_name = "_".join(name.split("|")[-1].split()[0:4])
    return safe_name

def fasta2nexus(fasta_filename, location):
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
        print("Generated Garli config file for {0} model \
                ({1}/{2})".format(model, counter, num_models))
        likelihood_dict[conf_name] = 0
    return likelihood_dict


def get_likelihood(garli_stdout):
    likelihood = float(garli_stdout.partition("Results:\nReplicate 1 : ")[2]
                    .partition("\n\nParameter estimates:")[0])
    return likelihood

def compute_likelihoods(model_dict):
    num_models = len(model_dict)
    for counter, model in zip(range(1, num_models + 1), model_dict):
        print("Computing tree with for {0}. ({1}/{2})".format(model, counter,
                num_models))
        garli_stdout, garli_stderr = run_cmd(['Garli',model], temp_directory)
        likelihood_score = get_likelihood(garli_stdout)
        model_dict[model] = likelihood_score
        print("Likelihood value of {0} for {1}".format(str(likelihood_score),
                model.replace(".conf","")))
    return model_dict


def write_results(final_scores, result_file):
    print("Writting results to {}.".format(result_file))
    with open(result_file, "a") as output_file:
        for model, score in final_scores.items():
            result_string = \
            "{0}\t{1}\n".format(model.replace(".conf",""),str(score))
            output_file.write((result_string))
    return "File saved!"


original_fasta = options.fasta_file
model_file  = options.model_list
temp_directory = options.temp_directory
result_file    = options.result_file


print("\nInput file: {} in format FASTA".format(original_fasta))
print("Creating temporary dir named: {}\n".format(temp_directory))

run_cmd(['mkdir', temp_directory])

print("Converting FASTA file to NEXUS file format")
nexus_filename, ntax, nchars = fasta2nexus(original_fasta, temp_directory)
print("NEXUS file written to {}".format(nexus_filename))

print("Processed a NEXUS file containing {0} taxa and {1} \
        chars.\n".format(ntax, nchars))

print("Reading model list from {}".format(model_file))
model_dictionary = model_reader(model_file)

print("Generating Garli conf files for specified models")

likelihood_init = garliconf_gen(model_dictionary, nexus_filename,
        temp_directory)

print("\n\nComputing likelihood scores for {} \
        models...\n".format(len(likelihood_init)))

likelihood_scores = compute_likelihoods(likelihood_init)

print("\nLikelihood calculations completed!\n")

print(write_results(likelihood_scores, result_file))

print("Done!\n")

end_time = time.monotonic()
time_delta = str(timedelta(seconds=end_time - start_time))

print("Task completed in: {}".format(time_delta))
print("... plus 5 hours of coding! >.< ")
