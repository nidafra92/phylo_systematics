#!/usr/bin/env python3


#imports
import sys, os
from optparse import OptionParser
from Bio import SeqIO
from Bio.Alphabet import IUPAC

exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
sys.path.insert(0,os.path.abspath(os.path.join(exe_path,"..", "..")))
exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]

parser = OptionParser(version="%prog 1.0", usage = "\n\n  %prog --fasta_file \
                    FASTA --model_list MODEL_LIST")
parser.add_option("-f", "--fasta_file", dest="fasta_file",
                  metavar="FILE",
                  type="string",
                  help="Name and path of input alignment in FASTA format. ")
parser.add_option("-m", "--model_list", dest="model_list",
                  metavar="FILE",
                  type="string",
                  #default="all",
                  help="Name and path of file contanining models to be tested. ")
parser.add_option("--result_file", dest="result_file",
                  metavar="FILE",
                  type="string",
                  default="final.phylo_results",
                  help="Name and path of where the files should end up. ")
parser.add_option("-t", "--temp_directory", dest="temp_directory",
                  metavar="PATH",
                  type="string",
                  default="tmp",
                  help="Name and path of temporary directory where calculations should take place. ")
parser.add_option("-j", "--jobs", dest="jobs",
                  default=1,
                  metavar="jobs",
                  type="int",
                  help="Specifies the number of jobs (operations) to run \
                  in parallel.")

(options, remaining_args) = parser.parse_args()

original_fasta = options.fasta_file
model_file  = options.model_list
temp_directory = options.temp_directory
result_file    = options.result_file




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
    gamma_dict = {True:"ratehetmodel = gamma\nnumratecats = 4\n", False:"ratehetmodel = none\nnumratecats = 1"}
    inv_dict = {True:"invariantsites = estimate\n",False:"invariantsites = none"}
    # parse model name
    model_info = model_string.split("+")
    model = model_info[0]
    inv_status = True if "I" in model_info else False
    gamma_status = True if "G" in model_info else False
    # compose Garli block
    composed_string = model_dict[model] + gamma_dict[gamma_status] + inv_dict[inv_status]
    return composed_string

def model_reader(model_file):
    """
    Reads a model file list and returns models with commands
    """
    models_to_eval = {}
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
    with open(fasta_filename) as fasta, open(os.path.join(location, nex_filename),"a") as nexus:
        sequences = SeqIO.parse(fasta, "fasta", alphabet=IUPAC.ambiguous_dna)
        new_seqs = []
        for sequence in sequences:
            new_name = name_reformat(sequence.description)
            sequence.id = new_name
            new_seqs.append(sequence)
        SeqIO.write(new_seqs, nexus, "nexus")
    return nex_filename

def garliconf_gen():
    pass

#print(exe_path)

os.mkdir(temp_directory)

#print(os.path.join(temp_directory, original_fasta.replace(".fasta",".nexus")))

nexus_filename = fasta2nexus(original_fasta, temp_directory)
print(nexus_filename)
