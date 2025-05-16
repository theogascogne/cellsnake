#!/usr/bin/env python
'''
Created on 02/12/2022
cellsnake main
@author: Sinan U. Umu, sinanugur@gmail.com
'''

# =============================================================================
# cellsnake CLI entrypoint
#
# Usage:
#   This script is invoked by users via the `cellsnake` CLI to run various
#   scRNA-seq workflows (e.g. minimal, standard, integrated). It parses CLI
#   arguments, validates input, builds a snakemake command, and optionally
#   logs the run.
#   Example of a command: cellsnake standard <input>
#   Please note that the input has to be named 'data' and follow a file path pattern.
#   The provided test sample in the Readme is an example of this
#   Or the data folder within cellsnake's test directory
# =============================================================================

#from __future__ import print_function
import re
import warnings
warnings.filterwarnings("ignore")
from docopt import docopt
import os
import sys
import subprocess
import shutil
import datetime
import random
#from fuzzywuzzy import fuzz
import timeit
import errno
import os
import yaml
from yaml.loader import SafeLoader
from subprocess import call
import pathlib
from pathlib import Path

#from schema import Schema, And, Or, Use, SchemaError

from collections import defaultdict

import cellsnake
cellsnake_path=os.path.dirname(cellsnake.__file__)
options = ["clustree","clusteringTree","minimal","standard","advanced", "test"] #and integration


__author__ = 'Sinan U. Umu'
__version__= '0.2.0.12'
__logo__="""
             _  _                     _           
            | || |                   | |          
  ___   ___ | || | ___  _ __    __ _ | | __  ___  
 / __| / _ \| || |/ __|| '_ \  / _` || |/ / / _ \ 
| (__ |  __/| || |\__ \| | | || (_| ||   < |  __/ 
 \___| \___||_||_||___/|_| |_| \__,_||_|\_\ \___| 
                                                  
"""  


__licence__="""
MIT License
Copyright (c) 2023 Sinan U. Umu (SUU) sinanugur@gmail.com
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

#================================================
# To add more functionality to cellsnake, please
# add them to this list below with a description
#================================================

__doc__=f"""Main cellsnake executable, version: {__version__}
{__logo__} 
Usage:
    cellsnake <command> <INPUT> [options] [--unlock|--remove] [--dry]
    cellsnake integrated <command> <INPUT> [options] [--unlock|--remove] [--dry]
    cellsnake --generate-template
    cellsnake --install-packages
    cellsnake (-h | --help)
    cellsnake --version

commands:
    minimal                                Run cellsnake with minimal workflow. 
    standard                               Run cellsnake with standard workflow.
    advanced                               Run cellsnake with advanced workflow.
    clustree                               Run cellsnake with clustree workflow.
    integrate                              Run cellsnake to integrate samples under analyses folder.
                                           This option expects you have already finished processing multiple samples.

main arguments:
    INPUT                                  Input directory or a file to process (if a directory given, batch mode is ON).
                                           After integration, provide the Seurat object file not a directory to work on
                                           integrated object, see tutorial if necessary (e.g. analyses_integrated/seurat/integrated.rds)
    --configfile <text>                    Config file name in YAML format, for example, "config.yaml". No default but can be created with --generate-template.
    --metadata <text>                      Metadata file name in CSV, TSV or Excel format, for example, "metadata.csv", header required, first column sample name. No default but can be created with --generate-template.
    --metadata_column <text>               Metadata column for differential expression analysis [default: condition].

other arguments:
    --gene <gene or filename>              Create publication ready plots for a gene or a list of genes from a text file.
    
main options:
    --percent_mt <double>                  Maximum mitochondrial gene percentage cutoff, 
                                           for example, 5 or 10, write "auto" for auto detection [default: 10].
    --resolution <double>                  Resolution for cluster detection, write "auto" for auto detection [default: 0.8].

other options:
    --doublet_filter <bool>                [default: True] #this may fail on some samples and on low memory. If you have a problem, try False.
    --percent_rp <double>                  [default: 0] #Ribosomal genes minimum percentage (0-100), default no filtering
    --min_cells <integer>                  [default: 3] #seurat default, recommended
    --min_features <integer>               [default: 200] #seurat default, recommended, nFeature_RNA
    --max_features <integer>               [default: Inf] #seurat default, nFeature_RNA, 5000 can be a good cutoff
    --min_molecules <integer>              [default: 0] #seurat default, nCount_RNA, min_features usually handles this so keep it 0
    --max_molecules <integer>              [default: Inf] #seurat default, nCount_RNA, to filter potential doublets, doublet filtering is already default, so keep this Inf
    --highly_variable_features <integer>   [default: 2000] #seurat defaults, recommended
    --variable_selection_method <text>     [default: vst] #seurat defaults, recommended
    
    --normalization_method <text>          [default: LogNormalize]
    --scale_factor <integer>               [default: 10000]
    --logfc_threshold <double>             [default: 0.25]
    --test_use <text>                      [default: wilcox]
    
    --persist <text>                       [default: False] #An experimental feature, runs the workflow on a persistent R session. It is unstable and only works for Minimal workflow.


    --mapping <text>                       [default: org.Hs.eg.db] #you may install others from Bioconductor, this is for human
    --organism <text>                      [default: hsa] #alternatives https://www.genome.jp/kegg/catalog/org_list.html
    --species <text>                       [default: human] for cellchat, #only human or mouse is accepted

plotting parameters:
    --min_percentage_to_plot <double>        [default: 5] #only show clusters more than % of cells on the legend
    --show_labels <bool>                     [default: True] #
    --marker_plots_per_cluster_n <integer>   [default: 20] #plot summary marker plots for top markers
    --umap_markers_plot <bool>               [default: True]
    --tsne_markers_plot <bool>               [default: False]

annotation options:
    --singler_ref <text>                    [default: BlueprintEncodeData] # https://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html#1_Overview
    --celltypist_model <text>               [default: Immune_All_Low.pkl] #refer to Celltypist for another model 

microbiome options:
    --kraken_db_folder <text>              No default, you need to provide a folder with kraken2 database
    --taxa <text>                          [default: genus] # available options "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"
    --microbiome_min_cells <integer>       [default: 1]
    --microbiome_min_features <integer>    [default: 3]
    --confidence <double>                  [default: 0.05] #see kraken2 manual
    --min_hit_groups <integer>             [default: 4] #see kraken2 manual

integration options:
    --dims <integer>                       [default: 30] #refer to Seurat for more details
    --reduction <text>                     [default: cca] #refer to Seurat for more details

others:
    --generate-template                    Generate config file template and metadata template in the current directory.
    --install-packages                     Install, reinstall or check required R packages.
    -j <integer>, --jobs <integer>         Total CPUs. [default: 1]
    -u, --unlock                           Rescue stalled jobs (Try this if the previous job ended prematurely or currently failing).
    -r, --remove                           Delete all output files (this won't affect input files).
    -d, --dry                              Dry run, nothing will be generated.
    -h, --help                             Show this screen.
    --version                              Show version.
    --lint                                 Linting.
    --stats                                Show workflow metrics.
"""

#=============================
# Function: check_command_line_arguments(arguments)
#
# Purpose:
#   Validates the command-line arguments passed to the `cellsnake` CLI tool.
#   This function ensures required files and options are correctly set up before
#   the pipeline begins execution. Returns False on invalid inputs.
#
#   Ensures an <INPUT> path exists unless the command is 'test'
#   If the 'integrated' option is set, confirms the input is a valid Seurat .rds file
#   Verifies that required files like --configfile and --metadata exist
#   If --kraken_db_folder is specified, confirms it exists and has the expected structure
#   Checks that --taxa, if given, is one of the accepted taxonomy levels
#
# Returns:
#   True if all validations pass, otherwise False (with a descriptive error message)
#=============================


def check_command_line_arguments(arguments):
    # Ensure <INPUT> exists, but only if it's not a 'test' mode
    if not arguments.get("<command>") == "test" and not os.path.exists(arguments["<INPUT>"]):
        print("File or input directory not found: ", arguments["<INPUT>"])
        return False

    # Validate 'integrated' option with directory input
    if arguments.get("integrated") and os.path.isdir(arguments["<INPUT>"]):
        print("You are running the integrated option, but you provided a directory, not a Seurat object file!")
        print("The default Seurat object is usually here: analyses_integrated/seurat/integrated.rds")
        print("You can try something like: cellsnake integrated standard analyses_integrated/seurat/integrated.rds")
        return False
    
    # Validate 'integrated' with invalid file type
    if arguments.get("integrated") and os.path.isfile(arguments["<INPUT>"]):
        file_extension = pathlib.Path(arguments["<INPUT>"])
        if file_extension.suffix.lower() != ".rds":
            print("You are running the integrated option, but the file is not a Seurat object file (.rds)!")
            print("The default Seurat object is usually here: analyses_integrated/seurat/integrated.rds")
            print("You can try something like: cellsnake integrated standard analyses_integrated/seurat/integrated.rds")
            return False

    # Validate config file
    if arguments.get("--configfile") and not os.path.isfile(arguments["--configfile"]):
        print("Config file not found: ", arguments["--configfile"])
        return False

    # Validate metadata file
    if arguments.get("--metadata") and not os.path.isfile(arguments["--metadata"]):
        print("Metadata file not found: ", arguments["--metadata"])
        return False

    # Validate Kraken DB folder
    if arguments.get("--kraken_db_folder"):
        kraken_db_folder = Path(arguments["--kraken_db_folder"])  # convert here
        if not kraken_db_folder.exists() and not (kraken_db_folder / "inspect.txt").is_file():
            print("KrakenDB directory not found: ", kraken_db_folder)
            print("You should download a proper DB from this link: https://benlangmead.github.io/aws-indexes/k2")
            return False

    # Validate taxa level
    valid_taxa = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    if arguments.get("--taxa") not in valid_taxa:
        print("Select a correct taxa level for microbiome analysis:", arguments["--taxa"])
        print("Possible options: ", valid_taxa)
        return False

    return True

#=============================
#
# CommandLine builds the full Snakemake command for running the Cellsnake pipeline.
#
# - Initializes with a base command and a random run ID.
# - Loads config from --configfile (if provided) or from command-line arguments.
# - Assembles Snakemake options:
#     CPU count (--jobs), Snakefile path, key=value config pairs
#     Special flags: --dry, --unlock, --remove
# - Handles integration-specific options and optional inputs (e.g., --gene, --kraken_db_folder).
# - Writes a log of the full command, parameters, and runtime.
#
# Example output:
# snakemake --rerun-incomplete -k -j 8 -s /path/to/Snakefile \
#   --config datafolder=my_data/ gene_to_plot=TP53 runid=abc12345 ...
#
#=============================

class CommandLine:
    def __init__(self):
        self.snakemake = "snakemake --rerun-incomplete -k "  # Base Snakemake command
        self.runid = "".join(random.choices("abcdefghisz", k=3) + random.choices("123456789", k=5))  # Random run ID
        self.config = []                         # Holds config key=value pairs
        self.configfile_loaded = False           # Flag to check if config file was loaded
        self.is_integrated_sample = False        # Whether sample is part of an integration
        self.is_this_an_integration_run = False  # Whether this is an integration run
        self.parameters = dict()                 # Stores all run parameters
        self.log = True                          # Logging enabled unless overridden by dry/unlock/remove/test

    def __str__(self):
        return self.snakemake

    def __repr__(self):
        return self.snakemake

    def add_config_argument(self):
        # Add accumulated config arguments to the snakemake command
        self.snakemake += " --config " + " ".join(self.config)

    def load_configfile_if_available(self, arguments):
        # Load YAML config file if provided and not already loaded
        if not self.configfile_loaded and arguments["--configfile"]:
            self.snakemake += f" --configfile={arguments['--configfile']}"
            configfile = arguments["--configfile"]
            with open(configfile) as f:
                self.parameters = yaml.load(f, Loader=SafeLoader)
            self.configfile_loaded = True

    def prepare_arguments(self, arguments):
        # Build Snakemake command line using arguments
        self.snakemake += f" -j {arguments['--jobs']} "  # Number of cores
        self.snakemake += f" -s {cellsnake_path}/scrna/workflow/Snakefile "  # Path to Snakefile

        self.load_configfile_if_available(arguments)  # Load config file if specified

        # Add input folder unless this is an integration run
        if not self.is_this_an_integration_run:
            self.config.append(f"datafolder={arguments['<INPUT>']}")

        self.config.append(f"cellsnake_path={cellsnake_path}/scrna/")  # Add pipeline path

        # Exclude these arguments from being passed as config key=value
        excluded_args = [
            "--jobs", "integrated", "--configfile", "--option", "--gene", "--kraken_db_folder",
            "--unlock", "--remove", "--dry", "--help", "--version", "<INPUT>", "<command>",
            "--install-packages", "--generate-template"
        ]

        # Handle general parameters and convert them to config entries
        for param, value in arguments.items():
            if param not in excluded_args:
                key = param.lstrip("--")
                if not self.configfile_loaded:
                    self.config.append(f"{key}={value}")
                    self.parameters[key] = str(value)
                elif self.parameters.get(key) and key not in sys.argv:
                    self.config.append(f"{key}={self.parameters.get(key)}")
                else:
                    self.config.append(f"{key}={value}")
                    self.parameters[key] = str(value)

        self.config.append(f"runid={self.runid}")  # Add random run ID

        # Handle gene input: either a file or single gene name
        if arguments["--gene"]:
            if os.path.isfile(arguments["--gene"]):
                self.config.append(f"selected_gene_file={arguments['--gene']}")
            else:
                self.config.append(f"gene_to_plot={arguments['--gene']}")

        # Handle Kraken DB path
        if arguments["--kraken_db_folder"]:
            self.config.append(f"kraken_db_folder={arguments['--kraken_db_folder']}")

        # Integration sample flag
        if self.is_integrated_sample:
            self.config.append("is_integrated_sample=True")

        # Add operation type (standard or integration)
        if not self.is_this_an_integration_run:
            self.config.append(f"option={arguments['<command>']}")
        else:
            self.config.append("option=integration")

        # Handle special Snakemake modes
        if arguments["--dry"]:
            self.snakemake += " -n "
            self.log = False
        if arguments["--unlock"]:
            self.snakemake += " --unlock "
            self.log = False
        if arguments["--remove"]:
            self.snakemake += " --delete-all-output "
            self.log = False
        if arguments["<command>"] == "test":
            self.config.append(f"test_mode={arguments['<INPUT>']}")
            self.log = False

        self.add_config_argument()  # Finalize config arguments

    def write_to_log(self, start):
        # Log run information and duration to a timestamped file
        logname = "_".join(["cellsnake", self.runid, datetime.datetime.now().strftime("%y%m%d_%H%M%S"), "runlog"])
        stop = timeit.default_timer()
        if self.log:
            with open(logname, "w") as f:
                f.write(__logo__ + "\n")
                f.write(f"Run ID : {self.runid}\n")
                f.write(f"Cellnake version : {__version__}\n")
                f.write(f"Cellsnake arguments : {' '.join(sys.argv)}\n\n")
                f.write("------------------------------\n")
                f.write(f"Snakemake arguments : {self.snakemake}\n\n")
                f.write("------------------------------\n")
                f.write(f"Run parameters: {self.parameters}\n\n")
                f.write("Total run time: {:.2f} mins\n".format((stop - start) / 60))
                
#This function, run_integration, is currently not in use.
def run_integration(arguments):

    start = timeit.default_timer()
    snakemake_argument=CommandLine()
    snakemake_argument.is_this_an_integration_run = True
    snakemake_argument.prepare_arguments(arguments)
    subprocess.check_call(str(snakemake_argument),shell=True)
    snakemake_argument.write_to_log(start)


    #then run workflow on integrated dataset
    #snakemake_argument=CommandLine()
    #snakemake_argument.is_this_an_integration_run = False
    #snakemake_argument.is_integrated_sample = True
    #snakemake_argument.config.append("datafolder=analyses_integrated/seurat/integrated.rds")
    #try:
    #    arguments["--option"].remove("integration")
    #except:
    #    pass
    #snakemake_argument.prepare_arguments(arguments)
    #subprocess.check_call(str(snakemake_argument),shell=True)
    #return snakemake_argument
    
# ======================================
#   run_workflow(arguments)

#   Given parsed CLI arguments:
#   Initializes CommandLine instance and sets flags.
#   Prepares snakemake command with all configs.
#   Executes command using os.system().
#   Checks for execution errors and writes run logs if successful.
# =======================================

def run_workflow(arguments):
    start = timeit.default_timer()  # Track start time
    snakemake_argument = CommandLine()  # Initialize command line builder

    if arguments["integrated"]:
        snakemake_argument.is_integrated_sample = True  # Mark as integrated sample

    snakemake_argument.prepare_arguments(arguments)  # Construct full command

    # Convert command object to string and run using shell
    command = str(snakemake_argument)
    exit_code = os.system(command)

    # Handle failed run
    if exit_code != 0:
        print(f"Error: Command failed with exit code {exit_code}")
        return

    # Log successful run
    snakemake_argument.write_to_log(start)
    
# =============================
# main()
#
#   Parses CLI arguments with docopt.
#   Handles special commands/flags:
#     -test modes (all, cover)
#     -lint check on Snakefile
#     -generate template config and metadata files
#     -install required R packages
#   Validates arguments.
#   Runs the workflow for valid commands.
# =============================

def main():
    cli_arguments = docopt(__doc__, version=__version__)
    
    # Handle test mode input "cellsnake test all"
    if cli_arguments['<command>'] == "test":
        if cli_arguments["<INPUT>"] is None:
            print("Please provide a test mode (e.g., `cellsnake test all` or `cellsnake test cover`)")
            return

        mode = cli_arguments["<INPUT>"]
        if mode not in ["all", "cover"]:
            print(f"Unknown test mode/input '{INPUT}', It is not supposed to be a directory if that was used , please use 'all' or 'cover'.")
            return
        print(f"Running test in '{mode}' mode...")
    
    # Run Snakemake's linting system
    if cli_arguments.get("--lint"):
        snakefile_path = f"{cellsnake_path}/scrna/workflow/Snakefile"
        result = os.system(f"snakemake --lint -s {snakefile_path}")
        if result == 0:
            print("Snakemake lint check passed.")
        else:
            print("Lint check found issues.")
        return
    
    # Generate template config and metadata files
    if cli_arguments["--generate-template"]:
        print("Generating config.yaml file...")
        print("You can use this as a template for a cellsnake run. You may change the settings.")
        shutil.copyfile(cellsnake_path + "/scrna/config.yaml", 'config.yaml')

        print("Generating metadata.csv file...")
        with open("metadata.csv", "w") as f:
            f.write("sample,condition\n")
            f.write("sample1,condition1\n")
            f.write("sample2,condition2\n")
        return

    if cli_arguments["--install-packages"]:
        subprocess.check_call(cellsnake_path + "/scrna/workflow/scripts/scrna-install-packages.R")
        return

    if not check_command_line_arguments(cli_arguments):
        print("""Please check your command line arguments. Use "cellsnake --help" for more information""")
        return

    # List of valid commands
    valid_commands = ['minimal', 'standard', 'advanced', 'clustree', 'integrate', 'test']
    
    # If the command is in the valid commands list, run the workflow
    if cli_arguments['<command>'] in valid_commands:
        run_workflow(cli_arguments)
