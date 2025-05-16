

command_line.py is the CLI wrapper for CellSnake. It handles argument parsing and forwards commands to Snakemake. 
For example, running cellsnake minimal data sends the request to the Snakefile. 
The Snakefile then builds a list of output files for Snakemake's rule all.

This list is generated using functions in create_plot_functions.py. That file lives in the workflow_utils folder. 
The same folder contains workflow_funcs.py, which detects sample names from the input directory. 
It also includes extra_functions.py, which supports the other two scripts.

Tests (Python side)
Tests for the Python logic are in tests/py-tests. Test file names match the script they check. 
For example, one test focuses on the scrna-celltypist.py script, which is used in the Snakemake workflow. 
It runs cell type predictions using the CellTypist library.

scripts folder (R side)
This folder holds the main computational logic for the workflow. 
Each rule in Snakemake usually maps to an R script here. Not all R scripts are refactored yet, but the important ones delegate their logic to helper files.
These helpers are inside the workflow_helper_functions folder inside scripts folder. 

Tests (R side)
R tests follow the same naming convention. They test the helper files directly. 
This is needed for the covr package, which checks test coverage in R. covr uses source() on the file being tested. 
If we tested full scripts, we’d accidentally run the entire workflow. 
That's why the logic had to be separated into helper files.

testData folder
This folder stores test outputs. It's essential for testing. 
The cellsnake test all command depends on it, especially on the datafolder and cmdTest inside.

Do not delete this 'data' folder. 
If you do, you can regenerate test data with generate_10X_data.R in the scripts folder. It requires real data, which is available from the CellSnake repository. The script expects samples named 10X_17_28 and 10X_17_29.




Key commands:
cellsnake test all: Runs all tests (Python and R). Acts like a full integration test.

cellsnake test covr: Runs test coverage analysis using covr_coverage.R.

cellsnake <mode> <input_dir>: Runs the workflow in a chosen mode (e.g., minimal, full, integration).

rules folder
This folder holds all .smk files for the Snakemake pipeline. Each rule defines:
-Inputs and outputs
-Expected parameters
-R script or shell command to run
-Rules cover everything from filtering to clustering to cell type annotation.



prototype/experimental:  persistent.sh

There is a bash script embedded inside each affected rule - those beloning to the minimal and standard one.
If you wish to try it out, be aware that it may crash or fail to work. More info is found inside the bash scrip.
you need to manually turn 'background' to TRUE to use its customized persistent r-session execution.

