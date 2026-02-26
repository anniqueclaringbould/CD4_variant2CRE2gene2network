## Add samples, optional cell numbers and read files per sample to snakemake config object

import os
import re
import pandas as pd

# Define functions ---------------------------------------------------------------------------------

# function to get file path without file extension and read identifier
def clean_path(filepath):
    filepath = re.sub(r'\.(fastq|fq)(\.gz)?$', '', filepath)
    filepath = re.sub(r'_R\d+$', '', filepath)
    return(filepath)
  
# create a unique run id for one group of files
def create_run_id(group):
    fileset_codes = pd.factorize(group['fileset'])[0] + 1
    return group['sample'] + "_run" + fileset_codes.astype(str)
  
# create unique run id for each pair of files, assuming files from different reads only differ by a
# _R1, _R2 extension and are located in the same directory
def add_run_id(read_files):
  read_files['fileset'] = read_files['file'].apply(clean_path)
  read_files['run'] = read_files.groupby('sample', group_keys=False)[['sample', 'fileset']].apply(create_run_id)
  return(read_files[['sample', 'file', 'run']])

# Load samples and read files ----------------------------------------------------------------------

# load samples file
samples_file = pd.read_csv(config['samples_file'], sep='\t', header=None, names=['sample', 'cells'],
                           on_bad_lines='skip')

# create simple array of sample ids
samples = list(samples_file['sample'])

# create dictionary with cell numbers per sample (NaN) if not provided
cell_numbers = samples_file.set_index('sample')['cells'].to_dict()

# load reads file
read_files = pd.read_csv(config['reads_file'], sep='\t')

# create run id if not already present
if 'run' not in read_files.columns:
    read_files = add_run_id(read_files)

# create symlink path based on run id and read file
read_files['symlink'] = read_files['run'] + "_" + read_files['file'].apply(os.path.basename)

# create a dictionary containing the original file for each symlink
symlinks = read_files[['symlink', 'file']].set_index('symlink')['file'].to_dict()

# create a dictionary of lists, containing all read files (symlinks) per sample
read_files = read_files[['sample', 'symlink']].groupby('sample')['symlink'].apply(list).to_dict()

# add samples, cell numbers and read files to snakemake config object
config["samples"] = samples
config["cell_numbers"] = cell_numbers
config["symlinks"] = symlinks
config["read_files"] = read_files
