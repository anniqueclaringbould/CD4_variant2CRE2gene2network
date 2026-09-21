## Get Hi-C interaction frequencies for all E-G pairs on the same chromosome

import hicstraw
import pandas as pd

# Define functions ---------------------------------------------------------------------------------

# look up contact frequency for one E-G pair
def get_int_freq_pair(enh_bin, tss_bin, matrix_chr):
  int_freq = matrix_chr.getRecords(int(enh_bin), int(enh_bin), int(tss_bin), int(tss_bin))
  if len(int_freq) > 0:
    int_freq = int_freq[0].counts
  else:
    int_freq = float("nan")
  return(int_freq)

# look up all contact frequencies for E-G pairs on one chromosome
def get_eg_int_freqs_chr(etp, hic, chrom, hic_res):
  matrix_chr =  hic.getMatrixZoomData(chrom, chrom, "observed", "SCALE", "BP", hic_res)
  etp_chr = etp.loc[(etp['pert_chr'] == chrom) & (etp['gene_chr'] == chrom)].copy()
  etp_chr['hic_int_freq'] = etp_chr.apply(
    lambda x: get_int_freq_pair(x['enh_hic_start'], x['tss_hic_start'], matrix_chr),
    axis = 1
  )
  return(etp_chr)

# get contact frequencies for all E-G pairs in a ETP table
def get_eg_int_freqs(etp, hic, hic_res):
  
  # calculate hic_bin for each enhancer based on enhancer center
  etp['enh_center'] = (etp['pert_start'] + etp['pert_end']) / 2
  etp['enh_hic_start'] = etp['enh_center'] // hic_res * hic_res
  etp['enh_hic_end'] = etp['enh_hic_start'] + hic_res
  
  # calculate hic_bin for each gene TSS
  etp['tss_hic_start'] = etp['gene_tss'] // hic_res * hic_res
  etp['tss_hic_end'] = etp['tss_hic_start'] + hic_res
  
  # create list for output and get all chromosomes in ETP table
  etp_hic = list()
  chroms = etp['pert_chr'].unique()
  
  # get interaction frequencies for all E-G pairs on the same chromosome (exclude trans-pairs)
  print("Getting interaction frequencies for:")
  for chrom in chroms:
    print("  " + chrom)
    etp_hic.append(get_eg_int_freqs_chr(etp, hic, chrom, hic_res))
  etp_hic = pd.concat(etp_hic)
  
  # set interaction frequency of trans interactions to NaN and add to dataframe
  etp_trans = etp.loc[etp['pert_chr'] != etp['gene_chr']].copy()
  etp_trans['hic_int_freq'] = float("nan")
  etp_hic = pd.concat([etp_hic, etp_trans])
  
  return(etp_hic)

# Get contact frequencies for E-G pairs ------------------------------------------------------------

# load ETP table
etp = pd.read_csv(snakemake.input[0], low_memory = False)

# open connection to .hic file
hic = hicstraw.HiCFile(snakemake.input[1])

# get interaction frequencies for all E-G pairs in ETP table and write output to file
etp = get_eg_int_freqs(etp, hic, hic_res = snakemake.params[0])
etp = etp.drop(columns=['enh_center', 'enh_hic_start', 'enh_hic_end', 'tss_hic_start', 'tss_hic_end'])
etp.to_csv(snakemake.output[0], index = False)
