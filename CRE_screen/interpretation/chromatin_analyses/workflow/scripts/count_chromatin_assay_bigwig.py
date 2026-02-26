## Compute assay signal at elements. Based on ABC code: https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/main/workflow/scripts/neighborhoods.py

import argparse
import pandas as pd
import pyBigWig as bw

# Define funtions ----------------------------------------------------------------------------------

# function to get total from bigwig file
def count_total_bigwig(assay_bigwig):
  assay = bw.open(assay_bigwig)
  result = sum(
    l * assay.stats(ch, 0, l, 'mean', exact = True)[0] for ch, l in assay.chroms().items()
  )
  return result

# function to count assay signal (sum) at elements in bed file and write output to new file
def count_signal_bigwig(elements_bed, assay_bigwig, output_file):
  
  # count total signal in assay file
  total = count_total_bigwig(assay_bigwig)
  
  # open elements bed file and assay bigwig file
  elements = pd.read_table(elements_bed, names = ['chr', 'start', 'end'], header = None)
  assay = bw.open(assay_bigwig)
  
  # count signal at every enhancer
  with open(output_file, 'wb') as outfile:
    header = ('\t'.join(['chr', 'start', 'end', 'signal.count', 'signal.cpm']) + '\n').encode('ascii')
    outfile.write(header)
    for chr, start, end, in elements.itertuples(index = False, name = None):
      try:
        if chr not in assay.chroms():
          signal = 0
          signal_norm = 0
        else:
          signal = (assay.stats(chr, int(start), int(end), type = 'sum', exact = True)[0] or 0)
          signal_norm = signal * 1e6 / total
      except RuntimeError:
        print('Failed on', chr, start, end)
        raise
      output = ('\t'.join([chr, str(start), str(end), str(signal), str(signal_norm)]) + '\n').encode('ascii')
      outfile.write(output)

# count assay signal at elements from provided files if run from command line ----------------------

if __name__ == "__main__":

  # parse command line arguments
  parser = argparse.ArgumentParser(description = ("Count assay signal at elements"))
  parser.add_argument("-e", "--elements_bed", type = str, required = True,
                      help = "Bed file containing elements")
  parser.add_argument("-a", "--assay_bigwig", type = str, required = True,
                      help = "BigWig containing assay signal")
  parser.add_argument("-o", "--output_file", type = str, required = True,
                      help = "Output file path")
  args = parser.parse_args()

  # count assay signal from bigwig at elements in bed file
  count_signal_bigwig(elements_bed = args.elements_bed, assay_bigwig = args.assay_bigwig,
                      output_file = args.output_file)
