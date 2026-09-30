#Submit python script for CHEERS normalisation
#https://github.com/TrynkaLab/CHEERS

conda activate cheers

#Input info
ROOT=/g/steinmetz/project/otar
DIR=$ROOT/select_enhancer_peaks
SCRIPTS=$ROOT/scripts
GWAS=$ROOT/data/gwas_selection
PEAKS=$ROOT/data/trynka_lab/CHEERS_peak_files/PEAKS

#Make output directory
#mkdir $DIR/normalised_peaks
mkdir $DIR/normalised_peaks_across_conditions

#Run script for ATAC and H3K27ac peaks separately
for TYPE in ATAC H3K27ac

do

#select the stimulation conditions we want to include
#change name of files that we're including so the python script recognises them
#they need to have "ReadsInPeaks.txt" at the end
#FILES=$(ls $PEAKS/$TYPE/*naive*TH0*)
#FILES=$(ls $PEAKS/$TYPE/*naive* | grep -E '(16H_UNS|D5_TH0|16H_TH0)'\.txt$)
FILES=$(ls $PEAKS/$TYPE/)

for FILE in $FILES
do
cp $PEAKS/$TYPE/$FILE $PEAKS/$TYPE/${FILE}_ReadsInPeaks.txt
done

#Run CHEERS normalisation

#--input: path to files with read counts per peak (default: None)
#--prefix: file prefix 
#--outdir: directory where to output results
#python $SCRIPTS/CHEERS_normalize.py --input $PEAKS/$TYPE --prefix $TYPE --outdir $DIR/normalised_peaks/
python $SCRIPTS/CHEERS_normalize.py --input $PEAKS/$TYPE --prefix $TYPE --outdir $DIR/normalised_peaks_across_conditions/

done
