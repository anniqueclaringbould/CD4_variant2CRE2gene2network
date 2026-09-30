#Submit python script for CHEERS LD block formation
#https://github.com/TrynkaLab/CHEERS

conda activate cheers
ml tabix

ROOT=/g/steinmetz/project/otar
SCRIPTS=$ROOT/scripts
GWAS=$ROOT/data/gwas_selection
LD=$ROOT/data/trynka_lab/LD_GRCH38

TRAITS=$(cat $GWAS/disease_names.tsv)

for TRAIT in $TRAITS

do

echo "
#####################################
NOW RUNNING: $TRAIT
#####################################
"

mkdir $GWAS/$TRAIT/LD

cd $GWAS/$TRAIT

python $SCRIPTS/create_LD_blocks.py $GWAS/$TRAIT/snplist_pos.tsv $GWAS/$TRAIT/LD $LD > $GWAS/$TRAIT/create_LD.log

done
