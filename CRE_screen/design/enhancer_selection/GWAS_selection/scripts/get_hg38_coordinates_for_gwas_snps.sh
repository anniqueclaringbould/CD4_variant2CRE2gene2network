#From: https://www.biostars.org/p/9534406/
ROOT=/g/steinmetz/project/otar
GWAS=$ROOT/data/gwas_selection

TRAITS=$(cat $GWAS/disease_names.tsv)

for TRAIT in $TRAITS
do

echo "
#####################################
NOW RUNNING: $TRAIT
#####################################
"

#header
echo -e "SNP\tChrom\tBP" > $GWAS/$TRAIT/snplist_pos.tsv

#for each SNP, look up the genome build 38 position and output that here into the format requested by Blagoje's pipeline
cat $GWAS/$TRAIT/snplist.tsv | while read rsid ;
do
  pos=$(curl -sX GET "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id=$rsid&retmode=text&rettype=text" | \
    sed 's/<\//\n/g' | \
    grep -o -P '\<CHRPOS\>.{0,15}' | \
    cut -f2 -d">") ;
  chr=$(echo $pos | cut -f1 -d":")
  bp=$(echo $pos | cut -f2- -d":") 
  echo -e "${rsid}""\tchr""${chr}""\t""${bp}" 
done >> $GWAS/$TRAIT/snplist_pos.tsv

done
