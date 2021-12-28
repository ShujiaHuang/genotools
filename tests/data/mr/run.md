```bash
python ../../../genotools/mr.py TTC -I test2.vcf.gz --fam test.fam| awk '$1!~/^##/ {print $1 "\t" $2 "\t" $10 "\t" $11 "\t" $12}' | awk 'BEGIN{FS="[\t:]"; OFS="\t";}{if(NR==1){print $0}else{print $1,$2,$3,$6,$9;}}'

python ../../../genotools/mr.py TTC -I test.vcf.gz --fam test.fam > t.p.vcf
python ../../../genotools/mr.py Split -I test.vcf.gz --fam test.fam --dosage 
python ../../../genotools/mr.py GeneticScore -I FPG.vcf.gz -b FPG.beta.txt --fam test.fam


python ../../../genotools/mr.py GeneticScore -I FPG.vcf.gz -b FPG.beta.txt
python ../../../genotools/mr.py GeneticScore -I FPG.vcf.gz -b FPG.beta.txt --dosage

```
