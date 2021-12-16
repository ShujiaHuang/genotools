```bash
python ../../../genotools/mr.py TTC -I test2.vcf.gz --fam test.fam| awk '$1!~/^##/ {print $1 "\t" $2 "\t" $10 "\t" $11 "\t" $12}' | awk 'BEGIN{FS="[\t:]"; OFS="\t";}{if(NR==1){print $0}else{print $1,$2,$3,$6,$9;}}'

python ../../../genotools/mr.py PRS -I test.vcf.gz --betavaluefile test.gwas.beta.tsv --famfile test.fam > tt
python ../../../genotools/mr.py TTC -I test2.vcf.gz --fam test.fam > t.p.vcf
```
