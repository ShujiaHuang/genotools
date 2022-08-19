python ../../genometools/dosage.py -I convert2.vcf -N --vcf
python ../../genometools/dosage.py -I convert.vcf.gz -N --vcf
awk '$1!~/^##/{print $4"\t"$5"\t"$11}' convert2.vcf | awk -F ':' '{print $1}'

cat samtools_depth.txt|  python ../../scripts/merge_samtools_depth.py

