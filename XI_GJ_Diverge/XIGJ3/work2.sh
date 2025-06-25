#conda activate r4
Rscript XIGJ_statplot3.R sample.lst
##clone genes
sed '1d' XIGJ_res/所有样本Binmap图谱.xls >binmap.bed
bedtools intersect -a binmap.bed -b cloneGenes.bed -wa -wb >XIGJ_res/binmap.cloneGene.xls

