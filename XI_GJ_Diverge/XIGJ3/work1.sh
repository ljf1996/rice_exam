## get all sample list
#bcftools query -l ../vcf/filtered_variant.vcf.gz >sample.lst

## each sample introgression
#sample=$1
cat sample.lst | while read sample;do
#cat sample.lst |sed -n '2p' | while read sample;do
echo $sample
mkdir $sample
## extract sample vcf
vcftools --gzvcf ../vcf/filtered_variant.vcf.gz --recode --recode-INFO-all --stdout --indv $sample >$sample/${sample}.vcf
## build haplotype
mkdir -p $sample/haplotype
perl /mnt/project/Zayou_center/Denovo_project/project/09.indica_japonica_mix/3KRG-HAP-master/scripts/gatk_vcf_to_haplotype_with_varlist.pl --vcf $sample/${sample}.vcf --var /mnt/project/Zayou_center/Denovo_project/project/09.indica_japonica_mix/3KRG-HAP-master/data/3K-SNP.varlist --out $sample/haplotype --nohet
cat $sample/haplotype/*.haplotype >$sample/${sample}.haplotype
## map NAF score
mkdir -p $sample/NAF_score
perl /mnt/project/Zayou_center/Denovo_project/project/09.indica_japonica_mix/3KRG-HAP-master/scripts/classify_sample_haplotype_score.pl $sample/${sample}.haplotype /mnt/project/Zayou_center/Denovo_project/project/09.indica_japonica_mix/3KRG-HAP-master/data/3K-HAP.haplotype.NAF_score $sample/NAF_score
perl /mnt/project/Zayou_center/Denovo_project/project/09.indica_japonica_mix/3KRG-HAP-master/scripts/scan_haplotype_stdratio.pl $sample/NAF_score/${sample}.hapratio /mnt/project/Zayou_center/Denovo_project/project/09.indica_japonica_mix/3KRG-HAP-master/data/window.10k.bin.bed $sample/sample.bin_NAF
## assign subpop
#perl /mnt/project/Zayou_center/Denovo_project/project/09.indica_japonica_mix/3KRG-HAP-master/scripts/dissect_rice_bin.pl $sample/sample.bin_NAF >$sample/${sample}.bin_NAF
perl dissect_rice_bin3.pl $sample/sample.bin_NAF >$sample/${sample}.bin_NAF
## plot
#Rscript /mnt/project/Zayou_center/Denovo_project/project/09.indica_japonica_mix/3KRG-HAP-master/scripts/draw_bin.rice.R /mnt/project/Zayou_center/Denovo_project/project/09.indica_japonica_mix/3KRG-HAP-master/data/chr.len $sample/${sample}.bin_NAF $sample/${sample}.bin_NAF.pdf
Rscript draw_bin.rice3.R /mnt/project/Zayou_center/Denovo_project/project/09.indica_japonica_mix/3KRG-HAP-master/data/chr.len $sample/${sample}.bin_NAF $sample/$sample

done

#conda activate r4
#script XIGJ_statplot.R sample.lst
