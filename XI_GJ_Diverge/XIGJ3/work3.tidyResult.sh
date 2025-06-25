mkdir -p 划分三类亚群结果/{单样本结果,群体比较}
cat sample.lst |while read id;
do
    mkdir -p 划分三类亚群结果/单样本结果/$id
    cp $id/${id}.binmap*  划分三类亚群结果/单样本结果/$id
    cp $id/${id}.bin_NAF  划分三类亚群结果/单样本结果/$id/${id}.binmap.txt
    cp $id/${id}.haplotype 划分三类亚群结果/单样本结果/$id/${id}.haplotype.txt
    cp XIGJ_res/${id}.籼粳比例* 划分三类亚群结果/单样本结果/$id
    cp XIGJ_res/binmap.cloneGene.xls XIGJ_res/籼粳成分数目及比例统计.xlsx 划分三类亚群结果/群体比较
    cp XIGJ_res/所有* 划分三类亚群结果/群体比较

done
