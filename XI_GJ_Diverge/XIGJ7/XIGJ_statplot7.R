rm(list = ls())
library(tidyverse)
library(writexl)
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(export) #输出pdf中文字体
library(ggpubr)

args <- commandArgs(T)


# 合并所有样本binmap ------------------------------------------------------------


#samplst <- read.delim("sample.lst",header = F)
samplst <- read.delim(args[1],header = F)
sampName <- sapply(samplst[,1],as.character)
d <- NULL
XIGJ_counts <- NULL
XIGJ_prop <- NULL
for(i in 1:nrow(samplst)){
  # i =1
  d[[i]] <- read.delim(paste0(samplst[i,1],"/",samplst[i,1],".bin_NAF"),header = F)
  #d[[i]] <- read.delim(paste0(samplst[i,1],".bin_NAF"),header = F)
  counts <- table(d[[i]]$V4) #籼粳成分数目
  XIGJ_counts[[i]] <- as.data.frame(counts)
  props <- prop.table(counts)*100 #籼粳成分比例
  XIGJ_prop[[i]] <- as.data.frame(props)
}

data <- Reduce(function(x,y)inner_join(x,y,by=c("V1","V2","V3")),d)
names(data) <- c("Chr","Start","End",sampName)
data$Chr <- factor(data$Chr,levels = c(paste0("Chr",1:9),"Chr10","Chr11","Chr12"))
data2 <- data[order(data$Chr),]

mainDir <- getwd()
subDir <-"XIGJ_res"
if (file.exists(subDir)){
  setwd(file.path(mainDir, subDir))
} else {
  dir.create(file.path(mainDir, subDir))
  setwd(file.path(mainDir, subDir))
}
write.table(data2,"所有样本Binmap图谱.xls",row.names = F,col.names = T,quote = F,sep = "\t")



# binmap热图绘制 --------------------------------------------------------------------


heatd <- data2[-c(1:3)]
rownames(heatd) <- paste0("Bin",1:nrow(heatd))
# head(heatd)
# Heatmap(heatd,name = "Origin",show_row_names = F,border = 'black')

chr_split <- as.data.frame(table(data2$Chr))
chr_name <- sapply(chr_split$Var1,as.factor)
chr_count <- sapply(chr_split$Freq,as.numeric)

# pdf("XIGJ.binmap.allsamples.pdf")
# png("XIGJ.binmap.allsamples.png")
Heatmap(heatd,name = "Origin",
        show_row_names = F,
        show_column_names = T,
        border = 'black',
        row_split = rep(chr_name,times=chr_count),
        # column_split = colnames(heatd),
        col = c("XI_aus"="#3C5488", 
                "XI_ind"="#4DBBD5", 
                "XI_admix"="#8491B4", 
                "GJ_tmp"="#F39B7F", 
                "GJ_trp"="#FFFFB6", 
                "GJ_admix"="#E64B35", 
                "Admixed"="#D9D9D9"))

# dev.off()
graph2pdf(file="所有样本籼粳Binmap图谱.pdf",height=12,width=8)
graph2png(file="所有样本籼粳Binmap图谱.png",height=12,width=8)



#  籼粳成分比例统计 ---------------------------------------------------------------

XIGJ_counts_data <- Reduce(function(x,y)inner_join(x,y,by=c("Var1")),XIGJ_counts)
XIGJ_prop_data <- Reduce(function(x,y)inner_join(x,y,by=c("Var1")),XIGJ_prop)
names(XIGJ_counts_data) <- c("Origin",sampName)
names(XIGJ_prop_data) <- c("Origin",sampName)
list_xlsx <- list("籼粳Bin比例统计"=XIGJ_prop_data, "籼粳Bin数目统计"=XIGJ_counts_data)
write_xlsx(list_xlsx,"籼粳成分数目及比例统计.xlsx")


rownames(XIGJ_prop_data) <- XIGJ_prop_data$Origin
# pdf("所有样本籼粳成分比例.pdf")
# Heatmap(XIGJ_prop_data[-1],cluster_rows = F,cluster_columns = F,name = "Prop(%)") #,name = "Log2(Prop)"
Heatmap(XIGJ_prop_data[-1],cluster_rows = F,cluster_columns = F,name = "Prop(%)",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", XIGJ_prop_data[-1][i, j]), x, y, gp = gpar(fontsize = 10))})

# dev.off()
graph2pdf(file="所有样本籼粳成分比例.pdf",height=6,width=6)
graph2png(file="所有样本籼粳成分比例.png",height=6,width=6)



# 籼粳成分比例单样本饼图绘制 --------------------------------------------------------------
#https://www.jianshu.com/p/6bc99acbbb57

sampName
for(i in 1:length(sampName)){
  # i=5
  pd <- inner_join(XIGJ_counts_data[c(1,(i+1))],XIGJ_prop_data[c(1,(i+1))],by="Origin")
  names(pd)[2:3] <- c("Freq","Prop")
  
  ##可能会没有出现aus比例或极低的情况。其他成分暂且认为都会出现
  if( "XI_aus" %in% pd$Origin  &  length(pd$Origin)==7){
    pd$Prop<-paste0(round(pd$Prop,2),"%")
    pd$Origin<-paste0(pd$Origin,"(",pd$Prop,")")
    p2 <- ggpie(pd,"Freq",label=pd$Prop,fill="Origin",
                color="white",lab.font = c(3,"black"),lab.pos="out",
                palette= c("#D9D9D9","#E64B35","#F39B7F","#FFFFB6","#8491B4","#3C5488","#4DBBD5"))+
      theme(legend.position = "right",legend.text = element_text(size=8))
    ggsave(p2,filename=paste0(sampName[i],".籼粳比例.png"),width = 6,height = 5,dpi = 300)
    ggsave(p2,filename=paste0(sampName[i],".籼粳比例.pdf"),width = 6,height = 5)   
  }else if(length(pd$Origin)==6){
    pd$Prop<-paste0(round(pd$Prop,2),"%")
    pd$Origin<-paste0(pd$Origin,"(",pd$Prop,")")
    p2 <- ggpie(pd,"Freq",label=pd$Prop,fill="Origin",
                color="white",lab.font = c(3,"black"),lab.pos="out",
                palette= c("#D9D9D9","#E64B35","#F39B7F","#FFFFB6","#8491B4","#4DBBD5"))+
      theme(legend.position = "right",legend.text = element_text(size=8))
    ggsave(p2,filename=paste0(sampName[i],".籼粳比例.png"),width = 6,height = 5,dpi = 300)
    ggsave(p2,filename=paste0(sampName[i],".籼粳比例.pdf"),width = 6,height = 5)  
  }else{
    stop("出现非常纯的籼稻或粳稻品种！导致2个及以上的其他成分比例近乎于0")
  }

}

