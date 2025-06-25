rm(list = ls())
library(tidyverse)
args <- commandArgs(T)

# 合并所有样本binmap ------------------------------------------------------------


# samplst <- read.delim("sample.lst",header = F)
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
write.table(data2,"merge_binmap.txt",row.names = F,col.names = T,quote = F,sep = "\t")


# binmap热图绘制 --------------------------------------------------------------------


# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

heatd <- data2[-c(1:3)]
rownames(heatd) <- paste0("Bin",1:nrow(heatd))
# head(heatd)
# Heatmap(heatd,name = "Origin",show_row_names = F,border = 'black')

chr_split <- as.data.frame(table(data2$Chr))
chr_name <- sapply(chr_split$Var1,as.factor)
chr_count <- sapply(chr_split$Freq,as.numeric)

pdf("XIGJ.binmap.allsamples.pdf")
# png("XIGJ.binmap.allsamples.png")
Heatmap(heatd,name = "Origin",
        show_row_names = F,
        show_column_names = F,
        border = 'black',
        row_split = rep(chr_name,times=chr_count),
        column_split = colnames(heatd),
        col = c("sp_aus"="#3C5488", 
                "sp_ind"="#4DBBD5", 
                "ss_INDx"="#8491B4", 
                "sp_tej"="#F39B7F", 
                "sp_trj"="#FFFFB6", 
                "ss_JAPx"="#E64B35", 
                "admixed"="#D9D9D9"))

dev.off()



#  籼粳成分比例统计 ---------------------------------------------------------------


# XIGJ_counts[[1]]
# XIGJ_prop[[1]]
XIGJ_counts_data <- Reduce(function(x,y)inner_join(x,y,by=c("Var1")),XIGJ_counts)
XIGJ_prop_data <- Reduce(function(x,y)inner_join(x,y,by=c("Var1")),XIGJ_prop)
names(XIGJ_counts_data) <- c("Origin",sampName)
names(XIGJ_prop_data) <- c("Origin",sampName)

write.csv(XIGJ_counts_data,"XIGJ_counts_data.csv",row.names = F)
write.csv(XIGJ_prop_data,"XIGJ_proportion_data.csv",row.names = F)

rownames(XIGJ_prop_data) <- XIGJ_prop_data$Origin
pdf("XIGJ.proportion.allsamples.pdf")
Heatmap(XIGJ_prop_data[-1],cluster_rows = F,cluster_columns = F,name = "Prop(%)")
# Heatmap(log2(XIGJ_prop_data[-1]),cluster_rows = F,cluster_columns = F,name = "Log2(Prop)")
dev.off()


# 籼粳成分比例单样本饼图绘制 --------------------------------------------------------------


library(ggpubr)
#https://www.jianshu.com/p/6bc99acbbb57

sampName
for(i in 1:length(sampName)){
  # i=1
  pd <- inner_join(XIGJ_counts_data[c(1,(i+1))],XIGJ_prop_data[c(1,(i+1))],by="Origin")
  names(pd)[2:3] <- c("Freq","Prop")
  pd$Prop<-paste0(round(pd$Prop,2),"%")
  pd$Origin<-paste0(pd$Origin,"(",pd$Prop,")")
  
  p2 <- ggpie(pd,"Freq",label=pd$Prop,fill="Origin",
              color="white",lab.font = c(3,"black"),lab.pos="out",
              palette=c("#D9D9D9","#3C5488","#4DBBD5","#F39B7F","#FFFFB6","#8491B4","#E64B35"))+
    theme(legend.position = "right",legend.text = element_text(size=8))
  # p2
  ggsave(p2,filename=paste0(sampName[i],".proportion.png"),width = 6,height = 5,dpi = 300)
  ggsave(p2,filename=paste0(sampName[i],".proportion.pdf"),width = 6,height = 5)
}

