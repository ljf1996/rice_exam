library(ggplot2)
library(scales)
library(stringr)

args <- commandArgs(TRUE)

chrlist <- args[1]
input <- args[2]
output <- args[3]

data <- read.delim(chrlist,header=F)
data$V2 <- (data$V2+1000)/1000000  #换成Mb
bin <- read.delim(input,header=F)
bin[2:3] <- bin[2:3]/1000000  #换成Mb

plot.height <- 6
plot.width <- length(data$V1) * 0.5 + 4

#pdf(output, width = plot.width, height = plot.height)

data$V1 <- factor(data$V1, levels = data$V1)
bin$V1 <- factor(bin$V1, levels = data$V1)

types <- length(levels(bin$V4))

p <-ggplot(data=data) +
  geom_rect(aes(xmin = as.numeric(V1) - 0.2, 
				xmax = as.numeric(V1) + 0.2 , 
				ymax = V2, ymin = 0),
				colour="black", fill = "white") +
  #coord_flip() +
  geom_rect(data=bin, aes(xmin = as.numeric(V1) - 0.18, 
				xmax = as.numeric(V1) + 0.18, 
				ymax = V2, ymin = V3, fill=V4)) +
  #scale_fill_manual(values = c("aus"="#FFFF00", "ind"="#556B2F", "INDICA"="#9ACD32", "tej"="#00BFFF", "trj"="#8B008B", "JAPONICA"="#0000FF")) +
  #scale_fill_manual(values = c("sp_aus"="#FF0000", "sp_ind"="#FFD700", "ss_INDx"="#FF8C00", "sp_tej"="#00FFFF", "sp_trj"="#800080", "ss_JAPx"="#4080C0", "admixed"="#C0C0C0")) +
  #scale_fill_manual(values = c("sp_aus"="#3C5488", "sp_ind"="#4DBBD5", "ss_INDx"="#8491B4", "sp_tej"="#F39B7F", "sp_trj"="#FFFFB6", "ss_JAPx"="#E64B35", "admixed"="#D9D9D9")) +
  scale_fill_manual(values = c("XI"="#4DBBD5", "GJ"="#E64B35", "Admixed"="#D9D9D9")) +
  guides(fill=guide_legend(title="Origin")) +
  theme(axis.text.x = element_text(colour = "black"),
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank()) + 
  scale_x_discrete(position = "top", name = "Chromosomes", limits = data$V1) +
  scale_y_continuous(trans="reverse", labels = comma) +
  ylab("Region (Mb)")
ggsave(plot=p,filename=paste0(output,".binmap.pdf"),width = plot.width, height = plot.height)
ggsave(plot=p,filename=paste0(output,".binmap.png"),width = plot.width, height = plot.height,dpi=300)
#dev.off()
