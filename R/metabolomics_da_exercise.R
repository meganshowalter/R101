## Writing our own functions
library(XLConnect)
install.packages("XLConnect")
library(ggplot2)
install.packages(plotly)
library(plotly)
library(mixOmics)
library(ReporteRs)
install.packages("ReporteRs")
library(svglite)
install.packages("svglite")
library(grid)
library(ggrepel)
library(gridExtra)

data.raw <- read.csv("./data/mx 107155 _lung cancer tissue_summer course_08-2015_submit.csv", header = F, stringsAsFactors = F)
data.raw <- data.raw[,-5] ## We get rid of the mass spec column. this is the one with mass spectra

head(data.raw)

data.raw[20,10] # inspecting the value of row number 20 and column number 10
data.raw[1:5,8:12] # inspecting the value of raw number 1 to 5 and column number 8 t 12

data.raw[,1] # first column – metabolite names
data.raw[1,] # first row – sample names

data.raw[10,] # class information

cpd.names <- data.raw[,1]
cpd.names[1:20]


#### Internal Standards
data.raw.is <- data.raw[grep("^z ",cpd.names),] # make a separate data frame for Internal standards. ^means look at begining of the word
head(data.raw.is)

###
data.raw <- data.raw[-grep("^z ",cpd.names),] # get rid of the internal standards. -sign means remmove

## Indexing information (from which row have metab or samples starting)
metstart.index <- (grep("BinBase",data.raw[,1])+1) # get the row from where we have metabolites are starting. do + bc intesnity starts after
row.ind <-  metstart.index:nrow(data.raw) # row numbers that has peak intensities. creates a vector
col.ind <- grep("bbdsa", data.raw[1,]) # column numbers that has peak intensities. this is GC instrument data file, change for LC

head(metstart.index)
head(row.ind)

### Header Details
sample.header <- data.raw[1:(metstart.index-1),(col.ind[1]-1)] #-2 is to grab above inchi kept as -1 for now but change later
metabolite.header <- data.raw [grep("BinBase",data.raw[,1]),1:(col.ind[1]-1)]

## Metabolites meta data
cpd.metadata <- data.raw[ metstart.index:nrow(data.raw)  , 1:(col.ind[1]-1) ]
#nrow means go to final row
#index for column until right b4 peak intenisty starts

## Subject MetabData
subject.metadata <- data.raw[ 1:(metstart.index-1), (col.ind[1]):ncol(data.raw)   ]

table(as.character(subject.metadata[10,])) ## Get the levels in the phenotype column.
#table takes the vector and looks for frequency of items in it

## Peak Intensities
numeric.data <- data.raw[row.ind,col.ind] # make a dataframe for only the metabolites.


#Principal Component Analysis
pcamat <- do.call("cbind",lapply(numeric.data,as.numeric))
#cbind to make columns together, lapply to make them all numeric to avoid issues with numeric if has them as character
#do.call is best way he found to do that
pcamat[which(is.na(pcamat)==TRUE)] <- min(pcamat[which(is.na(pcamat)==FALSE)])
#need to exclude

#to exclude 0s or make threshold for min feature value
pacmat[which(pcamat==0)] <- min(pcamat>=1)
pacmat[which(pcamat==0)] <- min(pcamat[-which(pcamat==0)])

library(mixOmics)
pca.res <- pca(pcamat, ncomp=2,max.iter=100)
pc.scores <- pca.res$loadings[[1]]
pca.res


# combine subject metadata and pc_scores in one data frame
subject.metadata.transpose <- as.data.frame(t(subject.metadata), stringsAsFactors = F)
## data frame need to be trasponsed so can be combined with PCA results.
#t is function for transposing
head(subject.metadata.transpose)

pcdf <- cbind(subject.metadata.transpose,pc.scores) #pcdf = PC data frame (df is also a list in R)
typeof(pcdf) #
is.matrix(pcdf) #
is.data.frame(pcdf) #
colnames(pcdf) <- c(sample.header,"PC1","PC2") #provide names for plot

p1 <- ggplot(pcdf, aes(PC1,PC2)) +
  scale_y_continuous(paste("PC2 - variance explained : ",signif(pca.res$explained_variance[2],digit=4),"(%) ",sep="")) +
  scale_x_continuous(paste("PC1 - variance explained : ",signif(pca.res$explained_variance[1],digit=4),"(%) ",sep="")) +
  theme_bw() +
  labs(title = "Principal component analysis (PCA)") +
  theme(
    plot.title = element_text(face="bold", size=14),
    axis.title.x = element_text(face="bold", size=12),
    axis.title.y = element_text(face="bold", size=12, angle=90),
    panel.grid.major = element_blank(), # switch off major gridlines
    panel.grid.minor = element_blank(), # switch off minor gridlines
    legend.position = c(0.1,0.1), # manually position the legend (numbers being from 0,0 at bottom left of whole plot to 1,1 at top right)
    legend.title = element_blank(), # switch off the legend title
    legend.text = element_text(size=12),
    legend.key.size = unit(1.5, "lines"),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.spacing = unit(.05, "cm"),
    axis.text.x = element_text(size=0,angle = 45, hjust = 1),
    axis.text.y = element_text(size=0,angle = 45, hjust = 1)
  )
p2 <- p1 + geom_point(aes(shape=phenotype, colour=phenotype ), size=3)
p3 <- p2 +  stat_ellipse( aes(colour=phenotype ), linetype = 1 )

plot(p2)
plot(p3)
ggsave("PCA.png", width = 12, height = 8, dpi = 84)
ggsave("PCA.pdf", width = 12, height = 8, dpi = 84)

# t.test

groupA.ind <- which(subject.metadata[10,]=="non-malignant") ## non-cancer tissues
groupB.ind <- which(subject.metadata[10,]=="tumor") ## tumor tissues
#location for where they are but not making new df or sorting

## we need to create a function for running ttest.

runttest <- function(index) {
  results <- t.test(pcamat[index,groupA.ind],pcamat[index,groupB.ind])
  foldchange <- median(pcamat[index,groupA.ind])/median(pcamat[index,groupB.ind])
  unlist(c(cpd.metadata[index,],results$p.value,foldchange))
}
#unlist take any df or list and convert to vector, and if single character will convert to character, if no character, all numeric
#likes to convert to character to keep from converting to factor which will mess it up
#make everything character, then convert pvalue and FC later (below) back to number
runttest(1) ## run the ttest for metabolite 1.

ttest.results <- do.call(rbind,lapply(1:nrow(pcamat), runttest)) ## run the ttest for all the compounds and combine that results into one dataframe.
colnames(ttest.results) <- c(metabolite.header,"PValue","FoldChange")
ttest.results <- as.data.frame(ttest.results, stringsAsFactors = F) #convert to df so can use for graphics
ttest.results$PValue <- as.numeric(ttest.results$PValue ) #convert character to numeric
ttest.results$FoldChange <- as.numeric(ttest.results$FoldChange )
pvalue.adjusted <- p.adjust(ttest.results$PValue,method = "fdr") ## FDR adjustment of PValues

ttest.results <- cbind(ttest.results,PValue.Adjusted=pvalue.adjusted)
write.table(ttest.results, file="ttestresults.txt", col.names = T, row.names = F,quote = F, sep = "\t") ## export the results to an external file.

# Volcano plot
volcanodf <- data.frame(name=ttest.results$`BinBase name`, pvalue=ttest.results$PValue, foldchange = ttest.results$FoldChange)
volcanodf$Significant <- ifelse(volcanodf$pvalue < 0.005, "PValue < 0.005", "Not Sig")
table(volcanodf$Significant)

pv1 <- ggplot(volcanodf, aes(x = log10(foldchange), y = -log10(pvalue))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 16) +
  geom_text_repel(
    data = subset(volcanodf, pvalue < 0.005),
    aes(label = name),
    size = 2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )

plot(pv1)

ggsave("volcano-1.png", width = 12, height = 8, dpi = 84)

# box-whisker plots

data_bw <- cbind(as.data.frame(t(pcamat)),phenotype =pcdf$phenotype)
#adding extra column phenotype from pcdf

k = 4 # valine

f1 = ggplot(data_bw, aes(phenotype, data_bw[[ colnames(data_bw) [ k ]]]) ) +
  geom_boxplot(fill = "grey", colour = "#3366FF",outlier.shape = 3) +
  geom_point(aes(color=phenotype), size=2, position = position_jitterdodge()) +
  scale_y_continuous("Metabolite levels (normalized ion count)") +
  scale_shape_manual(values=c(15,19,17,5)) +
  scale_fill_manual(values=c("#CC0000", "#006600", "#669999", "#00CCCC","#660099", "#CC0066", "#FF9999", "#FF9900","black", "black", "black", "black", "black")) +
  theme_bw() + # make the theme black-and-white rather than grey (do this before font changes, or it overrides them)
  labs(title = paste(" Levels of  ", cpd.metadata$V1[ k  ],sep=" : ")) +
  theme(
    plot.title = element_text(face="bold", size=14), # use theme_get() to see available options
    axis.title.x = element_text(face="bold", size=12),
    axis.title.y = element_text(face="bold", size=12, angle=90),
    panel.grid.major = element_blank(), # switch off major gridlines
    panel.grid.minor = element_blank(), # switch off minor gridlines
    legend.position = c(0.1,0.9), # manually position the legend (numbers being from 0,0 at bottom left of whole plot to 1,1 at top right)
    legend.title = element_blank(), # switch off the legend title
    legend.text = element_text(size=12),
    legend.key.size = unit(1.5, "lines"),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.margin = unit(.05, "cm"),
    axis.text.x = element_text(size=12,angle = 45, hjust = 1)
  )
plot(f1)

## let's make a loop to go over all the significant metabolites.

topcpds.index <- order(ttest.results$PValue) [1:10]

for (k in topcpds.index) {
  f1 = ggplot(data_bw, aes(phenotype, data_bw[[ colnames(data_bw) [ k ]]]) ) +
    geom_boxplot(fill = "grey", colour = "#3366FF",outlier.shape = 3) +
    geom_point(aes(color=phenotype), size=2, position = position_jitterdodge()) +
    scale_y_continuous("Metabolite levels (normalized ion count)") +
    scale_shape_manual(values=c(15,19,17,5)) +
    scale_fill_manual(values=c("#CC0000", "#006600", "#669999", "#00CCCC","#660099", "#CC0066", "#FF9999", "#FF9900","black", "black", "black", "black", "black")) +
    theme_bw() + # make the theme black-and-white rather than grey (do this before font changes, or it overrides them)
    labs(title = paste(" Levels of  ", cpd.metadata$V1[ k  ],sep=" : ")) +
    theme(
      plot.title = element_text(face="bold", size=14), # use theme_get() to see available options
      axis.title.x = element_text(face="bold", size=12),
      axis.title.y = element_text(face="bold", size=12, angle=90),
      panel.grid.major = element_blank(), # switch off major gridlines
      panel.grid.minor = element_blank(), # switch off minor gridlines
      legend.position = c(0.1,0.9), # manually position the legend (numbers being from 0,0 at bottom left of whole plot to 1,1 at top right)
      legend.title = element_blank(), # switch off the legend title
      legend.text = element_text(size=12),
      legend.key.size = unit(1.5, "lines"),
      legend.key = element_blank(), # switch off the rectangle around symbols in the legend
      legend.margin = unit(.05, "cm"),
      axis.text.x = element_text(size=12,angle = 45, hjust = 1)
    )
  plotname <- paste0("bw_plots_",cmpd.metadataV1[k],".jpg")
  ggsave(plotname, width = 12, height = 8, dpi = 84)
}


## add multiple graphs on one page of pdf.

plot.box.whisker <- function(k) {
  f1 = ggplot(data_bw, aes(phenotype, data_bw[[ colnames(data_bw) [ k ]]]) ) +
    geom_boxplot(fill = "grey", colour = "#3366FF",outlier.shape = 3) +
    geom_point(aes(color=phenotype), size=2, position = position_jitterdodge()) +
    scale_y_continuous("Metabolite levels (normalized ion count)") +
    scale_shape_manual(values=c(15,19,17,5)) +
    scale_fill_manual(values=c("#CC0000", "#006600", "#669999", "#00CCCC","#660099", "#CC0066", "#FF9999", "#FF9900","black", "black", "black", "black", "black")) +
    theme_bw() + # make the theme black-and-white rather than grey (do this before font changes, or it overrides them)
    labs(title = paste(" Levels of  ", cpd.metadata$V1[ k  ],sep=" : ")) +
    theme(
      plot.title = element_text(face="bold", size=14), # use theme_get() to see available options
      axis.title.x = element_text(face="bold", size=12),
      axis.title.y = element_text(face="bold", size=12, angle=90),
      panel.grid.major = element_blank(), # switch off major gridlines
      panel.grid.minor = element_blank(), # switch off minor gridlines
      legend.position = c(0.1,0.9), # manually position the legend (numbers being from 0,0 at bottom left of whole plot to 1,1 at top right)
      legend.title = element_blank(), # switch off the legend title
      legend.text = element_text(size=12),
      legend.key.size = unit(1.5, "lines"),
      legend.key = element_blank(), # switch off the rectangle around symbols in the legend
      legend.margin = unit(.05, "cm"),
      axis.text.x = element_text(size=12,angle = 45, hjust = 1)
    )
  f1
}

plotlist <- lapply(topcpds.index,plot.box.whisker) # we collect all the plots into a list object
ml <- marrangeGrob(plotlist, nrow=2,ncol=3) # then we arrange the list
ggsave("multipage2.pdf", ml, height = 10, width = 10)




















# write an if statement to get the names of all the internal standards
for (i in cpd.names) {
  if( length(grep('^z ',i)>0) ){
    print(i)
  }
}


# getting the row number and name of int standards.

for (i in 1:length(cpd.names)) {
  if( length(grep('^z ',cpd.names[i])>0) ){
    print(c(i,cpd.names[i]))
  }
}


# lets get all the samples ids that has a label “tumor”
data.raw[1, grep("tumor",data.raw[10,] )  ]



# lets get all the columns that has sample data,
data.raw[1,grep("bbdsa", data.raw[1,])]




