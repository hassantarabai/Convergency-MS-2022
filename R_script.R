### Convergency R scripts for the publication Microbiomes of blood feeding triatomines in the context of their predatory relatives and the environment ###

#load libraries
library(phyloseq) #microbiome analysis #version 1.38.0
library(decontam) #For identifying contaminant OTUs
library(MicEco) #Venn based analysis
library(readxl) #import excel sheet files
library(tibble) #converting columns to row names
library(dplyr) # data handling
library(microeco) #microbiome analysis tool #version 0.6.5
library(ggplot2) #plotting #version 3.3.5
library(tidytree) #tidy taxonomy #version 0.3.6
library(ape) #phylogenetic tree processing #version 5.5
library(agricolae) #duncan test in alpha diversity #version 1.3-5
library(igraph) #network analysis #version 1.2.11
library(file2meco) #import objects created in different packages
library(ggalluvial)
library(ggh4x)
library(FSA) #for KW Dunn test
library(ape) #phylogenetic tree processing #version 5.5
library(agricolae) #duncan test in alpha diversity #version 1.3-5
library(igraph) #network analysis #version 1.2.11
library(file2meco) #import objects created in different packages
library(viridis) #colour package
library(ggsci) #colour package
library(svglite) # produce scalable vector package

##Analysis of positive controls

#MCE controls
mce=read.csv("poscontrols_mce.csv",header=T) # load the table of MCE controls in the R environment. The table contains the taxa composition and the sequenced and expected frequencies for each taxon
mce_c<-rep(mce$Composition) # create a vector for the expected taxa composition (names of taxa)
mce_s<-c(rep("MCE",10),rep("MCE1",10),rep("Standard",10)) # create a vector of the samples being analysed
mce_v=c(rep(mce$MCE_per),rep(mce$MCE1_per),rep(mce$Standard)) # create a vector of the frequencies of the samples being analysed

# mce_s and mce_v need to be in the same order
data <- data.frame(mce_s,mce_c,mce_v) # join the vectors in a dataframe
MCEplot <- ggplot(data, aes(fill=mce_c, y=mce_v, x=mce_s)) + geom_bar(position="stack", stat="identity") + scale_fill_brewer(palette = "Spectral") # create the plot
ggsave("MCE_plot.svg",MCEplot) # save the plot in SVG format for graphical adjustment in InkScape

#MCS controls
mcs=read.csv("poscontrols_mcs.csv",header=T)
mcs_c<-rep(mcs$Composition)
mcs_s<-c(rep("MCS",10),rep("MCS1",10),rep("Standard",10))
mcs_v=c(rep(mcs$MCS_per),rep(mcs$MCS1_per),rep(mcs$Standard))
data2 <- data.frame(mcs_s,mcs_c,mcs_v)
data2
MCSplot <- ggplot(data2, aes(fill=mcs_c, y=mcs_v, x=mcs_s)) + geom_bar(position="stack", stat="identity") + scale_fill_brewer(palette = "Spectral")
ggsave("MCS_plot.svg",MCSplot)

#PC1, PC3, ZymoStd
zym=read.delim("PC1_PC3_ZymoStd.tab")
zym_c<-rep(zym$Composition)
zym_s<-c(rep("PC1",8), rep("PC3",8),rep("ZymoStandard",8))
zym_v=c(rep(zym$PC1),rep(zym$PC3),rep(zym$ZymoStandard))
datazym <- data.frame(zym_s,zym_c,zym_v)
zymplot <- ggplot(datazym, aes(fill=zym_c, y=zym_v, x=zym_s)) + geom_bar(position="stack", stat="identity") + scale_fill_brewer(palette = "Spectral")
ggsave("PC1_PC3_ZymoStd_plot.svg",zymplot)

#PC2, PC31, ZymoLog
zymlog=read.delim("PC2_PC31_ZymoLog.tab")
zymlog_c<-rep(zymlog$Composition)
zymlog_s<-c(rep("PC2",8),rep("PC31",8),rep("ZymoLog",8))
zymlog_v=c(rep(zymlog$PC2),rep(zymlog$PC31),rep(zymlog$ZymoLog))
datazymlog <- data.frame(zymlog_s,zymlog_c,zymlog_v)
zymlogplot <- ggplot(datazymlog, aes(fill=zymlog_c, y=zymlog_v, x=zymlog_s)) + geom_bar(position="stack", stat="identity") + scale_fill_brewer(palette = "Spectral")
ggsave("PC2_PC31_ZymoLog.svg",zymlogplot)

# all plots were modified for graphical improvements and joined in a single figure in InkScape

## Identifying contaminants in the obtained OTU table with taxonomy

#set working directory
setwd("C:/Users/hassa/Desktop/Convergency/Decontam/")

#Loading datasheets
otu_original<- read_excel("C:/Users/hassa/Desktop/Convergency/Decontam/otu_original.xlsx") #upload the OTU abundance file
tax_original<- read_excel("C:/Users/hassa/Desktop/Convergency/Decontam/tax_original.xlsx") #upload the taxonomy file
sample_df<- read_excel("C:/Users/hassa/Desktop/Convergency/Decontam/sample_df.xlsx") #Upload the metadata


# Generate a phyloseq object
#1-defining row names from the otu column
otu_original <- otu_original %>%
  tibble::column_to_rownames("OTU") 
#2-identifying metadata and taxonomy files	
samples_df <- sample_df %>% 
  tibble::column_to_rownames("SampleID") 
tax_original <- tax_original %>% 
  tibble::column_to_rownames("OTU")
#3-transforming otu_mat and tax_mat into matrixes
otu_original <- as.matrix(otu_original)
tax_original <- as.matrix(tax_original)

# Transform into phyloseq objects
OTU = otu_table(otu_original, taxa_are_rows = TRUE)
TAX = tax_table(tax_original)
samples = sample_data(samples_df)
alldata <- phyloseq(OTU, TAX, samples)

#Subset the taxa to include only OTUs from bacterial kingdom and exclude taxonomical assignments with family "mitochondria" and class "Chloroplast".
clean_data1 <- alldata %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Kingdom != "Archaea" &
      Kingdom != "Animalia" &
      Kingdom != "Eukaryota" &
      Kingdom != "ND" &
      Kingdom != "Plantae" &
      Family != "mitochondria" &
      Class != "Chloroplast"
  )

#Set seed 100
set.seed(100)


#identify Contaminants by Frequency
contamdf.freq <- isContaminant(clean_data1, method="frequency", conc="PCRp_conc")
head(contamdf.freq)

table(contamdf.freq$contaminant)

head(which(contamdf.freq$contaminant))

#identify negative control data in clean_data1
sample_data(clean_data1)$is.neg <- sample_data(clean_data1)$Sample_type == "negative"

#Identify contaminants by prevalence
contamdf.prev<- isContaminant(clean_data1, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant) #Show prevalence result table

#OTU table with its associated taxonomy were filtered manually for identified contaminants by decontam package
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###Microbiome analysis###

#Loading datasheets
otu_mat<- read_excel("C:/Users/hassa/Desktop/Convergency/Decontam/otu_mat.xlsx") #upload the OTU abundance file
tax_mat<- read_excel("C:/Users/hassa/Desktop/Convergency/Decontam/tax_mat.xlsx") #upload the taxonomy file
sample_df<- read_excel("C:/Users/hassa/Desktop/Convergency/Decontam/sample_df.xlsx") #Upload the metadata
OTU_tree <- read.tree("OTU_tree.nwk") #Upload the generated newick tree of assigned OTUs

# Generate a phyloseq object
#1-defining row names from the otu column
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("OTU") 
#2-identifying metadata and taxonomy files
samples_df <- sample_df %>% 
  tibble::column_to_rownames("SampleID") 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("OTU")
#3-transforming otu_mat and tax_mat into matrixes
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

# Transform into phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
tree = read_tree(OTU_tree)
alldata1 <- phyloseq(OTU, TAX, samples, tree)

#Subset the taxa to include only OTUs from bacterial kingdom annd exclude taxonomical assignments with family "mitochondria" and class "Chloroplast".
clean_data <- alldata1 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Kingdom != "Archaea" &
      Kingdom != "Animalia" &
      Kingdom != "Eukaryota" &
      Kingdom != "ND" &
      Kingdom != "Plantae" &
      Family != "mitochondria" &
      Class != "Chloroplast"
  )

##Rarefaction of data



#Rarefaction at 800 sample size for analysis of nest material
set.seed(5) #Set the seed for data reproducibility
env_800 <- rarefy_even_depth(clean_data, sample.size = 800)

#Rarfaction at 1000 sample size for analysis of sample types (Ticks, assassin bugs and Triatoma)
set.seed(5) #Set seed for data reproducibility
clean_data_1000 <- rarefy_even_depth(clean_data, sample.size = 1000)

##Analysis of microbiome composition of nest material

#Setting working directory for nest material analysis
setwd("C:/Users/hassa/Desktop/Convergency/Nest material/")

#Converting phyloseq object to meco dataset
meco_nest <- phyloseq2meco(env_800)

#make the OTU and sample information consistent across all files in the dataset object
meco_nest$tidy_dataset()

#clone the full data set
group1 <- clone(meco_nest)

#subsetting group1 to include only samples from nest material
group1$sample_table <- subset(group1$sample_table, Sample_type == "Nest material")

#make the OTU and sample information consistent across all files in the data set object
group1$tidy_dataset()

#calculate abundance for taxa
group1$cal_abund()

#calculate alpha-diversity for samples
group1$cal_alphadiv(measures = NULL, PD = FALSE)

#calculate beta-diversity for samples
group1$cal_betadiv(unifrac = FALSE)

#create objects for abundance analysis
t1 <- trans_abund$new(dataset = group1, taxrank = "Genus", ntaxa = 30)
t2 <- trans_abund$new(dataset = group1, taxrank = "Phylum", ntaxa = 20)

#Generate heat map for the relative abundance of top 30 genera for nest material with their distribution on nests and locations
t1$plot_heatmap(
    color_values = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")),
    facet = "Location",
    order_facet = c("CampBullis","LacklandAFB","Gainesville"),
    facet2 = "Nest",
    x_axis_name = NULL,
    order_x = NULL,
    withmargin = TRUE,
    plot_numbers = FALSE,
    plot_text_size = 4,
    plot_breaks = NULL,
    margincolor = "white",
    plot_colorscale = "log10",
    min_abundance = 0.01,
    max_abundance = NULL,
    strip_text = 9,
    xtext_size = 8,
    ytext_size = 11,
    xtext_keep = TRUE,
    xtitle_keep = TRUE,
    grid_clean = TRUE,
    xtext_type_hor = TRUE,
    pheatmap = FALSE,
)

#Generate heat map for the relative abundance of top 20 phylum for nest material with their distribution on nests and locations

t2$plot_heatmap(
    color_values = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")),
    facet = "Location",
    facet2 = "Nest",
    order_facet = c("CampBullis","LacklandAFB","Gainesville"),
    x_axis_name = NULL,
    order_x = NULL,
    withmargin = TRUE,
    plot_numbers = FALSE,
    plot_text_size = 4,
    plot_breaks = NULL,
    margincolor = "white",
    plot_colorscale = "log10",
    min_abundance = 0.01,
    max_abundance = NULL,
    strip_text = 9,
    xtext_size = 8,
    ytext_size = 11,
    xtext_keep = TRUE,
    xtitle_keep = TRUE,
    grid_clean = TRUE,
    xtext_type_hor = TRUE,
    pheatmap = FALSE,
)

#Statistical analysis of alpha diversity of nest material based on pair comparison of locations using Dunn's Kruskal-Wallis Multiple Comparisons
t3 <- trans_alpha$new(dataset = group1, group = "Location") #creating a trans_alpha object
t3$cal_diff(method = "KW_dunn") #Do the pair compirson analysis

#View result
View(t3[["res_diff"]])

#Plot the prepared alpha diversity based on Dunn's Kruskal-Wallis Multiple Comparisons method and Chao1 measure
t3$res_diff %<>% base::subset(Significance != "ns") #remove "non significant" sign from showing in the plot
t3$plot_alpha(measure = "Shannon", xtext_size = 15) #Plot the alpha diversity using Shannon measure
t3$plot_alpha(measure = "Chao1", xtext_size = 15) #Plot alpha diversity using Chao1 measure
t3$plot_alpha(measure = "Observed", xtext_size = 15) #Plot alpha diversity using Observed measure

#Beta-diversity analysis using Non-Metric Dimensional Scaling (NMDS) based on Bray Curtis distances for assessing microbial composition of nest material between locations
t4 <- trans_beta$new(dataset = group1, group = "Location", measure = "bray") #creating a trans_beta object
t4$cal_ordination(ordination = "NMDS") #Using NMDS for analysis
t4$plot_ordination(plot_color = "Location", plot_shape = "Location", plot_type = c("point", "chull", "centroid"), add_sample_label = "Nest")
t4$cal_manova(manova_all = TRUE) #Manova analysis for nest material samples based on group "Location"


##Identification of OTUs associated with nest material (Minimum threshold set for presence of OTUs in minimal 30% of nest material alone or in nest material 
#in combination with other sample types.

#subsetting env_800 phyloseq to filter out negative and positive controls
sample_types_800 <- subset_samples(env_800,Sample_type %in% c("Tick", "Triatoma", "assassin", "Nest material"))

#Identifying OTUs present in 0.3, 0.5 and 0.8 fraction of all sample types
#1-venn at 0.3 fraction
venn<-ps_venn(sample_types_800,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.3)

#2-venn at 0.5 fraction 
venn<-ps_venn(sample_types_800,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.5)

#3- venn at 0.8 fraction
venn<-ps_venn(sample_types_800,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8)

##Triatoma analysis

#Setting working directory for triatoma analysis
setwd("C:/Users/hassa/Desktop/Convergency/triatoma/")

#Subseting the dataset to include only Triatoma samples
triatoma <- clean_data_1000 %>% subset_samples(Sample_type == "Triatoma")

#Analysis of common OTUs shared between different species of Triatoma at venn fraction 1, 0.9 and 0.5

#1-Venn at different samples' population fractions
venn1<-ps_venn(triatoma,"Species", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 1)
venn2<-ps_venn(triatoma,"Species", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.9)
venn3<-ps_venn(triatoma,"Species", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8)
venn4<-ps_venn(triatoma,"Species", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.5)

#Abundance analysis for Triatoma species
meco_triat <- phyloseq2meco(triatoma) #converting phyloseq object "triatoma" to meco dataset

#make the OTU and sample information consistent across all files in the dataset object
meco_triat$tidy_dataset()

#clone the full data set
triat1 <- clone(meco_triat)
triat2 <-clone (meco_triat) #For triatoma instar analysis

#Subsetting "triat2" dataset to include only triatoma instars with a defined variable
triat2$sample_table <- subset(triat2$sample_table, Instar != "ND")

#make the OTU and sample information consistent across all files in the dataset "triat2"
triat2$tidy_dataset()

#calculate abundance for taxa
triat1$cal_abund()
triat2$cal_abund()


#create objects for abundance analysis
g1 <- trans_abund$new(dataset = triat1, taxrank = "Genus", ntaxa = 30)
g2 <- trans_abund$new(dataset = triat1, taxrank = "Phylum", ntaxa = 20)
g3 <- trans_abund$new(dataset = triat2, taxrank = "Genus", ntaxa = 20)
g4 <- trans_abund$new(dataset = triat2, taxrank = "Phylum", ntaxa = 20)

Generate heat map for the relative abundance of top 30 genera in Triatoma species
g1$plot_heatmap(
    color_values = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")),
    facet = "Species",
    x_axis_name = NULL,
    order_x = NULL,
    withmargin = TRUE,
    plot_numbers = FALSE,
    plot_text_size = 4,
    plot_breaks = NULL,
    margincolor = "white",
    plot_colorscale = "log10",
    min_abundance = 0.01,
    max_abundance = NULL,
    strip_text = 9,
    xtext_size = 7,
    ytext_size = 11,
    xtext_keep = TRUE,
    xtitle_keep = TRUE,
    grid_clean = TRUE,
    xtext_type_hor = TRUE,
    pheatmap = FALSE,
)


#Generate heat map for the relative abundance of top 20 phyla in Triatoma species
g2$plot_heatmap(
    color_values = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")),
    facet = "Species",
    x_axis_name = NULL,
    order_x = NULL,
    withmargin = TRUE,
    plot_numbers = FALSE,
    plot_text_size = 4,
    plot_breaks = NULL,
    margincolor = "white",
    plot_colorscale = "log10",
    min_abundance = 0.01,
    max_abundance = NULL,
    strip_text = 9,
    xtext_size = 7,
    ytext_size = 11,
    xtext_keep = TRUE,
    xtitle_keep = TRUE,
    grid_clean = TRUE,
    xtext_type_hor = TRUE,
    pheatmap = FALSE,
)


#Create a plot bar for the relative abundance of top 20 microbial genera in Triatoma species based on thier instars
g3$plot_bar(color_values = color_palette_20, others_color = "grey70", facet = "Species", facet2 = "Instar", xtext_keep = FALSE,xtext_size = 7, legend_text_italic = FALSE, barwidth = 1)

#Create a plot bar for the relative abundance of top 20 microbial phyla in Triatoma species based on thier instars
g4$plot_bar(color_values = color_palette_20, others_color = "grey70", facet = "Species", facet2 = "Instar", xtext_keep = FALSE,xtext_size = 7, legend_text_italic = FALSE, barwidth = 1)

#calculate alpha-diversity for samples in dataset "triat1"
triat1$cal_alphadiv(measures = NULL, PD = FALSE)

#calculate beta-diversity for samples in fataset "triat1"
triat1$cal_betadiv(unifrac = FALSE)

#create an object for beta diversity analyses of Triatoma species using bray measure
g5 <- trans_beta$new(dataset = triat1, group = "Species", measure = "bray")

#Generate "NMDS" based ordination of created object "c4"
g5$cal_ordination(ordination = "NMDS")

#Plot the NMDS based ordination for group "Genus"
g5$plot_ordination(plot_color = "Species", plot_shape = "Species", plot_type = c("point", "ellipse")) #NMDS for group "Species"

#Calculate manova within Triatoma "Species" group
g5$cal_manova(manova_all = TRUE)

#Calculate manova for pairs of Triatoma species
g5$cal_manova(manova_all = FALSE)

#Creating an object for alpha diversity analyses of triatoma species
g6 <- trans_alpha$new(dataset = triat1, group = "Species")

#Calculate alpha diversity ofor pairs of Triatoma species based on Dunn's Kruskal-Wallis Multiple Comparisons method
g6$cal_diff(method = "KW_dunn") 

#Plot the prepared alpha diversity based on Dunn's Kruskal-Wallis Multiple Comparisons method and Chao1 measure
g6$res_diff %<>% base::subset(Significance != "ns") #remove "non significant" sign from showing in the plot
g6$plot_alpha(measure = "Chao1", xtext_size = 15) #Plot the alpha diversity 

#Analysis of common OTUs shared between different species of Triatoma, ticks and assassin bugs in individual nests at venn fraction 0.8
#Setting working directory for the analysis
setwd("C:/Users/hassa/Desktop/Convergency/sample_types/")

#Subseting the dataset to include Triatoma, ticks and assassin bugs
sample_types_1000 <- subset_samples(clean_data_1000,Sample_type %in% c("Tick", "Triatoma", "assassin"))

#1-Nest N1 analysis
N1 <-sample_types_1000 %>% subset_samples(Nest == "N1") #Subsetting for nest N1
venn <- ps_venn(N1,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8

#2-Nest N2 analysis
N2 <- sample_types_1000 %>% subset_samples(Nest == "N2") #Subsetting for nest N2
venn <- ps_venn(N2,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8

#3-Nest N4 analysis
N4 <- sample_types_1000 %>% subset_samples(Nest == "N4") #Subsetting for nest N4
venn <- ps_venn(N4,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8

#4-Nest N5 analysis
N5 <- sample_types_1000 %>% subset_samples(Nest == "N5") #Subsetting for nest N5
venn <- ps_venn(N5,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8

#5-Nest N6 analysis
N6 <- sample_types_1000 %>% subset_samples(Nest == "N6") #Subsetting for nest N6
venn <- ps_venn(N6,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8

#6-Nest N7 analysis
N7 <- sample_types_1000 %>% subset_samples(Nest == "N7") #Subsetting for nest N7
venn <- ps_venn(N7,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8

#7-Nest N7d analysis
N7d <- sample_types_1000 %>% subset_samples(Nest == "N7d") #Subsetting for nest N7d
venn <- ps_venn(N7d,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8

#8-Nest N8 analysis
N8 <- sample_types_1000 %>% subset_samples(Nest == "N8") #Subsetting for nest N8
venn <- ps_venn(N8,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8

#9-Nest N12 analysis
N12 <- sample_types_1000 %>% subset_samples(Nest == "N12") #Subsetting for nest N12
venn <- ps_venn(N12,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8

#10-Nest N13 analysis
N13 <- sample_types_1000 %>% subset_samples(Nest == "N13") #Subsetting for nest N13
venn <- ps_venn(N13,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8

#11-Nest N15 analysis
N15 <- sample_types_1000 %>% subset_samples(Nest == "N15") #Subsetting for nest N15
venn <- ps_venn(N15,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8

#12-Nest N16 analysis
N16 <- sample_types_1000 %>% subset_samples(Nest == "N16") #Subsetting for nest N16
venn <- ps_venn(N16,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8

#13-Nest N18 analysis
N18 <- sample_types_1000 %>% subset_samples(Nest == "N18") #Subsetting for nest N18
venn <- ps_venn(N18,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8

#14-Nest N41 analysis
N41 <- sample_types_1000 %>% subset_samples(Nest == "N41") #Subsetting for nest N41
venn <- ps_venn(N41,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8

##Analysis of microbial composition based on primary blood meal of Triatoma and ticks

#Setting working directory for the analysis
setwd("C:/Users/hassa/Desktop/Convergency/blood_meal/")

#Subsetting for samples of triatoma and ticks
triatoma_ticks <- subset_samples(sample_types_1000,Sample_type %in% c("Tick", "Triatoma"))

#Venn analysis of microbial community in Triatoma and ticks samples with Neotoma as primary blood meal for venn fraction 0.9 and 0.8
neotoma <- triatoma_ticks %>% subset_samples(Primary_blood_meal == "Neotoma") #Subsetting for Neotoma
venn <- ps_venn(neotoma,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 1.0) #venn analysis for microbiome dstribution in group "Species" at fraction 1
venn <- ps_venn(neotoma,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.9) #venn analysis for microbiome dstribution in group "Species" at fraction 0.9
venn <- ps_venn(neotoma,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8
venn <- ps_venn(neotoma,"Sample_type", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.5) #venn analysis for microbiome dstribution in group "Species" at fraction 0.5

#Venn analysis of microbial community in Triatoma samples with Dasypys as primary blood meal for venn fraction 1

dasypus <- triatoma %>% subset_samples(Primary_blood_meal == "Dasypus") #Subsetting for Dasypus
venn <- ps_venn(dasypus,"Species", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 1.0) #venn analysis for microbiome dstribution in group "Species" at fraction 1
venn <- ps_venn(dasypus,"Species", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.9) #venn analysis for microbiome dstribution in group "Species" at fraction 0.9
venn <- ps_venn(dasypus,"Species", plot = FALSE ,WEIGHT = TRUE, relative = TRUE, fraction = 0.8) #venn analysis for microbiome dstribution in group "Species" at fraction 0.8

## Analyzing and comparing convergency in microbial communities of Triatoma, ticks and assassin bug species based on:
##1- Decontaminated taxonomic assingment and OTUs data that was filetered for environmental OTUs idnetified at venn fraction 0.5
##2- Entire decontaminated taxonomic assignment and OTUs data

#1- Analysis using Entire decontaminated taxonomic assignment and OTUs data

#Set working directory
setwd ("/Users/evanovakova/Dropbox/Triatomy/2021/Convergency_Anna/Convergency_Rready/FinalConvergency/f05")

#Loading datasheet of OTU table that was filtered for OTUs identified in nest material at venn fraction 0.5
otu_0.5 <-read_excel("/Users/evanovakova/Dropbox/Triatomy/2021/Convergency_Anna/Convergency_Rready/FinalConvergency/f05/filtered_f05_OTU.xlsx")

# Generate a phyloseq object
#1-defining row names from the otu column
otu_0.5 <- otu_0.5 %>%
  tibble::column_to_rownames("OTU") 

#2-transforming otu_0.5 into matrix
otu_0.5 <- as.matrix(otu_0.5)

#3-Transform into phyloseq objects
OTU0.5 = otu_table(otu_0.5, taxa_are_rows = TRUE)
data0.5 <- phyloseq(OTU0.5, TAX, samples, tree)

#Subset the taxa to include only OTUs from bacterial kingdom annd exclude taxonomical assignments with family "mitochondria" and class "Chloroplast".
clean_data0.5 <- data0.5 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Kingdom != "Archaea" &
      Kingdom != "Animalia" &
      Kingdom != "Eukaryota" &
      Kingdom != "ND" &
      Kingdom != "Plantae" &
      Family != "mitochondria" &
      Class != "Chloroplast"
  )

#Set seed 5
set.seed(5)

#Rarfaction at 1000 sample size
clean_data0.5_1000 <- rarefy_even_depth(clean_data0.5, sample.size = 1000)

#Converting phyloseq object to meco dataset
meco_0.5 <- phyloseq2meco(clean_data0.5_1000) #dataset with rarefaction of 1000

make the OTU and sample information in consistent across all files in the dataset object 
meco_0.5$tidy_dataset()

#Setting random number generator to make results reproducible
set.seed(5)

#Clone the dataset rarfied at 1000 sample size
clone0.5 <- clone(meco_0.5)

#Subsetting clone 0.5 for tick, Triatoma and assasin samples

clone0.5$sample_table <- subset(clone0.5$sample_table, Nest!= "ND") 
clone0.5$sample_table <- subset(clone0.5$sample_table, environment!= "Yes") 

#make the OTU and sample information in consistent across all files in clone1
clone0.5$tidy_dataset()

#Create object for convergency analysis
att0.5 <- clone(clone0.5) 

#Meging samples of create object "att" based on group "Species"
att0.5 <- clone0.5$merge_samples(use_group = "Species")

#calculate abundance for taxa
att0.5$cal_abund()

#calculate alpha-diversity for samples
att0.5$cal_alphadiv(measures = NULL, PD = TRUE)

#calculate beta-diversity for samples
att0.5$cal_betadiv(unifrac = TRUE, binary = TRUE)

#generating a trans beta object with measure "wei_unifrac"
a0.5 <- trans_beta$new(dataset = att0.5, measure = "wei_unifrac")

#Plot the distance between group "Species"
a0.5$plot_clustering()

#2- Entire decontaminated taxonomic assignment and OTUs data

#Converting phyloseq object to meco dataset
meco_1000 <- phyloseq2meco(clean_data_1000) #dataset with rarefaction of 1000

#make the OTU and sample information in consistent across all files in the dataset object 
meco_1000$tidy_dataset()

#Setting random number generator to make results reproducible
set.seed(5)

#Clone the dataset rarfied at 1000 sample size
clone0 <- clone(meco_1000)

#Subsetting clone 0 for tick, Triatoma and assasin samples

clone0$sample_table <- subset(clone0$sample_table, Nest!= "ND") 
clone0$sample_table <- subset(clone0$sample_table, environment!= "Yes") 

#make the OTU and sample information in consistent across all files in clone1
clone0$tidy_dataset()

#Create object for convergency analysis
att0 <- clone(clone0) 

#Meging samples of create object "att" based on group "Species"
att0 <- clone0$merge_samples(use_group = "Species")

#calculate abundance for taxa
att0$cal_abund()

#calculate alpha-diversity for samples
att0$cal_alphadiv(measures = NULL, PD = TRUE)

#calculate beta-diversity for samples
att0$cal_betadiv(unifrac = TRUE, binary = TRUE)

#generating a trans beta object with measure "wei_unifrac"
a1 <- trans_beta$new(dataset = att0, measure = "wei_unifrac")

#Plot the distance between group "Species"
a1$plot_clustering()

#Create an object for beta-diversity analysis for group "Sample_types"
att0.1 <- clone(clone0)

#calculate abundance for taxa
att0.1$cal_abund()

#calculate beta-diversity for samples
att0.1$cal_betadiv(unifrac = TRUE, binary = TRUE)

#generating a trans beta object with measure "wei_unifrac"
a1.1 <- trans_beta$new(dataset = att0.1, group = "Sample_type", measure = "wei_unifrac")

#Generate NMDS ordination object a1.1
a1.1$cal_ordination(ordination = "NMDS")

#Plot the generated NMDS ordination
a1.1$plot_ordination(plot_color = "Sample_type", plot_shape = "Sample_type", plot_type = c("point", "ellipse"))

#Calculate manova for group "Sample_type"
a1.1$cal_manova(manova_all = TRUE)

#Calculate manova for pairs within group "Sample_type"
a1.1$cal_manova(manova_all = FALSE)

#Clone the dataset rarfied at 1000 sample size
clone0.1 <- clone(meco_1000)

#Subsetting clone 0 for tick, Triatoma and assasin samples

clone0.1$sample_table <- subset(clone0$sample_table, Nest!= "ND") 
clone0.1$sample_table <- subset(clone0$sample_table, environment!= "Yes")
clone0.1$sample_table <- subset(clone0$sample_table, Diet!= "ND")

#make the OTU and sample information in consistent across all files in clone1
clone0.1$tidy_dataset()

#Create object for beta diversity analysis of group "Sample_types" based on group "Diet"
att0.2 <- clone(clone0.1)

#calculate abundance for taxa
att0.2$cal_abund()

#calculate beta-diversity for samples
att0.2$cal_betadiv(unifrac = TRUE, binary = TRUE)

#generating a trans beta object with measure "wei_unifrac"
a1.2 <- trans_beta$new(dataset = att0.2, group = "Sample_type", measure = "wei_unifrac")

#Calculate manova for group "Sample_type" based on group set "Diet"
a1.2$cal_manova(manova_set = "Diet", manova_all = TRUE) 

## Abundance analysis including the top 20 "Family" taxa within the group "Species"

#Create object for abundance analysis
gbar <- trans_abund$new(dataset = clone0, taxrank = "Family", ntaxa = 20, groupmean = "Species")

#Order the group "Species" output 
order <- c("longipes", "sonoraensis", "hirticornis", "arizonensis", "turicata", "rubida", "protracta", "gerstaeckeri", "lecticularia", "sanguisuga")

#Generate the abundance bar plot
g1 <- gbar$plot_bar(color_values =  color_palette_20, order_x = c("longipes", "sonoraensis", "hirticornis", "arizonensis", "turicata", "rubida", "protracta", "gerstaeckeri", "lecticularia", "sanguisuga"), bar_type = "full", barwidth = NULL, use_alluvium = TRUE, others_color = "grey90", legend_text_italic = FALSE, clustering = FALSE, xtext_type_hor = TRUE)
g1 + theme_classic() + theme(axis.title.y = element_text(size = 18))

##Analysis of assassin bugs

#Clone the dataset rarfied at 1000 sample size
clone1 <- clone(meco_1000)

# Subsetting "clone1" for assassin bugs
clone1$sample_table <- subset(clone1$sample_table, Sample_type == "assassin") 

#Make the OTU and sample information consistent across all files in the data set object
clone1$tidy_dataset()

#calculate alpha-diversity for samples
clone1$cal_alphadiv(measures = NULL, PD = FALSE)

a0 <- trans_alpha$new(dataset = clone1, group = "Species") #creating a trans_alpha object based on group "Species"
a00 <- trans_alpha$new(dataset = clone1, group = "Genus")#creating a trans_alpha object based on group "Genus"
a0$cal_diff(method = "KW_dunn") #Do the pair compirson analysis
a00$cal_diff(method = "KW_dunn") #Do the pair compirson analysis

#View result
View(a0[["res_diff"]])
View(a00[["res_diff"]])

#Plot the prepared alpha diversity based on Dunn's Kruskal-Wallis Multiple Comparisons method and Chao1 measure
a0$res_diff %<>% base::subset(Significance != "ns") #remove "non significant" sign from showing in the plot
a00$res_diff %<>% base::subset(Significance != "ns") #remove "non significant" sign from showing in the plot

a0$plot_alpha(measure = "Shannon", xtext_size = 15) #Plot the alpha diversity 
a00$plot_alpha(measure = "Shannon", xtext_size = 15) #Plot the alpha diversity 

# Create an object that include "clone1" dataset with top 30 "Genus" taxa
a1 <- trans_abund$new(dataset = clone1, taxrank = "Genus", ntaxa = 30)

#Generate heat map for the relative abundance of top 30 genera for nest material with thier destribution on nests and locations
a1$plot_heatmap(
    color_values = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")),
    facet = "Species",
    x_axis_name = NULL,
    order_x = NULL,
    withmargin = TRUE,
    plot_numbers = FALSE,
    plot_text_size = 4,
    plot_breaks = NULL,
    margincolor = "white",
    plot_colorscale = "log10",
    min_abundance = 0.01,
    max_abundance = NULL,
    strip_text = 11,
    xtext_size = 6,
    ytext_size = 11,
    xtext_keep = TRUE,
    xtitle_keep = TRUE,
    grid_clean = TRUE,
    xtext_type_hor = TRUE,
    pheatmap = FALSE,
)

#Calculate beta-diversity for assassin samples
clone1$cal_betadiv(unifrac = FALSE)

#create an object for beta diversity analyses based on group "Genus"
a2 <- trans_beta$new(dataset = clone1, group = "Genus", measure = "bray")

#Generate "NMDS" based ordination of created object "a2"
a2$cal_ordination(ordination = "NMDS")

#Plot the NMDS based ordination for group "Genus"
a2$plot_ordination(plot_color = "Genus", plot_type = c("point", "ellipse")) #NMDS for group "Genus"

#Evaluate statistital signifiscance of found dissimilarities among group "Genus" of assassin bugs based on manova method
a2$cal_manova(
  manova_all = TRUE, #TRUE for the overall test; FALSE represents test for all the paired groups.
  manova_set = NULL,
  group = "Genus",
  p_adjust_method = "fdr",
) #Manova analyis for group "Genus"

a2$res_manova #Check the generated manova statistical result

#Evaluate statistital signifiscance of found dissimilarities among paired groups based on genus of assassin bugs using manova method
a2$cal_manova(
  manova_all = FALSE, #TRUE for the overall test; FALSE represents test for all the paired groups.
  manova_set = NULL,
  group = "Genus",
  p_adjust_method = "fdr",
) #Manova analyis for group "Genus"


#create an object for beta diversity analyses based on group "Species"
a2.5 <- trans_beta$new(dataset = clone1, group = "Species", measure = "bray")

#Generate "NMDS" based ordination of created object "a2.5"
a2.5$cal_ordination(ordination = "NMDS")

#Plot the NMDS based ordination for group "Species"
a2.5$plot_ordination(plot_color = "Species", plot_type = c("point", "ellipse")) #NMDS for group "Genus"

#Evaluate statistital signifiscance of found dissimilarities among group "Species" of assassin bugs based on manova method
a2.5$cal_manova(
  manova_all = TRUE, #TRUE for the overall test; FALSE represents test for all the paired groups.
  manova_set = NULL,
  group = "Species",
  p_adjust_method = "fdr",
) #Manova analyis for group "Species"

a2.5$res_manova #Check the generated manova statistical result

#Evaluate statistital signifiscance of found dissimilarities among paired groups based on species of assassin bugs using manova method
a2.5$cal_manova(
  manova_all = FALSE, #TRUE for the overall test; FALSE represents test for all the paired groups.
  manova_set = NULL,
  group = "Species",
  p_adjust_method = "fdr",
) #Manova analyis for group "Species"

#Analysis of common OTUs shared between different Genus of assassin bugs at venn fraction 0.3
assassin <- clean_data_1000 %>% subset_samples(Sample_type == "assassin") #Subset phyloseq object "cleandata_1000" to include only assassin bugs

#Venn analysis of common OTUs between assassin species
venn1<-ps_venn(assassin,"Genus", plot = FALSE , WEIGHT = TRUE, relative = TRUE, fraction = 0.3) #Generate a list of common and unique OTUs between assassin Species at fraction 0.3
venn2<-ps_venn(assassin,"Genus", plot = FALSE , WEIGHT = TRUE, relative = TRUE, fraction = 0.5) #Generate a list of common and unique OTUs between assassin Species at fraction 0.5
venn3<-ps_venn(assassin,"Genus", plot = TRUE , WEIGHT = TRUE, relative = TRUE, fraction = 0.5, quantities = list(type=c("percent","counts"), font = 1), labels = list(cex = 1.5), col = "black", fill = c("#003E1F","#A89B8C","#F0DFAD","#8F5C38")) #Generate venn diagram for microbiome of assassins based on group "Genus"

##Analysis of ticks (Ornithodoros turicata) microbiome

#Clone the dataset rarfied at 1000 sample size
clone2 <- clone(meco_1000)

#Subsetting "clone2" for ticks
clone2$sample_table <- subset(clone2$sample_table, Sample_type == "Tick") 

#Make the OTU and sample information consistent across all files in the data set object
clone2$tidy_dataset()

##calculate abundance for taxa
clone2$cal_abund()

#Create object for abundance analysis
a3 <- trans_abund$new(dataset = clone2, taxrank = "Genus", ntaxa = 30)

#Generate heat map for abundance of top 30 Genus in ticks based on thier destribution in "State"
a3$plot_heatmap(facet = "State", xtext_keep = TRUE, xtext_type_hor = FALSE,xtext_size = 5, withmargin = FALSE)

#Calculate alpha-diversity for O.turicata samples in clone2 dataset
clone2$cal_alphadiv(measures = NULL, PD = FALSE)

#Alpha diversity analysis between ticks based on group "Location"
a4 <- trans_alpha$new(dataset = clone2, group = "Location") #Create a trans alpha object 
a4$cal_diff(method = "KW_dunn") #testing alpha diversity within samples in pairs of group "Location" using Dunn's Kruskal-Wallis Multiple Comparison method
a4$res_diff #Checking KW_dunn analysis results

#Remove non-significant (ns) sign from showing in the plot
a4$res_diff %<>% base::subset(Significance != "ns")

#Plot the alpha diversity results
a4$plot_alpha(
  color_values = RColorBrewer::brewer.pal(8, "Dark2"),
  measure = "Shannon",
  add_sig = TRUE,
  add_sig_label = "Significance",
  add_sig_text_size = 3.88,
  use_boxplot = TRUE,
  boxplot_color = TRUE,
  boxplot_add = "jitter",
  order_x_mean = FALSE,
  xtext_angle = 45,
  xtext_size = 10,
  ytitle_size = 12,
)

Calculate beta diversity of "clone2"dataset
clone2$cal_betadiv(unifrac = TRUE, binary = TRUE)

#Clusterning analysis for ticks dataset "clone2" with group "Location" based on Bray-Curtis distance
a5 <- trans_beta$new(dataset = clone2, group = "Location", measure = "bray") #Create a trans beta object
a5$cal_ordination(ordination = "PCoA") #Use PCoA ordination in the clustering analysis
a5$plot_ordination(plot_color = "Location", plot_shape = "Location", plot_type = c("point", "ellipse")) #plot the ordination

#Calculating sample distances within tick groups "Location"
a5$cal_group_distance(within_group = TRUE) #Perform Wilcoxon Rank Sum test

#Perform Wilcoxon Rank Sum test
a5$cal_group_distance_diff(method = "wilcox")

#Plot the distance of tick samples based on "Location" group
a5$plot_group_distance(boxplot_add = "mean") #generate the plot with mean value added in each boxplot

#Evaluate statistital significance of found dissimilarities among tick samples wthin group "Location" based on manova method
a5$cal_manova(
  manova_all = TRUE, #TRUE for the overall test; FALSE represents test for all the paired groups.
  manova_set = NULL,
  group = "Location",
  p_adjust_method = "fdr",
)

#View manova result
a5$res_manova

#Clusterning analysis for ticks dataset "clone2" with group "Location" based on Wei_unifrac measure 
a6 <- trans_beta$new(dataset = clone2, group = "Location", measure = "wei_unifrac")
a6$cal_ordination(ordination = "PCoA") #Use PCoA ordination in the clustering analysis
a6$plot_ordination(plot_color = "Location", plot_shape = "Location", plot_type = c("point", "ellipse")) #plot the ordination

#Calculating sample distances within tick groups "Location"
a6$cal_group_distance(within_group = TRUE) #Perform Wilcoxon Rank Sum test

#Perform Wilcoxon Rank Sum test
a6$cal_group_distance_diff(method = "wilcox")

#Plot the distance of tick samples based on "Location" group
a6$plot_group_distance(boxplot_add = "mean") #Using Wilcox test as default

#Evaluate statistital significance of found dissimilarities among tick samples wthin group "Location" based on manova method
a6$cal_manova(
  manova_all = TRUE, #TRUE for the overall test; FALSE represents test for all the paired groups.
  manova_set = NULL,
  group = "Location",
  p_adjust_method = "fdr",
)

#View manova result
a6$res_manova
