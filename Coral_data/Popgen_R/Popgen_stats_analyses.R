# Popgen and other statistical analyses associated with the manuscript:
# Rivera et al. 2022 - Palau’s warmest reefs harbor a thermally tolerant coral lineage that thrives across different habitats

## Setup ####
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ade4)
library(adegenet)
library(mmod)
library(hierfstat) # has a pairwise.fst function also, the one I use, if there are errors specify package
library(ape)
library(StAMPP)
library(pegas)
library(poppr)
library(Rcpp)
library(vcfR)
library(factoextra)
library(reshape2)
library(vegan)
library(viridis)
library(stringr)
library(PopGenReport)
library(ggpubr)

setwd("/Users/hannyrivera/Documents/MIT/Research_Papers/Palau_Microsat_Paper/Drafts/Comms_Bio_submission/Palau_porites/Coral_data/Popgen_R/")

## Importing rad data ####
# See snp_filtering.sh for step taken to produce this final vcf
rad_vcf<-read.vcfR("/input_files/Allsamps_maf05_final_snps.vcf")

# Convert vcf to genind and genlight objects for using adegenet package
rad_gl<-vcfR2genlight(rad_vcf, n.cores=5)
rad_gin<-vcfR2genind(rad_vcf)

# Import population strata data frame
# Read in master strata file
metadata<-read.table("/input_files/Strata_all_samples.txt", header=TRUE)

rad_samples<-as.data.frame(indNames(rad_gin))
colnames(rad_samples)<-c("Sample")
#convert to character to properly compare using all.equal below
rad_samples$Sample<-as.character(rad_samples$Sample)
#now merge in the info that corresponds with the samples in the current data set 
#only grabs info for samples that are in the rad datatypes 
rad_strata<-left_join(rad_samples,metadata, by="Sample")
#make sure the samples are all in order 
all.equal(rad_samples$Sample, rad_strata$Sample) 

# Including the structure results into the strata data frame
rad_K4<-read.table("/input_files/rad_K4_results.txt", header=TRUE)

rad_K4$K4_radpop<-colnames(rad_K4[,2:5])[max.col(rad_K4[,2:5],ties.method="first")]

# Fix a few samples that are almost perfectly admixed and assign to a cluster that makes more sense based on the PCA: 
# Palau236 and Palau265 are each 0.51 LB 0.49 RD but clearly very far away from the light blue cluster on the PCA 
rad_K4$K4_radpop[43]<-"RD"
rad_K4$K4_radpop[46]<-"RD"

rad_strata<-left_join(rad_strata, rad_K4, by="Sample")


strata(rad_gl)<-rad_strata
strata(rad_gin)<-rad_strata

setPop(rad_gl)<-~K4_radpop 
setPop(rad_gin)<-~K4_radpop

popNames(rad_gl) #check assignments
popNames(rad_gin)

#Check metadata
rad_gl
rad_gin

# Table S2
# select(metadata, -Cohort, -Lipid, -BoringVol, -TissThick)->tableS2
# write.csv(tableS2, file = "tableS2.csv",  
#           quote=FALSE, na = "-", row.names = FALSE)


## Import microsat data ####
msat_genid<-read.structure("/Users/hannyrivera/Documents/MIT/Research_Papers/Palau_Microsat_Paper/Drafts/Comms_Bio_submission/Palau_porites/Coral_data/Popgen_R/input_files/Palau_12M_input.stru") 
# 322 genos, 12 markers, col 1 genotypes, col 2 pop, no others, row 1 markers, single row yes 

#make strata for msat dataset 
msat_strata<-as.data.frame(indNames(msat_genid))
colnames(msat_strata)<-"Sample"
#pull metadata for samples with microsats
msat_strata<-left_join(msat_strata,metadata, by="Sample")

# Also include the STRUCTURE assignments for K=4 from microsatellite runs 
msat_K4<-read.table("/Users/hannyrivera/Documents/MIT/Research_Papers/Palau_Microsat_Paper/Drafts/Comms_Bio_submission/Palau_porites/Coral_data/Popgen_R/input_files/msat_K4_results.txt", header = TRUE)
msat_K4$K4_msatpop<-colnames(msat_K4[,2:5])[max.col(msat_K4[,2:5],ties.method="first")]

msat_strata<-left_join(msat_strata,select(msat_K4,Sample, K4_msatpop), by="Sample")
#set pop to lineages
strata(msat_genid)<-msat_strata
setPop(msat_genid)<-~K4_msatpop

# Incorporate K4 structure lineages to metadata 
metadata$K4_msatpop<-msat_strata[["K4_msatpop"]][match(metadata$Sample, msat_strata[["Sample"]])]
metadata$K4_radpop<-rad_strata[["K4_radpop"]][match(metadata$Sample, rad_strata[["Sample"]])]

metadata$K4_msatpop<-as.factor(metadata$K4_msatpop)
metadata$K4_radpop<-as.factor(metadata$K4_radpop)

# Combine across rad and msat assignments for downstream plotting of core data by lineage
metadata%>%mutate(final_lineage = case_when(is.na(K4_radpop) ~ K4_msatpop, TRUE ~ K4_radpop))->metadata

## Save initial data ####
#save.image(file="/Users/hannyrivera/Documents/MIT/Research_Papers/Palau_Microsat_Paper/Drafts/Comms_Bio_submission/Palau_porites/Coral_data/Popgen_R/Setup_genetic_data.RData")


## Lineages by site plots  (Figure 2) ####
# Calculate frequency of each lineage by site 
select(metadata, "Site", "final_lineage")%>%subset(!is.na(final_lineage))->out
out%>%group_by(Site)%>%count(final_lineage)%>%mutate(Freq=n/sum(n))->site_lineages
site_lineages<-as.data.frame(site_lineages)
#Some sites are missing lineages so just adding them in as zeros into the data frame
#this only matters if the site if missing something other than RD because that one is alphabetically last 
# if it's missing it just won't get filled which is fine, but otherwise if an earlier lineage is missing
# the colors will be wrong
rbind(site_lineages, data.frame("Site"= c("Mecherchar", "Ngermid", "Ngerchelong"), 
                                "final_lineage" = c("PI","PI", "LB"), 
                                "n" = c(0,0,0), 
                                "Freq"=c(0,0,0)))->site_lineages

for (reef in levels(site_lineages$Site)){
  subset(site_lineages, Site==reef)->test
  test$ymax <- cumsum(test$Freq)
  # Compute the bottom of each rectangle
  test$ymin <- c(0, head(test$ymax, n=-1))
  # Make the plot
  assign(paste("pie", reef, sep=""),ggplot(test, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=final_lineage)) +
           geom_rect() +
           scale_fill_manual(values=c("darkblue","cyan3","lightcoral","red2"))+
           coord_polar(theta="y") +
           xlim(c(2, 4)) +
           theme_void()+
           theme(legend.position = "none"))
  ggsave(paste(reef, ".png", sep=""), get(paste("pie", reef, sep="")), width=5, height=5, units="in", dpi=300, bg="transparent")
}




## Plot STRUCTURE assignments for RAD data (Figure 3A) ####
ggplot(rad_K4_meta_long, aes(x=Sample, y=Assignment, fill=Lineage))+geom_bar(stat="identity")+
  scale_fill_manual(values=c("darkblue","cyan3","lightcoral","red2"))+
  facet_grid(~factor(Site, levels=c('Mecherchar','Mecherchar_Channel', 'Taoch', 'Ngermid', 'Risong', 'Helen',
                                    'Ngelsible','Outer_Taoch', 'Ngerdiluches', 'Drop_Off','Ngerchelong', 'Melekeok', 'Kayangel')), 
             scales="free", space="free", 
             labeller = as_labeller(c(`Drop_Off`="Drop Off",`Helen` = "Helen",`Kayangel`= "Kayangel",`Mecherchar` = "Mecherchar",`Mecherchar_Channel`= "Mech. Channel", 
                                      `Melekeok` = "Melekeok",`Ngelsible` ="Ngelsible",`Ngerchelong` ="Ngerch", `Ngerdiluches`  ="Ngerd.",
                                      `Ngermid` ="Ngermid", `Outer_Taoch` ="Outer Taoch",`Risong`= "Risong",`Taoch` ="Taoch")) )+
  theme_bw()+guides(fill="none")+xlab("")+
  theme(axis.title = element_text(face="bold"),
        axis.text = element_text(face="bold"), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text =element_text(face="bold", size=6),
        panel.border = element_rect(colour="black", size=0.5),
        panel.spacing.x = unit(0.1, "lines")) 

ggsave("Figure3A.png", width=12, height=4, units="in", dpi=300)  

## RAD Principal components analyses and DAPC (Figure 3B-C) ####
#Calculate PCA
rad_scale<-scaleGen(rad_gin, NA.method="mean")
rad_pca<-dudi.pca(rad_scale, cent=FALSE, scale=FALSE, scannf = FALSE, nf = 3)
# gets the percent variance explained by the PCs
rad_pca_eig<-get_eigenvalue(rad_pca)
rad_pca_df<-rad_pca$li
rad_pca_df$Sample<-rownames(rad_pca_df)
rad_pca_df<-left_join(rad_pca_df,rad_strata, by="Sample")

## Figure 3B
ggplot(rad_pca_df, aes(x=Axis1, y=Axis2,colour=K4_radpop))+
          geom_point(size=0.85)+stat_ellipse(geom="polygon", aes(fill=K4_radpop),alpha=0.25)+
          scale_colour_manual(values=c("darkblue","powderblue", "lightcoral", "red1"))+
          scale_fill_manual(values=c("darkblue","powderblue", "lightcoral", "red1"))+
          guides(colour="none", fill="none")+
          xlab(paste("PC1 (", round(eval(rad_pca_eig$variance.percent[1]),2), "%)", sep=""))+
          ylab(paste("PC1 (", round(eval(rad_pca_eig$variance.percent[2]),2), "%)", sep=""))+
          theme_bw()+theme(axis.title = element_text(face="bold"),
                           axis.ticks = element_blank(),
                           axis.text = element_blank())
ggsave("Figure3B.png", width=4, height=4, units="in", dpi=300)  

#DAPC of results 
rad_clust<-find.clusters(rad_gin, n.pca = 140, n.clust = 4)

grp <- pop(rad_gin)
xval <- xvalDapc(rad_scale, grp, n.pca.max = 300, training.set = 0.9,result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 20, xval.plot = TRUE)
# 70 PCs achieves lowest MSE 
rad_dapc<-dapc(rad_gin, rad_clust$grp, n.pca = 70, n.da = 3)

# The colors will have to be reordered if rerun since it might cluster things in a different order
myCol <- c("cyan3","lightcoral","darkblue", "red1") #These are matched in order based on: 
#table(pop(rad_gin), rad_clust$grp) # double check order before plotting 

#2D DAPC (Figure3C) 
pdf("Figure3C.pdf", width=5, height=4)
scatter(rad_dapc, scree.da=FALSE, bg="white", pch=19, clab=0,
        cstar=0, col=myCol, scree.pca=TRUE, posi.pca="topright")
dev.off()

#1D DAPC (Figure 3C insert)
pdf("Figure3C_insert.pdf", width=5, height=4)
scatter(rad_dapc,1,1, col=myCol, bg="white",scree.da=FALSE, legend=FALSE, solid=.4)
dev.off()


## RAD Fst (Figure 3D) ####

# pairwise fst between struct pops 
rad_fst<-gl.fst.pop(rad_gl, nboots=100, percent = 95, nclusters=3)
# plot as heat map
rad_fst_pl<-data.frame("L1"=c(rep("DB", 3), rep ("LB",2), "PI"), 
                       "L2"= c("LB", "PI", "RD", "PI", "RD", "RD"),
                       "Fst"=c(0.24,0.67,0.40,0.64,0.35, 0.70))

ggplot(rad_fst_pl, aes(x=L1, y=L2, label=as.character(Fst)))+
  geom_tile(mapping=aes(fill=Fst))+
  scale_fill_gradient(low="thistle1", high="plum3")+
  geom_text(size=3, fontface="bold")+
  guides(fill=FALSE)+
  theme_bw()+
  theme(axis.text= element_text(face="bold"),
        axis.title=element_blank(),
        panel.grid = element_blank())

ggsave("Fig3D.png", width=3, height=3, units="in", dpi=300)



## Plot STRUCTURE assignments for msat data (Figure S1A) ####
ggplot(msat_K4_data_long, aes(x=Sample, y=Msat_val, fill=msat_group))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("darkblue","cyan3","lightcoral","red2"))+
  facet_grid(~factor(Site, levels=c('Mecherchar','Mecherchar_Channel', 'Taoch', 'Ngermid', 'Risong', 'Helen',
                                    'Ngelsible','Outer_Taoch', 'Ngerdiluches', 'Drop_Off','Ngerchelong', 'Melekeok', 'Kayangel')), 
             scales="free", space="free", 
             labeller = as_labeller(c(`Drop_Off`="Drop Off",`Helen` = "Helen",`Kayangel`= "Kayangel",`Mecherchar` = "Mecherchar",`Mecherchar_Channel`= "Mech. Channel", 
                                      `Melekeok` = "Melekeok",`Ngelsible` ="Ngelsible",`Ngerchelong` ="Ngerch", `Ngerdiluches`  ="Ngerd.",
                                      `Ngermid` ="Ngermid", `Outer_Taoch` ="Outer Taoch",`Risong`= "Risong",`Taoch` ="Taoch")) )+
  theme_bw()+
  guides(fill=FALSE)+
  xlab("")+ylab("Assignment")+
  theme(axis.title = element_text(face="bold"),
        axis.text = element_text(face="bold"), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text =element_text(face="bold", size=6),
        panel.border = element_rect(colour="black", size=0.5),
        panel.spacing.x = unit(0.1, "lines"))
ggsave("FigureS1A.png", width=12, height=3, units="in", dpi=300)  

## Msat PCA/DAPC (Figure S1B-C) ####
msat_scale<-scaleGen(msat_genid, NA.method="mean")
msat_pca<-dudi.pca(msat_scale, cent=FALSE, scale=FALSE, scannf = FALSE, nf = 3)
# gets the percent variance explained by the PCs
msat_pca_eig<-get_eigenvalue(msat_pca)
msat_pca_df<-msat_pca$li
msat_pca_df$Sample<-rownames(msat_pca_df)
msat_pca_df<-left_join(msat_pca_df,msat_strata, by="Sample")

FigS1B<-(ggplot(msat_pca_df, aes(x=Axis1, y=Axis2,colour=K4_msatpop))+
           geom_point(size=0.85)+stat_ellipse(geom="polygon", aes(fill=K4_msatpop),alpha=0.25)+
           scale_colour_manual(values=c("darkblue","cyan3", "lightcoral", "red1"))+
           scale_fill_manual(values=c("darkblue","cyan3", "lightcoral", "red1"))+
           guides(colour="none", fill="none")+
           xlab(paste("PC1 (", round(eval(msat_pca_eig$variance.percent[1]),2), "%)", sep=""))+
           ylab(paste("PC2 (", round(eval(msat_pca_eig$variance.percent[3]),2), "%)", sep=""))+
           theme_bw()+
           theme(axis.title = element_text(face="bold"),
                 axis.ticks = element_blank(),
                 axis.text = element_blank()))
#ggsave("FigS1B.png", width=5, height=4, units="in", dpi=300)

##DAPC 
xval_msat <- xvalDapc(msat_scale, pop(msat_genid), n.pca.max = 300, training.set = 0.9,
                      result = "groupMean", center = TRUE, scale = FALSE,
                      n.pca = NULL, n.rep = 30, xval.plot = TRUE)
#50 PCs achieve the lowest MSE

msat_clust<-find.clusters(msat_genid,n.pca = 300, n.clust = 4)
msat_dapc<-dapc(msat_genid, msat_clust$grp, n.pca = 50, n.da = 3)

myCol <- c( "red1","darkblue","cyan3","lightcoral") #These are matched in order based on: 
# table(pop(msat_genid), msat_dapc$grp) 
#Same comment here as the RAD section above. The colors may need reording based on the way the clusters order themselves

#2D DAPC (Fig S1C)
pdf("msat_dpc.pdf", width=5, height=4)
scatter(msat_dapc, scree.da=FALSE, bg="white", pch=19, clab=0,
        cstar=0, col=myCol, scree.pca=TRUE, posi.pca="bottomleft")
dev.off()

#1D DAPC (DF1) (Fig S1C insert)
pdf("/Users/hannyrivera/Documents/MIT/Research_Papers/Palau_Microsat_Paper/Post_review_reanalyses/git_repo/Cleaned/Figures/Msat_popgen//msat_dpc_1dim.pdf", width=5, height=4)
scatter(msat_dapc,1,1, col=myCol, bg="white",scree.da=FALSE, legend=FALSE, solid=.4)
dev.off()

## Msat Fst (Figure S1D) ####

msat_fst<-genet.dist(msat_genid, method="Nei87")
# plot as heat map
msat_fst_pl<-data.frame("L1"=c(rep("DB", 3), rep ("LB",2), "PI"), 
                        "L2"= c("LB", "PI", "RD", "PI", "RD", "RD"),
                        "Fst"=c(0.07,0.16,0.10,0.14,0.10, 0.16))

ggplot(msat_fst_pl, aes(x=L1, y=L2, label=as.character(Fst)))+
  geom_tile(mapping=aes(fill=Fst))+
  scale_fill_gradient(low="thistle1", high="plum3")+
  geom_text(size=3, fontface="bold")+
  guides(fill=FALSE)+
  theme_bw()+
  theme(axis.text= element_text(face="bold"),
        axis.title=element_blank(),
        panel.grid = element_blank())

ggsave("FigS1D.png", width=3, height=3, units="in", dpi=300)


## Compare Msat-RAD STRUCTURE results (Figure S2) ####
#setup data frame with both the msat and rad structure admixture results

msat_K4_data_long<-left_join(msat_K4, select(metadata, Sample, Site, Region), by="Sample")
colnames(msat_K4_data_long)[2:5]<-c("DB_msat","LB_msat","RD_msat","PI_msat")

rad_K4_copy<-rad_K4
colnames(rad_K4_copy)[2:5]<-c("LB_rad","PI_rad","DB_rad","RD_rad")

struct_comp<-inner_join(msat_K4_data_long, rad_K4_copy, by="Sample")

pivot_longer(struct_comp, cols=c("PI_msat", "RD_msat","DB_msat","LB_msat", "LB_rad","PI_rad","DB_rad","RD_rad"), names_to = c("Lineage", "DataType"), names_sep = "_", values_to = "Assignment")->struct_comp_long

## Make plot of assignments for each (Figure S2A)
FigS2A<-ggplot(struct_comp_long, aes(x=Sample, y=Assignment, fill=Lineage))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("darkblue","cyan3","lightcoral","red2"), na.value="grey80")+
  facet_wrap(~DataType,nrow = 2, ncol=1, 
             labeller = as_labeller(c(`msat`= "Microsatellites",`rad` = "RAD_seq")))+
  theme_bw()+guides(fill="none")+xlab("")+
  theme(axis.title = element_text(face="bold"),
        axis.text = element_text(face="bold"), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text =element_text(face="bold"),
        panel.grid = element_blank())

FigS2B<-ggplot(metadata, aes(x=K4_radpop, y=K4_msatpop, colour=K4_radpop))+
  scale_colour_manual(values=c("darkblue","cyan3","lightcoral","red2"), na.value="grey80")+
  geom_jitter(width=0.1, height=0.1)+
  guides(colour="none")+xlab("RAD Population")+
  ylab("Msat Population")+theme_bw()+
  theme(axis.title = element_text(face="bold"),
        axis.text = element_text(face="bold"))

ggarrange(FigS2A, FigS2B, ncol=1, nrow=2, labels = c("A", "B"))

## Core data by lineage (Figure 4) ####

# Skeletal Density
dens_tuk<-TukeyHSD(x=aov(lm(metadata$Density~metadata$final_lineage)), 'metadata$final_lineage', conf.level=0.95)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# metadata$final_lineage  3 0.4776 0.15919    9.04 4.07e-05 ***
# Residuals              68 1.1975 0.01761                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

dens_tuk # RD sig from all others
# $`metadata$final_lineage`
#              diff         lwr         upr     p adj
# LB-DB -0.08148990 -0.17607208  0.01309228 0.1155268
# PI-DB  0.04972222 -0.10439660  0.20384104 0.8304373
# RD-DB -0.24569444 -0.38230618 -0.10908271 0.0000667
# PI-LB  0.13121212 -0.02975974  0.29218398 0.1489216
# RD-LB -0.16420455 -0.30850317 -0.01990592 0.0194375
# RD-PI -0.29541667 -0.48417290 -0.10666043 0.0005924

Fig4A<-ggplot(subset(metadata,!is.na(final_lineage)), 
           aes(x=final_lineage, y=Density, colour=final_lineage))+
    geom_boxplot(na.shape=NA)+
    geom_jitter(height=0, width=0.1)+
    scale_colour_manual(values=c("darkblue","cyan3","lightcoral", "red2"))+
    ylab(expression(Skeletal~Density~(g~cm^-3)))+
    xlab("")+
    annotate(geom="text", x=4, y=1.4, label="*", size=10)+
    theme_bw()+
    guides(colour=FALSE)+
    theme(axis.text.y = element_text(face="bold"),
          axis.title.y = element_text(face="bold"),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank())

# Calcification
calc_tuk<-TukeyHSD(x=aov(lm(metadata$CalcRate~metadata$final_lineage)), 'metadata$final_lineage', conf.level=0.95)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# metadata$final_lineage  3  2.292  0.7639   5.959 0.00114 **
# Residuals              68  8.717  0.1282   

calc_tuk # RD sig from DB and LB, p =0.07 for PI
# $`metadata$final_lineage`
#              diff        lwr         upr     p adj
# LB-DB -0.06815657 -0.3233392  0.18702605 0.8953851
# PI-DB -0.11194444 -0.5277568  0.30386796 0.8931790
# RD-DB -0.58861111 -0.9571894 -0.22003281 0.0004438
# PI-LB -0.04378788 -0.4780898  0.39051401 0.9933856
# RD-LB -0.52045455 -0.9097721 -0.13113701 0.0042074
# RD-PI -0.47666667 -0.9859308  0.03259744 0.0747279

Fig4B<-ggplot(subset(metadata,!is.na(final_lineage)), 
          aes(x=final_lineage, y=CalcRate, colour=final_lineage))+
  geom_boxplot(na.shape=NA)+
  geom_jitter(height=0, width=0.1)+
  scale_colour_manual(values=c("darkblue","cyan3","lightcoral", "red2"))+
  ylab(expression(Calcification~(g~cm^-2~y^-1)))+
  xlab("")+
  ylim(c(0,2.4))+
  geom_segment(aes(x = 1, y = 2.25, xend = 4, yend = 2.25), colour="black")+
  geom_segment(aes(x = 2, y = 1.9, xend = 4, yend = 1.9), colour="black")+               
  annotate(geom="text", x=3, y=1.95, label="*", size=10)+
  annotate(geom="text", x=2.5, y=2.3, label="*", size=10)+
  theme_bw()+
  guides(colour=FALSE)+
  theme(axis.text.y = element_text(face="bold"),
        axis.text.x=element_blank(),
        axis.title.y = element_text(face="bold"),
        axis.ticks.x = element_blank())



ext_tuk<-TukeyHSD(x=aov(lm(metadata$ExtRate~metadata$final_lineage)), 'metadata$final_lineage', conf.level=0.95)
#                       Df Sum Sq Mean Sq F value Pr(>F)  
# metadata$final_lineage  3  0.897 0.29884   3.102 0.0323 *
# Residuals              68  6.550 0.09632    

ext_tuk # RD diff from DB and LB 
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# Fit: aov(formula = lm(metadata$ExtRate ~ metadata$final_lineage))
# $`metadata$final_lineage`
# diff        lwr         upr     p adj
# LB-DB  0.008813131 -0.2123882  0.23001444 0.9995823
# PI-DB -0.125277778 -0.4857187  0.23516311 0.7967123
# RD-DB -0.345277778 -0.6647745 -0.02578108 0.0291333
# PI-LB -0.134090909 -0.5105591  0.24237732 0.7845046
# RD-LB -0.354090909 -0.6915651 -0.01661671 0.0361161
# RD-PI -0.220000000 -0.6614481  0.22144813 0.5582455

Fig4C<-ggplot(subset(metadata,!is.na(final_lineage)), 
          aes(x=final_lineage, y=ExtRate, colour=final_lineage))+
  geom_boxplot(na.shape=NA)+
  geom_jitter(height=0, width=0.1)+
  scale_colour_manual(values=c("darkblue","cyan3","lightcoral", "red2"))+
  ylab(expression(Extension~Rate~(cm~y^-1)))+
  xlab("")+
  ylim(c(0,2.4))+
  geom_segment(aes(x = 1, y = 2.25, xend = 4, yend = 2.25), colour="black")+
  geom_segment(aes(x = 2, y = 1.9, xend = 4, yend = 1.9), colour="black")+               
  annotate(geom="text", x=3, y=1.95, label="*", size=10)+
  annotate(geom="text", x=2.5, y=2.3, label="*", size=10)+
  theme_bw()+
  guides(colour=FALSE)+
  theme(axis.text.y = element_text(face="bold"),
        axis.text.x=element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        axis.ticks.x = element_blank())

subset(metadata, !is.na(final_lineage))%>%group_by(final_lineage)%>%
  summarise(Freq98SB=sum(SB_98,na.rm=TRUE)/length(na.omit(SB_98)),
            Freq10SB=sum(SB_10,na.rm=TRUE)/length(na.omit(SB_10)),
            Count98SB=sum(SB_98,na.rm=TRUE), 
            Count10SB=sum(SB_10,na.rm=TRUE),
            TotalN_98=length(na.omit(SB_98)),
            TotalN_10=length(na.omit(SB_10)))%>%
  mutate(Count98NoSB=TotalN_98-Count98SB,
         Count10NoSB=TotalN_10-Count10SB)->SB_data

#Chi-squared test for 1998 stress bands
prop.test(SB_data$Count98SB, SB_data$TotalN_98)
# 4-sample test for equality of proportions without continuity correction
# 
# data:  SB_data$Count98SB out of SB_data$TotalN_98
# X-squared = 12.349, df = 3, p-value = 0.006278
# alternative hypothesis: two.sided
# sample estimates:
#   prop 1    prop 2    prop 3    prop 4 
# 0.6800000 0.2222222 0.1666667 0.2500000 

pivot_longer(SB_data, cols = c(Count98SB, Count10SB,Count98NoSB,Count10NoSB), names_to = "Band_type", values_to = "Band_counts")->SB_counts
SB_counts%>%mutate(year = case_when(c(Band_type=="Count98SB" |Band_type=="Count98NoSB") ~ "1998", TRUE ~ "2010"))->SB_counts

Fig4D<-ggplot(subset(SB_counts, year==1998 & Band_type=="Count98SB"), aes(x=final_lineage, y=Freq98SB, fill=final_lineage, colour=final_lineage))+
    geom_bar(stat="identity", aes(alpha=0.9))+
    scale_fill_manual(values=c("darkblue","cyan3","lightcoral", "red2"))+
    scale_colour_manual(values=c("darkblue","cyan3","lightcoral", "red2"))+
    theme_minimal()+ylab("Stress Band Prevalence")+xlab("")+
    theme_bw()+theme(panel.grid.minor=element_blank(), 
                     panel.grid.major.x = element_blank(),
                     axis.text = element_text(face="bold"),
                     axis.text.x = element_text(face="bold"), 
                     axis.ticks.x = element_blank())+
    guides(colour=FALSE, fill=FALSE, alpha=FALSE)+
    scale_y_continuous(limits= c(0,0.75), breaks=c(0,0.25,0.5,0.75), labels=c("0%", "25%", "50%", "75%"))+
    annotate(geom="text", label="*", x=1, y=0.73, size=10, fontface="bold")+
    annotate(geom="text", label="N=25", x=1, y=0.025, fontface="italic")+
    annotate(geom="text", label="N=18", x=2, y=0.025, fontface="italic")+
    annotate(geom="text", label="N=6", x=3, y=0.025, fontface="italic")+
    annotate(geom="text", label="N=8", x=4, y=0.025, fontface="italic")

ggarrange(Fig4A,Fig4B,Fig4C,Fig4D, nrow=2, ncol=2, labels = c("A", "B", "C", "D"))
ggsave("Figure4.png", width=6, height=6, units="in")

## Core data by for DB and LB lineages by Region (Figure 5) ####
subset(metadata, final_lineage=="DB" | final_lineage =="LB")->metadata_dblb
droplevels(metadata_dblb)->metadata_dblb

Fig5A<-ggplot(subset(metadata_dblb,!is.na(final_lineage)), 
            aes(x=final_lineage, y=Density, colour=final_lineage))+
    geom_boxplot(na.shape=NA)+
    geom_jitter(height=0, width=0.1)+
    scale_colour_manual(values=c("darkblue","cyan3"))+
    ylab(expression(Skeletal~Density~(g~cm^-3)))+
    xlab("")+
    #annotate(geom="text", x=4, y=1.4, label="*", size=10)+
    theme_bw()+
    guides(colour=FALSE)+
    theme(axis.text.y = element_text(face="bold"),
          axis.title.y = element_text(face="bold"),
          axis.text.x=element_text(face="bold"))+
    facet_wrap(~Region, labeller = as_labeller(c(`OR`= "Outer Reefs",`RI` = "Rock Islands")))

t.test(subset(metadata, final_lineage=="DB" & Region=="OR")$Density,subset(metadata, final_lineage=="DB" & Region=="RI")$Density)

# data:  subset(metadata, final_lineage == "DB" & Region == "OR")$Density and subset(metadata, final_lineage == "DB" & Region == "RI")$Density
# t = 6.0749, df = 12.18, p-value = 5.203e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  0.1517095 0.3209572
# sample estimates:
# mean of x mean of y 
#  1.286333  1.050000 


Fig5B<-ggplot(subset(metadata_dblb,!is.na(final_lineage)), 
            aes(x=final_lineage, y=CalcRate, colour=final_lineage))+
    geom_boxplot(na.shape=NA)+
    geom_jitter(height=0, width=0.1)+
    scale_colour_manual(values=c("darkblue","cyan3"))+
    ylab(expression(Calcification~(g~cm^-2~y^-1)))+
    xlab("")+
    ylim(c(0,2.4))+
    #geom_segment(aes(x = 1, y = 2.25, xend = 4, yend = 2.25), colour="black")+
    #geom_segment(aes(x = 2, y = 1.9, xend = 4, yend = 1.9), colour="black")+               
    #annotate(geom="text", x=3, y=1.95, label="*", size=10)+
    #annotate(geom="text", x=2.5, y=2.3, label="*", size=10)+
    theme_bw()+
    guides(colour=FALSE)+
    theme(axis.text.y = element_text(face="bold"),
          axis.text.x=element_text(face="bold"),
          axis.title.y = element_text(face="bold"))+
    facet_wrap(~Region, labeller = as_labeller(c(`OR`= "Outer Reefs",`RI` = "Rock Islands")))

t.test(subset(metadata, final_lineage=="DB" & Region=="OR")$ExtRate,subset(metadata, final_lineage=="DB" & Region=="RI")$ExtRate)

# data:  subset(metadata, final_lineage == "DB" & Region == "OR")$ExtRate and subset(metadata, final_lineage == "DB" & Region == "RI")$ExtRate
# t = 2.7921, df = 12.908, p-value = 0.01534
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  0.06643278 0.52223388
# sample estimates:
# mean of x mean of y 
#  1.109333  0.815000 


Fig5C<-ggplot(subset(metadata_dblb,!is.na(final_lineage)), 
            aes(x=final_lineage, y=ExtRate, colour=final_lineage))+
    geom_boxplot(na.shape=NA)+
    geom_jitter(height=0, width=0.1)+
    scale_colour_manual(values=c("darkblue","cyan3"))+
    ylab(expression(Extension~Rate~(cm~y^-1)))+
    xlab("")+
    ylim(c(0,2.4))+
    #geom_segment(aes(x = 1, y = 2.25, xend = 4, yend = 2.25), colour="black")+
    #geom_segment(aes(x = 2, y = 1.9, xend = 4, yend = 1.9), colour="black")+               
    #annotate(geom="text", x=3, y=1.95, label="*", size=10)+
    #annotate(geom="text", x=2.5, y=2.3, label="*", size=10)+
    theme_bw()+
    guides(colour=FALSE)+
    theme(axis.text.y = element_text(face="bold"),
          axis.text.x=element_text(face="bold"),
          axis.title.y = element_text(face="bold"))+
    facet_wrap(~Region,labeller = as_labeller(c(`OR`= "Outer Reefs",`RI` = "Rock Islands")))


t.test(subset(metadata, final_lineage=="DB" & Region=="OR")$CalcRate,subset(metadata, final_lineage=="DB" & Region=="RI")$CalcRate)

# data:  subset(metadata, final_lineage == "DB" & Region == "OR")$CalcRate and subset(metadata, final_lineage == "DB" & Region == "RI")$CalcRate
# t = 4.5354, df = 12.101, p-value = 0.000669
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  0.2809981 0.7996686
# sample estimates:
# mean of x mean of y 
# 1.3986667 0.8583333

subset(metadata, !is.na(final_lineage))%>%group_by(final_lineage, Region)%>%
  summarise(Freq98SB=sum(SB_98,na.rm=TRUE)/length(na.omit(SB_98)),
            Freq10SB=sum(SB_10,na.rm=TRUE)/length(na.omit(SB_10)),
            Count98SB=sum(SB_98,na.rm=TRUE), 
            Count10SB=sum(SB_10,na.rm=TRUE),
            TotalN_98=length(na.omit(SB_98)),
            TotalN_10=length(na.omit(SB_10)))%>%
  mutate(Count98NoSB=TotalN_98-Count98SB,
         Count10NoSB=TotalN_10-Count10SB)->SB_data_region

pivot_longer(SB_data_region, cols = c(Count98SB, Count10SB,Count98NoSB,Count10NoSB), names_to = "Band_type", values_to = "Band_counts")->SB_counts_region

SB_counts_region%>%mutate(year = case_when(c(Band_type=="Count98SB" |Band_type=="Count98NoSB") ~ "1998", TRUE ~ "2010"))->SB_counts_region

subset(SB_counts_region, final_lineage=="DB" | final_lineage=="LB")->SB_counts_region
droplevels(SB_counts_region)->SB_counts_region

sb_labs<-data.frame(Region=c("OR","OR","RI", "RI"), 
                    final_lineage=c("DB","LB","DB", "LB"),
                    labels=c("N=19", "N=9", "N=6","N=9"))

Fig5D<-ggplot(subset(SB_counts_region, year==1998 & Band_type=="Count98SB"), aes(x=final_lineage, y=Freq98SB, fill=final_lineage, colour=final_lineage))+
    geom_bar(stat="identity", aes(alpha=0.9))+
    scale_fill_manual(values=c("darkblue","cyan3","lightcoral", "red2"))+
    scale_colour_manual(values=c("darkblue","cyan3","lightcoral", "red2"))+
    theme_minimal()+ylab("Stress Band Prevalence")+
    xlab("")+facet_wrap(~Region, labeller = as_labeller(c(`OR`= "Outer Reefs",`RI` = "Rock Islands")))+
    theme_bw()+theme(panel.grid.minor=element_blank(), 
                     panel.grid.major.x = element_blank(),
                     axis.text = element_text(face="bold"),
                     axis.text.x = element_text(face="bold"), 
                     axis.ticks.x = element_blank())+
    guides(colour="none", fill="none", alpha="none")+
    scale_y_continuous(limits= c(0,0.75), breaks=c(0,0.25,0.5,0.75), labels=c("0%", "25%", "50%", "75%"))+
    geom_text(data=sb_labs, 
              mapping = aes(x = final_lineage, y = 0.05, label = labels), 
              colour="black", fontface="italic")

ggarrange(Fig5A,Fig5B,Fig5C,Fig5D, nrow=2, ncol=2, labels = c("A", "B", "C", "D"))
ggsave("Fig5.png", width=8, height=6, units="in")




## Stress band data for 2010 by lineage (Figure S4) ####
prop.test(SB_data$Count10SB, SB_data$TotalN_10)
# 4-sample test for equality of proportions without continuity correction
# 
# data:  SB_data$Count10SB out of SB_data$TotalN_10
# X-squared = 1.464, df = 3, p-value = 0.6906
# alternative hypothesis: two.sided
# sample estimates:
#   prop 1    prop 2    prop 3    prop 4 
# 0.1875000 0.1739130 0.0000000 0.2222222 

ggplot(subset(SB_counts, year==2010 & Band_type=="Count10SB"), aes(x=final_lineage, y=Freq10SB, fill=final_lineage, colour=final_lineage))+
  geom_bar(stat="identity", aes(alpha=0.9))+
  scale_fill_manual(values=c("darkblue","powderblue","lightcoral", "red2"))+
  scale_colour_manual(values=c("darkblue","powderblue","lightcoral", "red2"))+
  theme_minimal()+ylab("Stress Band Prevalence")+
  xlab("")+
  theme_minimal()+theme(panel.grid.minor=element_blank(), 
                        panel.grid.major.x = element_blank(),
                        axis.text = element_text(face="bold"),
                        axis.text.x = element_text(face="bold"), 
                        axis.ticks.x = element_blank(), 
                        axis.title = element_text(face="bold"))+
  guides(colour=FALSE, fill=FALSE, alpha=FALSE)+
  scale_y_continuous(limits= c(0,0.25), breaks=c(0,0.125,0.25), labels=c("0%","12.5%", "25%"))+
  annotate(geom="text", label="N=32", x=1, y=0.015, fontface="italic")+
  annotate(geom="text", label="N=23", x=2, y=0.015, fontface="italic")+
  annotate(geom="text", label="N=6", x=3, y=0.015, fontface="italic")+
  annotate(geom="text", label="N=9", x=4, y=0.015, fontface="italic")

ggsave("FigS4.png", width=5, height=3, units="in")

