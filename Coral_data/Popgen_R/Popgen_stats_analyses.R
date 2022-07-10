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


setwd("/Users/hannyrivera/Documents/MIT/Research_Papers/Palau_Microsat_Paper/Drafts/Comms_Bio_submission/Palau_porites/Coral_data/Popgen_R")

## Importing rad data ####
# See snp_filtering.sh for step taken to produce this final vcf
rad_vcf<-read.vcfR("input_files/Allsamps_maf05_final_snps.vcf")

# Convert vcf to genind and genlight objects for using adegenet package
rad_gl<-vcfR2genlight(rad_vcf, n.cores=5)
rad_gin<-vcfR2genind(rad_vcf)

# Import population strata data frame
# Read in master strata file
metadata<-read.table("input_files/Strata_all_samples.txt", header=TRUE)

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
rad_K4<-read.table("input_files/rad_K4_results.txt", header=TRUE)

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
save.image(file="/Users/hannyrivera/Documents/MIT/Research_Papers/Palau_Microsat_Paper/Drafts/Comms_Bio_submission/Palau_porites/Coral_data/Popgen_R/Setup_genetic_data.RData")


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



#make long dataframe for plotting
pivot_longer(rad_K4, cols=c("LB", "PI", "DB", "RD"), names_to = "Lineage", values_to = "Assignment")->rad_K4_meta_long
left_join(rad_K4_meta_long, select(metadata, c("Sample", "Site")), by="Sample")->rad_K4_meta_long

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
rad_fst<-stamppFst(rad_gl, nboots=100, percent = 95, nclusters=3)
# plot as heat map
rad_fst_pl<-data.frame("L1"=c(rep("DB", 3), rep ("LB",2), "PI"), 
                       "L2"= c("LB", "PI", "RD", "PI", "RD", "RD"),
                       "Fst"=c(0.24,0.67,0.40,0.64,0.35, 0.70))

ggplot(rad_fst_pl, aes(x=L1, y=L2, label=as.character(Fst)))+
  geom_tile(mapping=aes(fill=Fst))+
  scale_fill_gradient(low="thistle1", high="plum3")+
  geom_text(size=3, fontface="bold")+
  guides(fill="none")+
  theme_bw()+
  theme(axis.text= element_text(face="bold"),
        axis.title=element_blank(),
        panel.grid = element_blank())

ggsave("Fig3D.png", width=3, height=3, units="in", dpi=300)



## Plot STRUCTURE assignments for msat data (Figure S1A) ####
msat_K4_data_long<-left_join(msat_K4, select(metadata, Sample, Site, Region), by="Sample")
colnames(msat_K4_data_long)[2:5]<-c("DB_msat","LB_msat","RD_msat","PI_msat")
pivot_longer(msat_K4_data_long, cols=c("PI_msat", "RD_msat","DB_msat","LB_msat"), names_to = "msat_group", values_to = "Msat_val")->msat_K4_data_long

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
  guides(fill="none")+
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

(FigS1B<-ggplot(msat_pca_df, aes(x=Axis1, y=Axis2,colour=K4_msatpop))+
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
  guides(fill="none")+
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
(FigS2A<-ggplot(struct_comp_long, aes(x=Sample, y=Assignment, fill=Lineage))+
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
        panel.grid = element_blank()))

(FigS2B<-ggplot(metadata, aes(x=K4_radpop, y=K4_msatpop, colour=K4_radpop))+
  scale_colour_manual(values=c("darkblue","cyan3","lightcoral","red2"), na.value="grey80")+
  geom_jitter(width=0.1, height=0.1)+
  guides(colour="none")+xlab("RAD Population")+
  ylab("Msat Population")+theme_bw()+
  theme(axis.title = element_text(face="bold"),
        axis.text = element_text(face="bold")))

ggarrange(FigS2A, FigS2B, ncol=1, nrow=2, labels = c("A", "B"))

## Core data by lineage (Figure 4) ####
cores<-read.table("input_files/core_data_correct_coreP13.txt", header=TRUE)
cores<-left_join(cores, select(metadata, Sample, Region, final_lineage), by="Sample")
cores<-subset(cores, !is.na(final_lineage))

# Using only growth data from 1994-2013 -minus 1998 and 2010 
cores_sub<-subset(cores, Year %in% c(seq(1994,1997),seq(1999,2009),seq(2011,2013)))
cores_sub %>% group_by(Sample, final_lineage, Region) %>% summarize(ExtRate=mean(ExtRate), 
                                                                    Density=mean(Density), 
                                                                    CalcRate=mean(CalcRate)) -> cores_run


#normality test
shapiro.test(cores_run$Density)  # okay 
shapiro.test(cores_run$CalcRate) # okay 
shapiro.test(log(cores_run$ExtRate))  # okay 

# Skeletal Density
(dens_tuk<-TukeyHSD(x=aov(cores_run$Density~cores_run$final_lineage), 'cores_run$final_lineage', conf.level=0.95))
#                 Df Sum Sq Mean Sq F value   Pr(>F)    
# final_lineage  3  0.607 0.20233   10.73 8.09e-06 ***
# Residuals     65  1.226 0.01886     


# $`metadata$final_lineage`
# diff         lwr          upr     p adj
# LB-DB -0.1003861 -0.20356310  0.002790979 0.0594065  close 
# PI-DB  0.0358956 -0.12408941  0.195880609 0.9342721
# RD-DB -0.2694528 -0.40477467 -0.134130888 0.0000105  *** 
# PI-LB  0.1362817 -0.03327489  0.305838216 0.1578327
# RD-LB -0.1690667 -0.31558033 -0.022553112 0.0173663  *** 
# RD-PI -0.3053484 -0.49617799 -0.114518777 0.0004423  ***

(Fig4A<-ggplot(cores_run, 
           aes(x=final_lineage, y=Density, colour=final_lineage))+
    geom_boxplot(na.shape=NA)+
    geom_jitter(height=0, width=0.1)+
    scale_colour_manual(values=c("darkblue","cyan3","lightcoral", "red2"))+
    ylab(expression(Skeletal~Density~(g~cm^-3)))+
    xlab("")+
    annotate(geom="text", x=4, y=1.4, label="*", size=10)+
    theme_bw()+
    guides(colour="none")+
    theme(panel.grid=element_blank(),
          axis.text.y = element_text(face="bold"),
          axis.title.y = element_text(face="bold"),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank()))

# Calcification
(calc_tuk<-TukeyHSD(x=aov(cores_run$CalcRate~cores_run$final_lineage), 'cores_run$final_lineage', conf.level=0.95))
# Df Sum Sq Mean Sq F value  Pr(>F)   
#  cores_run$final_lineage  3  2.133  0.7111   5.336 0.00239 **
#   Residuals               65  8.661  0.1333   

# $`metadata$final_lineage`
# diff        lwr         upr     p adj
# LB-DB -0.09232838 -0.3666089  0.18195210 0.8112940
# PI-DB -0.07261920 -0.4979150  0.35267663 0.9693546
# RD-DB -0.54449654 -0.9042292 -0.18476384 0.0009556  *** 
# PI-LB  0.01970918 -0.4310312  0.47044952 0.9994461
# RD-LB -0.45216816 -0.8416523 -0.06268399 0.0165015  *** 
# RD-PI -0.47187734 -0.9791689  0.03541418 0.0773336  close 


(Fig4B<-ggplot(cores_run, 
          aes(x=final_lineage, y=CalcRate, colour=final_lineage))+
  geom_boxplot(na.shape=NA)+
  geom_jitter(height=0, width=0.1)+
  scale_colour_manual(values=c("darkblue","cyan3","lightcoral", "red2"))+
  ylab(expression(Calcification~(g~cm^-2~y^-1)))+
  xlab("")+
  scale_y_continuous(limits=c(0.4,2.35), breaks = c(0.4,0.7,1.0,1.3,1.6,1.9,2.2))+
  geom_segment(aes(x = 1, y = 2.25, xend = 4, yend = 2.25), colour="black")+
  geom_segment(aes(x = 2, y = 1.9, xend = 4, yend = 1.9), colour="black")+               
  annotate(geom="text", x=3, y=1.95, label="*", size=10)+
  annotate(geom="text", x=2.5, y=2.3, label="*", size=10)+
  theme_bw()+
  guides(colour="none")+
  theme(panel.grid=element_blank(),
        axis.text.y = element_text(face="bold"),
        axis.text.x=element_blank(),
        axis.title.y = element_text(face="bold"),
        axis.ticks.x = element_blank()))

(ext_tuk<-TukeyHSD(x=aov(log(cores_run$ExtRate)~cores_run$final_lineage), 'cores_run$final_lineage', conf.level=0.95))
# $`cores_run$final_lineage
# diff        lwr         upr     p adj
# LB-DB  0.02139733 -0.2076724  0.25046710 0.9946991
# PI-DB -0.05040897 -0.4056016  0.30478371 0.9819734
# RD-DB -0.31884791 -0.6192845 -0.01841135 0.0333217  ***
# PI-LB -0.07180630 -0.4482494  0.30463677 0.9580917
# RD-LB -0.34024524 -0.6655292 -0.01496125 0.0369434  ***
# RD-PI -0.26843894 -0.6921116  0.15523375 0.3473754


(Fig4C<-ggplot(cores_run,aes(x=final_lineage, y=ExtRate, colour=final_lineage))+
  geom_boxplot()+
  geom_jitter(height=0, width=0.1)+
  scale_colour_manual(values=c("darkblue","cyan3","lightcoral", "red2"))+
  ylab(expression(Extension~Rate~(cm~y^-1)))+
  xlab("")+
  scale_y_continuous(limits=c(0.4,2.3), breaks = c(0.4,0.7,1.0,1.3,1.6,1.9,2.2))+
  theme_bw()+
  geom_segment(aes(x = 1, y = 2.15, xend = 4, yend = 2.15), colour="black")+
  geom_segment(aes(x = 2, y = 1.9, xend = 4, yend = 1.9), colour="black")+               
  annotate(geom="text", x=3, y=1.95, label="*", size=10)+
  annotate(geom="text", x=2.5, y=2.2, label="*", size=10)+
  guides(colour="none")+
  theme(panel.grid=element_blank(),
        axis.text.y = element_text(face="bold"),
        axis.text.x=element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        axis.ticks.x = element_blank()))

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

(Fig4D<-ggplot(subset(SB_counts, year==1998 & Band_type=="Count98SB"), aes(x=final_lineage, y=Freq98SB, fill=final_lineage, colour=final_lineage))+
    geom_bar(stat="identity", aes(alpha=0.9))+
    scale_fill_manual(values=c("darkblue","cyan3","lightcoral", "red2"))+
    scale_colour_manual(values=c("darkblue","cyan3","lightcoral", "red2"))+
    theme_minimal()+ylab("1998 Stress Band Prevalence")+xlab("")+
    theme_bw()+theme(panel.grid.minor=element_blank(), 
                     panel.grid.major.x = element_blank(),
                     axis.text = element_text(face="bold"),
                     axis.text.x = element_text(face="bold"), 
                     axis.ticks.x = element_blank())+
    guides(colour="none", fill="none", alpha="none")+
    scale_y_continuous(limits= c(0,0.75), breaks=c(0,0.25,0.5,0.75), labels=c("0%", "25%", "50%", "75%"))+
    annotate(geom="text", label="*", x=1, y=0.73, size=10, fontface="bold")+
    annotate(geom="text", label="N=25", x=1, y=0.025, fontface="italic")+
    annotate(geom="text", label="N=18", x=2, y=0.025, fontface="italic")+
    annotate(geom="text", label="N=6", x=3, y=0.025, fontface="italic")+
    annotate(geom="text", label="N=8", x=4, y=0.025, fontface="italic"))

ggarrange(Fig4A,Fig4B,Fig4C,Fig4D, nrow=2, ncol=2, labels = c("A", "B", "C", "D"))
ggsave("Figure4.png", width=6, height=6, units="in")

## Core data by for DB and LB lineages by Region (Figure 5) ####
subset(cores_run, final_lineage=="DB" | final_lineage =="LB")->cores_dblb
droplevels(cores_dblb)->cores_dblb
cores_dblb$reg_line<-paste(cores_dblb$final_lineage, cores_dblb$Region, sep="_")
cores_dblb$reg_line<-as.factor(cores_dblb$reg_line)

## Stat tests for differences 
shapiro.test(cores_dblb$Density)  # okay 
shapiro.test(cores_dblb$CalcRate) # okay 
shapiro.test(log(cores_dblb$ExtRate))  # okay

# Power analysis: 
cohens_d(Density~Region, data=subset(cores_dblb, final_lineage=="DB")) #2.06
cohens_d(Density~Region, data=subset(cores_dblb, final_lineage=="LB")) #0.29

cohens_d(CalcRate~Region, data=subset(cores_dblb, final_lineage=="DB")) #1.41
cohens_d(CalcRate~Region, data=subset(cores_dblb, final_lineage=="LB")) #1.13

cohens_d(ExtRate~Region, data=subset(cores_dblb, final_lineage=="DB")) #0.89
cohens_d(ExtRate~Region, data=subset(cores_dblb, final_lineage=="LB")) #1.12


#Using those effect sizes to calculate power for our post hoc pairwise comparisons
pwr.t2n.test(n1=29,n2=6,d=2.06,sig.level=0.05,power=NULL) #DB_RI v DB_OR density = 0.99
pwr.t2n.test(n1=29,n2=6,d=1.41,sig.level=0.05,power=NULL) #DB_RI v DB_OR calc rate = 0.86
pwr.t2n.test(n1=29,n2=6,d=0.89,sig.level=0.05,power=NULL) #DB_RI v DB_OR ext rate = 0.48

pwr.t2n.test(n1=11,n2=8,d=0.29,sig.level=0.05,power=NULL) #LB_RI v LB_OR density = 0.09
pwr.t2n.test(n1=11,n2=8,d=1.13,sig.level=0.05,power=NULL) #LB_RI v LB_OR calc rate = 0.63
pwr.t2n.test(n1=11,n2=8,d=1.12,sig.level=0.05,power=NULL) #LB_RI v LB_OR ext rate = 0.622


(dens_tuk_dblb<-TukeyHSD(x=aov(cores_dblb$Density~cores_dblb$reg_line), 'cores_dblb$reg_line', conf.level=0.95))
# Df Sum Sq Mean Sq F value  Pr(>F)    
# cores_dblb$reg_line  3 0.4218 0.14062   9.212 5.9e-05 ***
#   Residuals           50 0.7633 0.01527  

# $`cores_dblb$reg_line`
# diff         lwr          upr     p adj
# DB_RI-DB_OR -0.24186897 -0.38913262 -0.094605326 0.0003623 *** 
# LB_OR-DB_OR -0.12557910 -0.24184994 -0.009308265 0.0295989 ***
# LB_RI-DB_OR -0.16422085 -0.29534812 -0.033093580 0.0086585 ***
# LB_OR-DB_RI  0.11628987 -0.05035373  0.282933471 0.2606588
# LB_RI-DB_RI  0.07764812 -0.09968079  0.254977034 0.6522288
# LB_RI-LB_OR -0.03864175 -0.19121261  0.113929121 0.9067702

mean(subset(cores_dblb, reg_line=="DB_OR")$Density) #1.290815
mean(subset(cores_dblb, reg_line=="DB_RI")$Density) #1.048946

(calc_tuk_dblb<-TukeyHSD(x=aov(cores_dblb$CalcRate~cores_dblb$reg_line), 'cores_dblb$reg_line', conf.level=0.95))
# Df Sum Sq Mean Sq F value  Pr(>F)   
# cores_dblb$reg_line  3  1.993  0.6645   5.498 0.00242 **
#   Residuals           50  6.043  0.1209 

# diff         lwr          upr     p adj
# DB_RI-DB_OR -0.52873694 -0.94309697 -0.114376905 0.0072442  *** 
# LB_OR-DB_OR -0.04482659 -0.37198125  0.282328068 0.9832964
# LB_RI-DB_OR -0.37291481 -0.74187144 -0.003958169 0.0466957  *** 
# LB_OR-DB_RI  0.48391034  0.01502036  0.952800321 0.0406950  ***
# LB_RI-DB_RI  0.15582213 -0.34313342  0.654777678 0.8399883
# LB_RI-LB_OR -0.32808821 -0.75738133  0.101204901 0.1904759

mean(subset(cores_dblb, reg_line=="DB_OR")$CalcRate) #1.35132
mean(subset(cores_dblb, reg_line=="DB_RI")$CalcRate) #0.822583
mean(subset(cores_dblb, reg_line=="LB_OR")$CalcRate) #1.306493
mean(subset(cores_dblb, reg_line=="LB_RI")$CalcRate) #0.9784051

(ext_tuk_dblb<-TukeyHSD(x=aov((log(cores_dblb$ExtRate))~cores_dblb$reg_line), 'cores_dblb$reg_line', conf.level=0.95))
# Df Sum Sq Mean Sq F value Pr(>F)  
# cores_dblb$reg_line  3  0.807  0.2690   3.229 0.0301 *
#   Residuals           50  4.165  0.0833   
#anova is sig but not the pairwise contrasts

# $`cores_dblb$reg_line`
# diff          lwr        upr     p adj
# DB_RI-DB_OR -0.29709349 -0.641093904 0.04690693 0.1128501
# LB_OR-DB_OR  0.08827254 -0.183330254 0.35987533 0.8234117
# LB_RI-DB_OR -0.19151558 -0.497822239 0.11479109 0.3544935
# LB_OR-DB_RI  0.38536602 -0.003904993 0.77463704 0.0532720
# LB_RI-DB_RI  0.10557791 -0.308653450 0.51980927 0.9051929
# LB_RI-LB_OR -0.27978811 -0.636185932 0.07660971 0.1716124

(Fig5A<-ggplot(cores_dblb, aes(x=Region, y=Density, colour=final_lineage))+
    geom_boxplot()+
    geom_jitter(height=0, width=0.1)+
    scale_colour_manual(values=c("darkblue","cyan3"))+
    ylab(expression(Skeletal~Density~(g~cm^-3)))+
    scale_y_continuous(limits=c(0.8, 1.62))+
    xlab("")+
    theme_bw()+
    guides(colour="none")+
    theme(panel.grid=element_blank(),
          axis.text.y = element_text(face="bold"),
          axis.title.y = element_text(face="bold"),
          axis.text.x=element_text(face="bold"))+
    facet_wrap(~final_lineage))

(Fig5B<-ggplot(cores_dblb, aes(x=Region, y=CalcRate, colour=final_lineage))+
    geom_boxplot()+
    geom_jitter(height=0, width=0.1)+
    scale_colour_manual(values=c("darkblue","cyan3"))+
    ylab(expression(Calcification~(g~cm^-2~y^-1)))+
    xlab("")+
    theme_bw()+
    guides(colour="none")+
    scale_y_continuous(limits=c(0.45,2.3))+
    theme(panel.grid=element_blank(),
          axis.text.y = element_text(face="bold"),
          axis.title.y = element_text(face="bold"),
          axis.text.x=element_text(face="bold"))+
    facet_wrap(~final_lineage))

(Fig5C<-ggplot(cores_dblb,aes(x=Region, y=CalcRate, colour=final_lineage))+
    geom_boxplot()+
    geom_jitter(height=0, width=0.1)+
    scale_colour_manual(values=c("darkblue","cyan3"))+
    ylab(expression(Extension~Rate~(cm~y^-1)))+
    xlab("")+
    theme_bw()+
    guides(colour="none")+
    theme(panel.grid=element_blank(),
          axis.text.y = element_text(face="bold"),
          axis.title.y = element_text(face="bold"),
          axis.text.x=element_text(face="bold"))+
    facet_wrap(~final_lineage))


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

prop.test(subset(SB_data_region,final_lineage=="DB" | final_lineage=="LB")$Count98SB, subset(SB_data_region,final_lineage=="DB" | final_lineage=="LB")$TotalN_98)

pairwise.prop.test(subset(SB_data_region,final_lineage=="DB" | final_lineage=="LB")$Count98SB, subset(SB_data_region,final_lineage=="DB" | final_lineage=="LB")$TotalN_98)
# 1    2    3   
# 2 1.00 -    -   
# 3 0.18 1.00 -   
# 4 0.18 1.00 1.00

sb_labs<-data.frame(Region=c("OR","OR","RI", "RI"), 
                    final_lineage=c("DB","LB","DB", "LB"),
                    labels=c("N=19", "N=9", "N=6","N=9"))

(Fig5D<-ggplot(subset(SB_counts_region, year==1998 & Band_type=="Count98SB"), aes(x=Region, y=Freq98SB, fill=final_lineage, colour=final_lineage))+
    geom_bar(stat="identity", aes(alpha=0.9))+
    scale_fill_manual(values=c("darkblue","cyan3","lightcoral", "red2"))+
    scale_colour_manual(values=c("darkblue","cyan3","lightcoral", "red2"))+
    theme_minimal()+ylab("Stress Band Prevalence")+
    xlab("")+facet_wrap(~final_lineage)+
    theme_bw()+theme(panel.grid.minor=element_blank(), 
                     panel.grid.major.x = element_blank(),
                     axis.text = element_text(face="bold"),
                     axis.text.x = element_text(face="bold"), 
                     axis.ticks.x = element_blank())+
    guides(colour="none", fill="none", alpha="none")+
    scale_y_continuous(limits= c(0,0.8), breaks=c(0,0.25,0.5,0.75), labels=c("0%", "25%", "50%", "75%"))+
    geom_text(data=sb_labs,
              mapping = aes(x = Region, y = 0.05, label = labels),
              colour="black", fontface="italic"))


ggarrange(Fig5A,Fig5B,Fig5C,Fig5D, nrow=2, ncol=2, labels = c("A", "B", "C", "D"))
ggsave("Figure5.png", width=8, height=6, units="in")


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
  guides(colour="none", fill="none", alpha="none")+
  scale_y_continuous(limits= c(0,0.25), breaks=c(0,0.125,0.25), labels=c("0%","12.5%", "25%"))+
  annotate(geom="text", label="N=32", x=1, y=0.015, fontface="italic")+
  annotate(geom="text", label="N=23", x=2, y=0.015, fontface="italic")+
  annotate(geom="text", label="N=6", x=3, y=0.015, fontface="italic")+
  annotate(geom="text", label="N=9", x=4, y=0.015, fontface="italic")

ggsave("FigS4.png", width=5, height=3, units="in")


## Within lineage pop-structure across sites gen struct ####

# Separate the lineages into different objects 
# RAD data 
LB_index<-which(rad_strata$K4_radpop=="LB")
PI_index<-which(rad_strata$K4_radpop=="PI")
RD_index<-which(rad_strata$K4_radpop=="RD")
DB_index<-which(rad_strata$K4_radpop=="DB")

rad_lg_LB<-rad_gl[LB_index,]
rad_lg_DB<-rad_gl[DB_index,]
rad_lg_PI<-rad_gl[PI_index,]
rad_lg_RD<-rad_gl[RD_index,]

setPop(rad_lg_LB)<-~Site
setPop(rad_lg_DB)<-~Site
setPop(rad_lg_PI)<-~Site
setPop(rad_lg_RD)<-~Site

rad_fst_LB<-stamppFst(rad_lg_LB, nboots=100, percent = 95, nclusters=3)
rad_fst_DB<-stamppFst(rad_lg_DB, nboots=100, percent = 95, nclusters=3)
rad_fst_PI<-stamppFst(rad_lg_PI, nboots=100, percent = 95, nclusters=3)
rad_fst_RD<-stamppFst(rad_lg_RD, nboots=100, percent = 95, nclusters=3)

# PCA 
# LB 
LB_rad_PCA <- glPca(rad_lg_LB)
LB_rad_PCA_res<-as.data.frame(LB_rad_PCA$scores)
LB_rad_PCA_res$Sample<-row.names(LB_rad_PCA_res)
LB_rad_PCA_res<-left_join(LB_rad_PCA_res,select(rad_strata, Sample, Site, Region), by="Sample")

(A1<-ggplot(LB_rad_PCA_res, aes(x=PC1, y=PC2,shape=Region,  colour=Region))+
  geom_point(size=2)+stat_ellipse(geom="polygon", aes(colour=Region),alpha=0)+
  scale_shape_manual(values=c(16,17), labels=c("Outer Reefs","Rock Islands"))+
  scale_colour_manual(values=c("blue", "red2"), labels=c("Outer Reefs","Rock Islands"))+
  ggtitle("LB lineage, RAD data")+
  theme_bw()+theme(axis.title = element_text(face="bold"),
                   axis.ticks = element_blank(),
                   axis.text = element_blank(),
                   plot.title = element_text(face="bold")))  

# DB 
DB_rad_PCA <- glPca(rad_lg_DB)
DB_rad_PCA_res<-as.data.frame(DB_rad_PCA$scores)
DB_rad_PCA_res$Sample<-row.names(DB_rad_PCA_res)
DB_rad_PCA_res<-left_join(DB_rad_PCA_res,select(rad_strata, Sample, Site, Region), by="Sample")

(B1<-ggplot(DB_rad_PCA_res, aes(x=PC1, y=PC2,shape=Region,  colour=Region))+
  geom_point(size=2)+stat_ellipse(geom="polygon", aes(colour=Region),alpha=0)+
  scale_shape_manual(values=c(16,17), labels=c("Outer Reefs","Rock Islands"))+
  scale_colour_manual(values=c("blue", "red2"), labels=c("Outer Reefs","Rock Islands"))+
  ggtitle("DB lineage, RAD data")+
  theme_bw()+theme(axis.title = element_text(face="bold"),
                   axis.ticks = element_blank(),
                   axis.text = element_blank(),
                   plot.title = element_text(face="bold")))  



# PI 
toRemove <- is.na(glMean(rad_lg_PI, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove) # position of entirely non-typed loci
rad_lg_PI_clean <- rad_lg_PI[, !toRemove]

PI_rad_PCA <- glPca(rad_lg_PI_clean)
PI_rad_PCA_res<-as.data.frame(PI_rad_PCA$scores)
PI_rad_PCA_res$Sample<-row.names(PI_rad_PCA_res)
PI_rad_PCA_res<-left_join(PI_rad_PCA_res,select(rad_strata, Sample, Site, Region), by="Sample")

(C1<-ggplot(PI_rad_PCA_res, aes(x=PC1, y=PC2,shape=Region,  colour=Region))+
  geom_point(size=2)+stat_ellipse(geom="polygon", aes(colour=Region),alpha=0)+
  scale_shape_manual(values=c(16,17), labels=c("Outer Reefs","Rock Islands"))+
  scale_colour_manual(values=c("blue", "red2"), labels=c("Outer Reefs","Rock Islands"))+
  ggtitle("PI lineage, RAD data")+
  theme_bw()+theme(axis.title = element_text(face="bold"),
                   axis.ticks = element_blank(),
                   axis.text = element_blank(),
                   plot.title = element_text(face="bold")))  


# RD
toRemove <- is.na(glMean(rad_lg_RD, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove) # position of entirely non-typed loci
rad_lg_RD_clean <- rad_lg_RD[, !toRemove]

RD_rad_PCA <- glPca(rad_lg_RD_clean)
RD_rad_PCA_res<-as.data.frame(RD_rad_PCA$scores)
RD_rad_PCA_res$Sample<-row.names(RD_rad_PCA_res)
RD_rad_PCA_res<-left_join(RD_rad_PCA_res,select(rad_strata, Sample, Site, Region), by="Sample")

(D1<-ggplot(RD_rad_PCA_res, aes(x=PC1, y=PC2,shape=Region,  colour=Region))+
  geom_point(size=2)+stat_ellipse(geom="polygon", aes(colour=Region),alpha=0)+
  scale_shape_manual(values=c(16,17), labels=c("Outer Reefs","Rock Islands"), limits=c("OR", "RI"))+
  scale_colour_manual(values=c("blue", "red2"), labels=c("Outer Reefs","Rock Islands"), limits=c("OR", "RI"))+
  ggtitle("RD lineage, RAD data")+
  theme_bw()+theme(axis.title = element_text(face="bold"),
                   axis.ticks = element_blank(),
                   axis.text = element_blank(),
                   plot.title = element_text(face="bold")) )

# microsat data
pop(msat_genid)
msat_genid_sep<-seppop(msat_genid)

setPop(msat_genid_sep$LB)<-~Site
setPop(msat_genid_sep$DB)<-~Site
setPop(msat_genid_sep$PI)<-~Site
setPop(msat_genid_sep$RD)<-~Site


msat_fst_LB<-genet.dist(msat_genid_sep$LB, method="Nei87")
msat_fst_DB<-genet.dist(msat_genid_sep$DB, method="Nei87")
msat_fst_PI<-genet.dist(msat_genid_sep$PI, method="Nei87")
msat_fst_RD<-genet.dist(msat_genid_sep$RD, method="Nei87")


# PCA 
# LB 
msat_scale_LB<-scaleGen(msat_genid_sep$LB, NA.method="mean")
msat_pca_LB<-dudi.pca(msat_scale_LB, cent=FALSE, scale=FALSE, scannf = FALSE, nf = 3)
# gets the percent variance explained by the PCs
msat_pca_eig_LB<-get_eigenvalue(msat_pca_LB)
msat_pca_df_LB<-msat_pca_LB$li
msat_pca_df_LB$Sample<-rownames(msat_pca_df_LB)
msat_pca_df_LB<-left_join(msat_pca_df_LB, select(msat_strata,Sample,Site, Region), by="Sample")

(A2<-ggplot(msat_pca_df_LB, aes(x=Axis1, y=Axis2, colour=Region, shape=Region))+
           geom_point(size=2)+stat_ellipse(geom="polygon", aes(colour=Region),alpha=0)+
           scale_shape_manual(values=c(16,17), labels=c("Outer Reefs","Rock Islands"))+
           scale_colour_manual(values=c("blue", "red2"), labels=c("Outer Reefs","Rock Islands"))+
    xlab("PC1")+
    ylab("PC2")+
          theme_bw()+
           ggtitle("LB lineage, msat data")+
           theme(axis.title = element_text(face="bold"),
                 axis.ticks = element_blank(),
                 axis.text = element_blank()))


# DB 
msat_scale_DB<-scaleGen(msat_genid_sep$DB, NA.method="mean")
msat_pca_DB<-dudi.pca(msat_scale_DB, cent=FALSE, scale=FALSE, scannf = FALSE, nf = 3)
# gets the percent variance explained by the PCs
msat_pca_eig_DB<-get_eigenvalue(msat_pca_DB)
msat_pca_df_DB<-msat_pca_DB$li
msat_pca_df_DB$Sample<-rownames(msat_pca_df_DB)
msat_pca_df_DB<-left_join(msat_pca_df_DB, select(msat_strata,Sample,Site, Region), by="Sample")

(B2<-ggplot(msat_pca_df_DB, aes(x=Axis1, y=Axis2, colour=Region, shape=Region))+
  geom_point(size=2)+stat_ellipse(geom="polygon", aes(colour=Region),alpha=0)+
  scale_shape_manual(values=c(16,17), labels=c("Outer Reefs","Rock Islands"))+
  scale_colour_manual(values=c("blue", "red2"), labels=c("Outer Reefs","Rock Islands"))+
    xlab("PC1")+
    ylab("PC2")+
    theme_bw()+
  ggtitle("DB lineage, msat data")+
  theme(axis.title = element_text(face="bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank()))

# PI 
msat_scale_PI<-scaleGen(msat_genid_sep$PI, NA.method="mean")
msat_pca_PI<-dudi.pca(msat_scale_PI, cent=FALSE, scale=FALSE, scannf = FALSE, nf = 3)
# gets the percent variance explained by the PCs
msat_pca_eig_PI<-get_eigenvalue(msat_pca_PI)
msat_pca_df_PI<-msat_pca_PI$li
msat_pca_df_PI$Sample<-rownames(msat_pca_df_PI)
msat_pca_df_PI<-left_join(msat_pca_df_PI, select(msat_strata,Sample,Site, Region), by="Sample")

(C2<-ggplot(msat_pca_df_PI, aes(x=Axis1, y=Axis2, colour=Region, shape=Region))+
  geom_point(size=2)+stat_ellipse(geom="polygon", aes(colour=Region),alpha=0)+
  scale_shape_manual(values=c(16,17), labels=c("Outer Reefs","Rock Islands"))+
  scale_colour_manual(values=c("blue", "red2"), labels=c("Outer Reefs","Rock Islands"))+
    xlab("PC1")+
    ylab("PC2")+
    theme_bw()+
  ggtitle("PI lineage, msat data")+
  theme(axis.title = element_text(face="bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank()))

# RD 
msat_scale_RD<-scaleGen(msat_genid_sep$RD, NA.method="mean")
msat_pca_RD<-dudi.pca(msat_scale_RD, cent=FALSE, scale=FALSE, scannf = FALSE, nf = 3)
# gets the percent variance explained by the PCs
msat_pca_eig_RD<-get_eigenvalue(msat_pca_RD)
msat_pca_df_RD<-msat_pca_RD$li
msat_pca_df_RD$Sample<-rownames(msat_pca_df_RD)
msat_pca_df_RD<-left_join(msat_pca_df_RD, select(msat_strata,Sample,Site, Region), by="Sample")

(D2<-ggplot(msat_pca_df_RD, aes(x=Axis1, y=Axis2, colour=Region, shape=Region))+
  geom_point(size=2)+stat_ellipse(geom="polygon", aes(colour=Region),alpha=0)+
  scale_shape_manual(values=c(16,17), labels=c("Outer Reefs","Rock Islands"))+
  scale_colour_manual(values=c("blue", "red2"), labels=c("Outer Reefs","Rock Islands"))+
  xlab("PC1")+
  ylab("PC2")+
  theme_bw()+
  ggtitle("RD lineage, msat data")+
  theme(axis.title = element_text(face="bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank()))

ggarrange(A1,A2, B1, B2, C1, C2, D1, D2, ncol=2, nrow=4, common.legend = TRUE, legend = "bottom")
ggsave("within_lineage_pcas.png", width=6, height=10, units="in", dpi=300)

## Symbiont popgen ####
# This contains 218 bi-allelic SNPs from 137 individuals 
snp_gen_p<-read.vcfR("input_files/sym_snps_final.recode.vcf")

# Convert vcf to genind and genlight object for using adegenet package
snp_gen_sym<-vcfR2genlight(snp_gen_p, n.cores=5)
snp_genid_sym<-vcfR2genind(snp_gen_p)

sym_meta<-read.table("input_files/sym_meta.txt", header=TRUE)

metadata_sym<-left_join(sym_meta, metadata, by="Sample")

strata(snp_genid_sym)<-metadata_sym
strata(snp_gen_sym)<-metadata_sym

#Set populations to Site first 
setPop(snp_genid_sym)<-~final_lineage
setPop(snp_gen_sym)<-~final_lineage

snp_fst<-stamppFst(snp_gen_sym, nboots=100, percent = 95, nclusters=3)

# $Fsts
#            LB         DB        PI RD
# LB         NA         NA        NA NA
# DB 0.04467465         NA        NA NA
# PI 0.21344820 0.20770618        NA NA
# RD 0.06600322 0.07450906 0.2535492 NA


sym_fst_pl<-data.frame("L1"=c(rep("DB", 3), rep ("LB",2), "PI"), 
                       "L2"= c("LB", "PI", "RD", "PI", "RD", "RD"),
                       "Fst"=c(0.04,0.21,0.07,0.21,0.06, 0.25))

ggplot(sym_fst_pl, aes(x=L1, y=L2, label=as.character(Fst)))+
  geom_tile(mapping=aes(fill=Fst))+
  scale_fill_gradient(low="thistle1", high="plum3")+
  geom_text(size=3, fontface="bold")+
  guides(fill=FALSE)+
  theme_bw()+
  theme(axis.text= element_text(face="bold"),
        axis.title=element_blank(),
        panel.grid = element_blank())

ggsave("Sym_fst.png", width=3, height=3, units="in", dpi=300)

Snp_pca_data_sym<-scaleGen(snp_genid_sym, NA.method="mean")
Snp_pca_sym<-dudi.pca(Snp_pca_data_sym)
Snp_pca_sym<-dudi.pca(Snp_pca_data_sym, cent=FALSE, scale=FALSE,scannf = FALSE, nf = 3)
#dudi.pca(df = Snp_pca_data_sym, scannf = FALSE, nf = 6)
col<-funky(15)
s.class(Snp_pca_sym$li, )
s.class(Snp_pca_sym$li, as.factor(pop(snp_genid_sym)), col=transp(col, 0.8), axesell = FALSE, cstar = 0, clabel=0.5, cpoint=3)
s.class(Snp_pca_sym$li, pop(snp_genid_sym), col=transp(col, 0.8), axesell = FALSE, cstar = 0, xax=1, yax=3, clabel=0.5, cpoint=3)
#setPop(snp_genid_sym)<-~MsatPop

fviz_pca_var(Snp_pca_sym, repel=TRUE, select.var = list(cos2 = 0.4), label="none")

(inv.pca.sym<-fviz_pca_ind(Snp_pca_sym, geom.ind = "point",col.ind =  pop(snp_genid_sym), # color by groups
                           palette= c( "cyan3","darkblue","lightcoral", "red1"),
                           addEllipses = TRUE, # Concentration ellipses
                           ellipse.type = "confidence",
                           legend.title = "Host Lineage",
                           repel = TRUE))

inv.pca.sym+scale_shape_manual(values=c(rep(16,12)))+
  guides(fill="none", shape="none")+ggtitle("Symbiont RAD Snps PCA")

ggsave("PCA_indvs_sym.png", dpi=300, width=6, height=4, units="in")


xval_sym <- xvalDapc(Snp_pca_data_sym, pop(snp_genid_sym), n.pca.max = 300, training.set = 0.9,
                     result = "groupMean", center = TRUE, scale = FALSE,
                     n.pca = NULL, n.rep = 30, xval.plot = TRUE)
#60 PCs achieve the lowest MSE

Snp_grp_p<-find.clusters(snp_genid_sym, max.n.clust = 10)
sym_dapc<-dapc(snp_genid_sym, Snp_grp_p$grp, n.pca = 60, n.da = 3)

myCol <- c( "red1","darkblue","cyan3","lightcoral") #These are matched in order based on: 
# table(pop(snp_genid_sym), Snp_grp_p$grp) 
# Same comment here as the RAD section above. The colors may need reording based on the way the clusters order themselves

#2D DAPC (Fig S1C)
pdf("sym_dpc.pdf", width=5, height=4)
scatter(sym_dapc, scree.da=FALSE, bg="white", pch=19, clab=0,
        cstar=0, col=myCol, scree.pca=TRUE, posi.pca="topright")
dev.off()

#1D DAPC (DF1) (Fig S1C insert)
pdf("sym_dpc_1dim.pdf", width=5, height=4)
scatter(sym_dapc,1,1, col=myCol, bg="white",scree.da=FALSE, legend=FALSE, solid=.4)
dev.off()
