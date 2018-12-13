library(phyloseq) #1.22.3
library(ggplot2) #2.2.1
library(ggrepel) #0.7.0
library(reshape2) #1.4.3
library(vegan) #2.4.6
library(plyr) #1.8.4
library(dplyr) #0.7.4
library(lme4) #1.1.15
library(nlme) #3.1.131
library(tibble) #1.4.2
library(cowplot) #0.9.2

otufile.t   <- "p303_run34_map97_ab5_utax.tab"
mapfile.t       <- "MAPFILE.txt"
treefile.t    <- "p303_run34_OTU_ab5.tre"
refseqfile.t  <- "p303_run34_OTU_ab5.fa"
combdata <- import_qiime(otufilename = otufile.t, mapfilename = mapfile.t, treefilename = treefile.t, refseqfilename = refseqfile.t)

#merge samples from two runs
mdata<-merge_samples(combdata, "Sample")
#fix factors in merged samples
sample_data(mdata)$sampletype <- levels(sample_data(combdata)$sampletype)[get_variable(mdata, "sampletype")]
sample_data(mdata)$controltype <- levels(sample_data(combdata)$controltype)[get_variable(mdata, "controltype")]
sample_data(mdata)$cl <- levels(sample_data(combdata)$cl)[get_variable(mdata, "cl")]
sample_data(mdata)$clone <- levels(sample_data(combdata)$clone)[get_variable(mdata, "clone")]
sample_data(mdata)$treatment <- levels(sample_data(combdata)$treatment)[get_variable(mdata, "treatment")]
sample_data(mdata)$mother <- levels(sample_data(combdata)$mother)[get_variable(mdata, "mother")]
sample_data(mdata)$umother <- levels(sample_data(combdata)$umother)[get_variable(mdata, "umother")]
sample_data(mdata)$extraction_day <- levels(sample_data(combdata)$extraction_day)[get_variable(mdata, "extraction_day")]
sample_data(mdata)$pcr_batch <- levels(sample_data(combdata)$pcr_batch)[get_variable(mdata, "pcr_batch")]
sample_data(mdata)$sedtype <- levels(sample_data(combdata)$sedtype)[get_variable(mdata, "sedtype")]
sample_data(mdata)$seddate <- levels(sample_data(combdata)$seddate)[get_variable(mdata, "seddate")]
sample_data(mdata)$X.SampleID <- levels(sample_data(combdata)$seddate)[get_variable(mdata, "X.SampleID")]
sample_data(mdata)$sampletype <- as.factor(sample_data(mdata)$sampletype)
sample_data(mdata)$controltype <- as.factor(sample_data(mdata)$controltype)
sample_data(mdata)$clone <- as.factor(sample_data(mdata)$clone)
sample_data(mdata)$treatment <- as.factor(sample_data(mdata)$treatment)
sample_data(mdata)$birthday <- as.factor(sample_data(mdata)$birthday)
sample_data(mdata)$mother <- as.factor(sample_data(mdata)$mother)
sample_data(mdata)$box <- as.factor(sample_data(mdata)$box)
sample_data(mdata)$sediment_bottle <- as.factor(sample_data(mdata)$sediment_bottle)
sample_data(mdata)$extraction_day <- as.factor(sample_data(mdata)$extraction_day)
sample_data(mdata)$extraction_round <- as.factor(sample_data(mdata)$extraction_round)
sample_data(mdata)$pcr_batch <- as.factor(sample_data(mdata)$pcr_batch)
sample_data(mdata)$sedtype <- as.factor(sample_data(mdata)$sedtype)
sample_data(mdata)$seddate <- as.factor(sample_data(mdata)$seddate)
sample_data(mdata)$eggs <- as.factor(sample_data(mdata)$eggs)

data_all<-mdata
sampdata_all<-cbind(sample_data(data_all), sample_sums(data_all)) #table of all sample data with numbers of reads
#only samples with over 5000 reads
data<-subset_samples(data_all, sample_sums(data_all)>5000)
data<-subset_taxa(data, taxa_sums(data)>0)

#animal data
adata<-subset_samples(data, sampletype=="animal")
adata<-subset_taxa(adata, taxa_sums(adata)>0)
adata_rar<-rarefy_even_depth(adata, sample.size=min(sample_sums(adata)), rngseed=3, replace=FALSE)
#data frame of variables and diversity data for animal samples
adf<-cbind(sample_data(adata), sample_sums(adata), estimate_richness(adata))
names(adf)[names(adf)=="sample_sums(adata)"]<-"totalreads"
adf$rrich<-estimate_richness(adata_rar, measures="Observed")$Observed #adding species richness measure based on rarefaction
adf$clone=factor(adf$clone, 
                 levels=c("TR-EG-1", "BE-OHZ-T10", 
                          "ES-DO1-1", "F2-82", "DE-K35-Mu10", 
                          "BE-WE-G59", "IXF1", "NO-V-7", 
                          "DE-KA-F28", "CZ-N2-6", "CZ-N1-1", "F2-918")) #clone factor levels arranged in order of avg browsing intensity

#sediment samples data
mud<-subset_samples(data, sampletype=="mud")
mud_df<-cbind(sample_data(mud), estimate_richness(mud, measures="Observed"))
mud_df$sedtype<-factor(mud_df$sedtype, labels=c("Autoclaved", "Natural"))
mud_rar<-rarefy_even_depth(mud, sample.size=min(sample_sums(mud)), rngseed=3, replace=FALSE)
plot_bar(subset_samples(mud_rar, seddate=="start"), fill="Class") + theme(legend.text=element_text(size=8))
#Figure S4
sediment_sample_richness<-plot(Observed~sedtype, mud_df, xlab="Sediment type", ylab="Bacterial species richness")

#relative abundances
adata_rel<-transform_sample_counts(adata, function(x) x/sum(x))
abundances<-cbind(adf, otu_table(adata_rel)) #table with relative abundances of each OTU for each sample
#OTU 1 and 2 in QTL panel vs diversity panel
mean(abundances$OTU_1)
sd(abundances$OTU_1)/sqrt(length(abundances$OTU_1))
abundances_qtl<-abundances[(abundances$clone=="IXF1")|(abundances$clone=="F2-82")|(abundances$clone=="F2-918"),]
mean(abundances_qtl$OTU_2)
sd(abundances_qtl$OTU_2)/sqrt(length(abundances_qtl$OTU_2))
abundances_diversitypanel<-abundances[!((abundances$clone=="IXF1")|(abundances$clone=="F2-82")|(abundances$clone=="F2-918")),]
mean(abundances_diversitypanel$OTU_2)
sd(abundances_diversitypanel$OTU_2)/sqrt(length(abundances_diversitypanel$OTU_2))

#examining frequency distribution of taxa across samples
adata_otus<-data.frame(t(otu_table(adata)))
adata_otus_pa<-decostand(adata_otus, method="pa") #presence-absence table
adata_otus_nsamples<-rowSums(adata_otus_pa)
adata_otus_freq<-adata_otus_nsamples/nsamples(adata)
hist(adata_otus_freq) #vast majority of OTUs are found in <10% of samples
otus_maj<-adata_otus_freq[adata_otus_freq>.3]
otus_min<-adata_otus_freq[adata_otus_freq<.1]

#examining alpha diversity
#individual jar-mate behavior vs. Shannon index (Fig S6)
diversity1<-ggplot(adf, aes(x=behavior, y=Shannon, color=treatment)) + geom_point() + scale_color_manual(values=c("gray","#E69F00","#56B4E9")) + stat_smooth(method=lm)
#summarize by clone and treatment
summary<-ddply(adf, .(clone, clonerank, cloneavgbehavior, treatment), summarise, 
               n=length(clone), trtbeh=mean(behavior, na.rm=TRUE), 
               InvSimp=mean(InvSimpson), Shan=mean(Shannon), rich=mean(rrich), 
               trtbeh_se=sd(behavior, na.rm=TRUE)/sqrt(n), 
               simp_se=sd(InvSimpson)/sqrt(n), shan_se=sd(Shannon)/sqrt(n), 
               rich_se=sd(rrich)/sqrt(n), median_shan=median(Shannon))
#order clones by increasing average behavior
summary$clone<-factor(summary$clone, levels=c("TR-EG-1", "BE-OHZ-T10", "ES-DO1-1", "F2-82", "DE-K35-Mu10", "BE-WE-G59", "IXF1", "NO-V-7", "DE-KA-F28", "CZ-N2-6", "CZ-N1-1", "F2-918"))

#alpha diversity measurements by clone and treatment - boxplots
#Figure 3
clone_diversity_plot_shan_box<-ggplot(adf, aes(x=clone, y=Shannon, fill=treatment))+
  theme_bw() + scale_y_continuous(limits=c(0,3.5)) + scale_fill_manual(values=c("gray","#E69F00","#56B4E9")) +
  geom_point(position=position_jitterdodge(), shape=21, alpha=0.5)+
  geom_boxplot(position=position_dodge(), outlier.size=0) +
  geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5), size=0.2) +
  theme(axis.text.x=element_text(angle=90), 
        legend.position=c(0.07,0.85), legend.title=element_blank(), 
        legend.margin=margin(t = 0, unit='cm'),
        legend.text=element_text(size=8)) + 
  labs(y="Shannon Index")
treatment_diversity_plot_shan_box<-ggplot(adf, aes(x=treatment, y=Shannon, color=treatment))+
  theme_bw() + scale_y_continuous(limits=c(0,3)) + scale_color_manual(values=c("gray","#E69F00","#56B4E9")) +
  geom_boxplot(position=position_dodge(), alpha=0) +
  geom_point(position=position_jitterdodge())+
  theme(axis.text.x=element_text(angle=90), legend.position=c(0.05,0.1)) + labs(y="Shannon Index +/- s.e.m.")

#alpha diversity measurements by clone and treatment - mean +/- SEM
clone_diversity_plot_simp<-ggplot(summary, aes(x=clone, y=InvSimp, color=treatment))+geom_point(stat="identity", position=position_dodge(width=0.9)) +
  theme_bw() + scale_y_continuous(limits=c(0,9)) + scale_color_manual(values=c("gray","#E69F00","#56B4E9")) +
  geom_errorbar(aes(ymin=(InvSimp-simp_se), ymax=(InvSimp+simp_se)), position=position_dodge(width=0.9)) + geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5), size=0.2) +
  theme(axis.text.x=element_text(angle=90), axis.title.y=element_text(size=8),
        legend.position="right", 
        legend.text=element_text(size=8), legend.title=element_blank(), 
        legend.margin=margin(t = 0, unit='cm')) + 
  labs(y="Simpson Index +/- s.e.m.")
clone_diversity_plot_rrich<-ggplot(summary, aes(x=clone, y=rich, color=treatment))+geom_point(stat="identity", position=position_dodge(width=0.9)) +
  theme_bw() + scale_y_continuous(limits=c(0,25)) + scale_color_manual(values=c("gray","#E69F00","#56B4E9")) +
  geom_errorbar(aes(ymin=(rich-rich_se), ymax=(rich+rich_se)), position=position_dodge(width=0.9)) + geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5), size=0.1) +
  theme(axis.text.x=element_text(angle=90), axis.title.y=element_text(size=8),
        legend.position="right", 
        legend.text=element_text(size=8), legend.title=element_blank(), 
        legend.margin=margin(t = 0, unit='cm')) + 
  labs(y="Richness +/- s.e.m.")
clone_diversity_plot_shan<-ggplot(summary, aes(x=clone, y=Shan, color=treatment))+geom_point(stat="identity", position=position_dodge(width=0.9)) +
  theme_bw() + scale_y_continuous(limits=c(0,2.5)) + scale_color_manual(values=c("gray","#E69F00","#56B4E9")) +
  geom_errorbar(aes(ymin=(Shan-shan_se), ymax=(Shan+shan_se)), position=position_dodge(width=0.9)) + geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5), size=0.1) +
  theme(axis.text.x=element_text(angle=90), axis.title.y=element_text(size=8),
        legend.position="right", 
        legend.text=element_text(size=8), legend.title=element_blank(),
        legend.margin=margin(t = 0, unit='cm')) + 
  labs(y="Shannon Index +/- s.e.m.")
#Figure S5
clone_diversity_plot_supplement<-plot_grid(clone_diversity_plot_rrich, clone_diversity_plot_simp, clone_diversity_plot_shan, labels=c("A","B","C"), nrow=3)

#treatment, clone and size
size_clntrt<-aov(size~clone*treatment, adf)
#behavior and size
sizebeh<-ggplot(adf, aes(x=cloneavgsize, y=cloneavgbehavior, label=clone)) + geom_point(shape=20) + geom_text(aes(labels=clone), check_overlap=FALSE, angle=45, vjust="bottom", size=3) +
  labs(x="Mean body size (mm)", y="Mean browsing intensity") + theme_bw()
#Figure S2
sizebeh_ind<-ggplot(adf, aes(x=size, y=behavior)) + geom_point(shape=20)  +
  labs(x="Body size (mm)", y="Browsing intensity") + theme_bw()
sizebeh_lm<-aov(cloneavgsize~cloneavgbehavior, adf) #clonal average behavior and clonal average size are uncorrelated
sizebeh_i_lm<-lm(size~behavior, adf) #individual behavior and size also uncorrelated
anova(sizebeh_i_lm)

#Examining differences between average microbial diversity in NET and SED treatment animals
summary_ns<-summary[summary$treatment!="AUT",]
diffs<-data.frame(clone=levels(summary_ns$clone))
diffs$cloneavgbehavior<-summary_ns[summary_ns$treatment== "NET", ][match(diffs$clone, summary_ns[summary_ns$treatment== "NET", ]$clone), 3]
diffs$ShanNet <- summary_ns[summary_ns$treatment== "NET", ][match(diffs$clone, summary_ns[summary_ns$treatment== "NET", ]$clone), 8]
diffs$ShanSed <- summary_ns[summary_ns$treatment== "SED", ][match(diffs$clone, summary_ns[summary_ns$treatment== "SED", ]$clone), 8]
diffs$shandiff<-diffs$ShanSed-diffs$ShanNet

d<-ggplot(diffs, aes(x=cloneavgbehavior, y=shandiff, label=clone)) + theme_bw(base_size = 14) + geom_point() + 
  labs(x="Clonal mean browsing intensity", y="Mean SED-Mean NET Shannon index difference") 
diffs_lm<-lm(shandiff~cloneavgbehavior, diffs)

agg_plot_shan <- ggplot(summary_ns, aes(x=trtbeh, y=Shan, color=treatment), ymin=0) + 
  theme_bw()   + 
  scale_colour_manual(values=c("#E69F00","#56B4E9")) + 
  scale_x_continuous(limits=c(3.2,5.25)) + scale_y_continuous(limits=c(0,2.2)) +
  geom_point(position=position_jitter(), stat="identity") +
  geom_errorbar(aes(ymin=(Shan-shan_se),ymax=(Shan+shan_se))) + 
  geom_errorbarh(aes(xmin=(trtbeh-trtbeh_se), xmax=(trtbeh+trtbeh_se))) + 
  geom_label_repel(aes(label=clone)) + theme_bw(base_size=14) +
  labs(x="Browsing intensity +/- s.e.m", y="Shannon Index +/- s.e.m.") +
  theme(legend.position=c(.1,.8)) 
#without labels
agg_plot_shan1 <- ggplot(summary_ns, aes(x=trtbeh, y=Shan, color=treatment), ymin=0) + 
  theme_bw()   + 
  scale_colour_manual(values=c("#E69F00","#56B4E9")) + 
  scale_x_continuous(limits=c(3.2,5.25)) + scale_y_continuous(limits=c(0,2.2)) +
  geom_point(position=position_jitter(), stat="identity") +
  geom_errorbar(aes(ymin=(Shan-shan_se),ymax=(Shan+shan_se))) + 
  geom_errorbarh(aes(xmin=(trtbeh-trtbeh_se), xmax=(trtbeh+trtbeh_se))) + 
  labs(x="Browsing intensity +/- s.e.m", y="Shannon Index +/- s.e.m.") +
  theme(legend.position=c(.1,.8)) 

#Fig 4
shannon_behavior_figure<-plot_grid(agg_plot_shan, d, labels=c("A","B"))

#Clone as fixed effect on different alpha diversity measurements
#exclude NOV-7 due to insufficient replication
adf_noNOV<-adf[adf$clone!="NO-V-7",]
#testing for batch effects
batchfx_rrich<-aov(rrich~extraction_day, adf_noNOV) #no batch effect on richness
batchfx_shan<-aov(Shannon~extraction_day, adf_noNOV) #no batch effect on Shannon index
batchfx_simp<-aov(InvSimpson~extraction_day, adf_noNOV) #batch has effect on Simpson index
#Analysis of variance
rich_aov<-aov(rrich~clone*treatment, adf_noNOV)
shan_aov<-aov(Shannon~clone*treatment, adf_noNOV)
simp_aov<-aov(InvSimpson~clone*treatment + Error(extraction_day), adf_noNOV)
#All 3 alpha diversity measures have clone or treatment as significant factors
#Continue analyses using Shannon index

#To examine Shannon in NET and SED groups only
adf_ns<-adf[adf$treatment!="AUT",]
adf_ns1<-adf_noNOV[adf_noNOV$treatment!="AUT",]

#Clone behavior/treatment effects on Shannon index
#For weighting by sample sizes: making a factor for sample sizes of each group
get_vals<-summary
get_vals$trtcl<-paste(get_vals$clone, get_vals$treatment, sep=".")
adf_ns$trtcl<-paste(adf_ns$clone, adf_ns$treatment, sep=".")
r<-match(adf_ns$trtcl, get_vals$trtcl)
adf_ns$trtcl.n<-get_vals[r,5]
adf_ns$trtbeh<-get_vals[r, "trtbeh"]
adf$trtcl<-paste(adf$clone, adf$treatment, sep=".")
r<-match(adf$trtcl, get_vals$trtcl)
adf$trtcl.n<-get_vals[r,5]
adf$trtbeh<-get_vals[r, "trtbeh"]

#basic analysis: effect of treatment and clonal average browsing intensity on Shannon
shan_behavior_aov<-aov(Shannon~treatment*cloneavgbehavior, adf_ns) #same effect when NOV7 is excluded
shan.i<-aov(Shannon~treatment*behavior, adf_ns)
#no interaction effect when individual jar-mate behavior is used as behavior proxy

#linear mixed-effect models with treatment, behavior; clone as random effect
#not weighted
shannon_model_t <- lme(Shannon~treatment+cloneavgbehavior+cloneavgsize+treatment:cloneavgbehavior+treatment:cloneavgsize, 
                       data=adf_ns, na.action=na.omit, random = ~ 1 | clone)
shan.mod_t<-anova(shannon_model_t)
#weighted by sample size
shannon_model_s <- lme(Shannon~treatment+cloneavgbehavior+cloneavgsize+treatment:cloneavgbehavior+treatment:cloneavgsize, 
                     data=adf_ns, weights= ~1/trtcl.n, 
                     na.action=na.omit, random = ~ 1 | clone)
shan.mod_s<-anova(shannon_model_s)
#model comparison
anova(shannon_model_t, shannon_model_s) 

#using treatment-specific behavior index averages
shan_t_behavior_aov<-aov(Shannon~treatment*trtbeh, adf_ns)
shannon_t_behavior_model<- lme(Shannon~treatment+trtbeh+cloneavgsize+treatment:trtbeh+treatment:cloneavgsize, data=adf_ns, na.action=na.omit, random = ~ 1 | clone)
anova(shannon_t_behavior_model)
#still have significant effects of treatment and treatment*behavior

#multidimensional comparisons of composition
adata_hell<-transform_sample_counts(adata_rel, function(x) sqrt(x)) #Hellinger transformation
adata_hell_bray<-phyloseq::distance(adata_hell, method="bray") #Bray-Curtis distance matrix
disp<-betadisper(adata_hell_bray, adf$treatment, type="centroid")
anova(disp)
#Fig 5
boxplot(disp) #NET group is less dispersed than others; exclude it - compare AUT and SED to see if microbiotas are different in different environments
adata_autsed_hell<-subset_samples(adata_hell, treatment!="NET")
sample_data(adata_autsed_hell)$clone<-factor(sample_data(adata_autsed_hell)$clone, levels=c("TR-EG-1", "BE-OHZ-T10", "ES-DO1-1", "F2-82", "DE-K35-Mu10", "BE-WE-G59", "IXF1", "NO-V-7", "DE-KA-F28", "CZ-N2-6", "CZ-N1-1", "F2-918"))
adata_autsed_hell_bray<-phyloseq::distance(adata_autsed_hell, method="bray")
adf_as<-data.frame(sample_data(adata_autsed_hell))
adf_as$clone<-factor(adf_as$clone, levels=c("TR-EG-1", "BE-OHZ-T10", "ES-DO1-1", "F2-82", "DE-K35-Mu10", "BE-WE-G59", "IXF1", "NO-V-7", "DE-KA-F28", "CZ-N2-6", "CZ-N1-1", "F2-918"))
disp_autsed<-betadisper(adata_autsed_hell_bray, adf_as$treatment, type="centroid")
anova(disp_autsed)
boxplot(disp_autsed)

#Principal Coordinates Analysis
ordi<-ordinate(adata_autsed_hell, "PCoA", "bray")
#Fig 6
autsed_ordinate_all<-plot_ordination(adata_autsed_hell, ordi, type="samples", color="treatment") + 
  scale_color_manual(values=c("gray","#56B4E9")) + theme_bw() + theme(legend.position=c(0.9,0.8))
autsed_ordinate<-plot_ordination(adata_autsed_hell, ordi, type="samples", color="treatment") + 
  scale_color_manual(values=c("gray","#56B4E9")) + theme_bw() + facet_wrap(~clone) + theme(legend.position="none")
ordination_figure<-plot_grid(autsed_ordinate_all, autsed_ordinate, labels=c("A","B"), nrow=2)
#biplot
p<-plot_ordination(adata_autsed_hell, ordi, type="split", color="Phylum", shape="treatment")

#Adonis comparison, stratified by batch
autsed_bray_adonis<-adonis(adata_autsed_hell_bray~treatment*clone, adf_as, strata=adf_as$pcr_batch)
disp_batch<-betadisper(adata_hell_bray, adf$pcr_batch, type="centroid")
plot(disp_batch, main="Samples by batch") #spider plot of samples by batch
disp_clone<-betadisper(adata_hell_bray, adf$clone, type="centroid")
plot(disp_clone, main="Samples by clone") #spider plot of samples by clone
boxplot(disp_clone)
permutest(disp_clone)
TukeyHSD(disp_clone)

#Examining the contribution of environmental bacteria to microbiota
library("DESeq2") #1.10.1
seddata<-subset_samples(data, sampletype=="mud")
seddata<-subset_samples(seddata, seddate=="start")
seddata<-subset_taxa(seddata, taxa_sums(seddata)>0)
seddiff <- phyloseq_to_deseq2(seddata, ~sedtype)
sedparam <- DESeq(seddiff)
resparam <- results(sedparam)
alpha = 0.05
sigtab = resparam[which(resparam$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(seddata)[rownames(sigtab), ], "matrix"))
head(sigtab)
res<-sigtab[order(-sigtab$log2FoldChange),]
topsedtaxa<-rownames(res[1:10,])
#the top 10 OTUs overrepresented in natural sediment
sedtaxa.8fold<-rownames(sigtab[sigtab$log2FoldChange > 8,])
sedtaxa.5fold<-rownames(sigtab[sigtab$log2FoldChange > 5,])
sedtaxa.10fold<-rownames(sigtab[sigtab$log2FoldChange > 10,])
sedtaxa.aut<-rownames(sigtab[sigtab$log2FoldChange < -5,])

adata.sedtaxa<-prune_taxa(sedtaxa.8fold, adata)
adata.autsedtaxa<-prune_taxa(sedtaxa.aut, adata)
adata.rar.sedtaxa<-prune_taxa(sedtaxa.8fold, adata_rar)
sednum<-sample_sums(adata.sedtaxa)
autsednum<-sample_sums(adata.autsedtaxa)
adf$sedreads<-sednum
adf$sedprop<-adf$sedreads/adf$totalreads
adf$sedproprar<-sample_sums(adata.rar.sedtaxa)/min(sample_sums(adata))
adf$sedtaxrich<-estimate_richness(adata.sedtaxa, measures="Observed")$Observed
adf$sedtaxrich.rar<-estimate_richness(adata.rar.sedtaxa, measures="Observed")$Observed
adf$autsedreads<-autsednum
adf$autprop<-adf$autsedreads/adf$totalreads

max(adata_otus_nsamples[sedtaxa.8fold], na.rm=TRUE)
median(adata_otus_nsamples[sedtaxa.8fold], na.rm=TRUE) #looking at how many animal samples sediment-derived bacteria tend to be found in

#Figure S7b (a and c are the same but using seprop with 5fold or 10fold)
prop<-ggplot(adf, aes(x=behavior, y=sedprop, color=treatment)) + 
  geom_point() + scale_color_manual(values=c("gray","#E69F00","#56B4E9")) + theme_bw() +
  stat_smooth(method=lm) +
  theme(legend.position=c(0.25,0.75)) +
  labs(x="Jar-mate browsing intensity", y="Proportion sediment-derived bacteria") +
  annotate("text", x=3.2, y=.75, label="Log2FC>8")

prop_sum<-ddply(adf, .(clone, clonerank, cloneavgbehavior, treatment), summarise, 
                n=length(clone), trtbeh=mean(behavior, na.rm=TRUE), proportion=mean(sedprop),
                proportion_se=sd(sedprop)/sqrt(n))
prop_sum$clone<-factor(prop_sum$clone, levels=c("TR-EG-1", "BE-OHZ-T10", "ES-DO1-1", "F2-82", "DE-K35-Mu10", "BE-WE-G59", "IXF1", "NO-V-7", "DE-KA-F28", "CZ-N2-6", "CZ-N1-1", "F2-918"))

#Fig 7
prop_box<-ggplot(adf, aes(x=clone, y=sedprop, fill=treatment)) + 
  scale_fill_manual(values=c("gray","#E69F00","#56B4E9")) + 
  geom_boxplot(outlier.size=0) +
  geom_point(position=position_jitterdodge(), shape=21, alpha=0.5) + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=90), legend.position=c(0.1,0.8),
        legend.title=element_blank()) + 
  geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5), size=0.1) +
  labs(y="Proportion sediment-derived bacteria")

prop_bar<-ggplot(prop_sum, aes(x=clone, y=proportion, fill=treatment)) + 
  geom_bar(stat="identity", position=position_dodge()) + theme_bw() +
  scale_fill_manual(values=c("gray","#E69F00","#56B4E9")) + geom_errorbar(aes(ymin=(proportion-proportion_se), ymax=(proportion+proportion_se)), position=position_dodge(width=0.9), width=0) +
  theme(axis.text.x=element_text(angle=90), legend.position=c(0.1,0.8)) + 
  labs(y="Proportion sediment-derived bacteria")


p_batch<-aov(sedprop~extraction_day, adf) #no effect of batch on prop sed-derived bacteria
prop_size<-lm(sedprop~cloneavgsize, adf) #no effect of size on prop. sed bacteria

#Examining effect of behavior on sediment-derived bacteria using different behavior proxies
#Using jar-mate behavior
sedprop_i_model <- lme(sedprop~treatment+behavior+treatment:behavior, 
                     data=adf, na.action=na.omit, 
                     random = list(~ 1 | clone) )
sedprop.i.mod<-anova(sedprop_i_model)
#using clonal average behavior
sedprop_model <- lme(sedprop~treatment+cloneavgbehavior+treatment:cloneavgbehavior, 
                     data=adf, na.action=na.omit, 
                     random = list(~ 1 | clone) )
sedprop.mod<-anova(sedprop_model)
#using treatment-specific behavior
sedprop_model_trtbeh <- lme(sedprop~treatment+trtbeh+treatment:trtbeh, 
                            data=adf, na.action=na.omit, 
                            random = list(~ 1 | clone))
sedprop.mod.tb<-anova(sedprop_model_trtbeh)

#big jump in sedprop around clonal average behavior of 4.4
low<-adf[(adf$cloneavgbehavior<4.4),]
high<-adf[adf$cloneavgbehavior>4.4,]
mean(low[low$treatment=="SED","sedprop"])
sd(low[low$treatment=="SED","sedprop"])/sqrt(length(low[low$treatment=="SED","sedprop"]))
mean(high[high$treatment=="SED","sedprop"])
sd(high[high$treatment=="SED","sedprop"])/sqrt(length(high[high$treatment=="SED","sedprop"]))
t.test(low[low$treatment=="SED","sedprop"],high[high$treatment=="SED","sedprop"])
mean(low[low$treatment!="SED","sedprop"])
mean(high[high$treatment!="SED","sedprop"])

#does sharp increase in proportion of sediment-derived bacteria 
#as a function of behavior hold if we look at individual jar-mate behavior?
adf_a<-adf[!is.na(adf$behavior),]
low_i<-adf_a[adf_a$behavior<4.4,]
high_i<-adf_a[adf_a$behavior>4.4,]
mean(low_i[low_i$treatment=="SED","sedprop"])
sd(low_i[low_i$treatment=="SED","sedprop"])/sqrt(length(low_i[low_i$treatment=="SED","sedprop"]))
mean(high_i[high_i$treatment=="SED","sedprop"])
sd(high_i[high_i$treatment=="SED","sedprop"])/sqrt(length(high_i[high_i$treatment=="SED","sedprop"]))
t.test(low_i[low_i$treatment=="SED","sedprop"], high_i[high_i$treatment=="SED","sedprop"])
mean(low_i[low_i$treatment!="SED","sedprop"])
mean(high_i[high_i$treatment!="SED","sedprop"])
t.test(low_i[low_i$treatment!="SED","sedprop"], high_i[high_i$treatment!="SED","sedprop"])
#yes

prop_i<-ggplot(adf, aes(x=behavior, y=sedprop, color=treatment)) + 
  geom_point() + scale_color_manual(values=c("gray","#E69F00","#56B4E9")) + theme_bw() +
  stat_smooth(method=lm) +
  theme(legend.position=c(0.25,0.75)) +
  labs(x="Jar-mate browsing intensity", y="Proportion sediment-derived bacteria") +
  annotate("text", x=3.2, y=.75, label="Log2FC>8") #+ facet_wrap(~pcr_batch)

#metacoder for creating taxonomic heat trees
devtools::install_github("ropensci/taxa")
devtools::install_github("grunwaldlab/metacoder")
library(metacoder) #0.3.0.1 # https://github.com/grunwaldlab/metacoder
obj<-parse_phyloseq(adata)

obj$data$otu_props <- calc_obs_props(obj, "otu_table")
obj$data$tax_abund <- calc_taxon_abund(obj, "otu_props")
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = sample_data(adata)$treatment)
obj$data$diff_table <- compare_groups(obj, dataset = "tax_abund", cols=obj$sam_data$sample_ids, groups = sample_data(adata)$treatment)
obj$data$medians<-calc_group_median(obj, dataset="tax_abund", groups=sample_data(adata)$treatment)

#Show taxonomic trees with presence/absence and median relative abundance 
#of taxa in each treatment
#Figures S3A-C
ht_AUT<-heat_tree(obj,
                  node_size = obj$data$tax_occ$AUT,
                  node_color = obj$data$medians$AUT,
                  node_label = taxon_names,
                  node_size_axis_label = "Samples with reads",
                  node_color_axis_label = "Median proportion")

ht_NET<-heat_tree(obj,
                  node_size = obj$data$tax_occ$NET,
                  node_color = obj$data$medians$NET,
                  node_label = taxon_names,
                  node_size_axis_label = "Samples with reads",
                  node_color_axis_label = "Median proportion")
ht_SED<-heat_tree(obj,
                  node_size = obj$data$tax_occ$SED,
                  node_color = obj$data$medians$SED,
                  node_label = taxon_names,
                  node_size_axis_label = "Samples with reads",
                  node_color_axis_label = "Median proportion")

ttrees<-plot_grid(ht_AUT, ht_NET, ht_SED, labels=c("AUT","NET", "SED"), nrow=3)

#Speculating: do SED-specific bacterial taxa outcompete "native" microbiota?
#Which taxa are different between AUT and SED animals?
adata_autsed<-subset_samples(adata, treatment!="NET")
adata_autsed<-subset_taxa(adata_autsed, taxa_sums(adata_autsed)>0)
animaldiff <- phyloseq_to_deseq2(adata_autsed, ~treatment)
sedparam1 <- DESeq(animaldiff)
resparam1 <- results(sedparam1)
alpha = 0.05
sigtab1 = resparam1[which(resparam1$padj < alpha), ]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(adata_autsed)[rownames(sigtab1), ], "matrix"))
head(sigtab1)
write.table(sigtab1, "autsedanimaltaxa.txt", sep=",") #a table of all the significantly differentially expressed OTUs

anres<-sigtab1[order(sigtab1$log2FoldChange),]
auttaxa<-rownames(anres[anres$log2FoldChange< -5,])
autt<-prune_taxa(auttaxa, adata_rar)
adf$au<-sample_sums(autt)
adf$auprop<-adf$au/adf$totalreads
adf$aurich<-estimate_richness(autt)$Observed
#in SED treatment group, is diversity of AUT-specific taxa 
#negatively correlated with higher relative abundance 
#of sediment-derived bacteria? 
sed.au.rich<-ggplot(adf[adf$treatment=="SED",], aes(x=sedprop, y=aurich)) + 
  geom_point(shape=20) + 
  labs(x="Relative abundance of sediment-derived bacteria", y="Number of AUT-specific taxa present") +
  annotate("text", x=.45, y=6.5, label="SED animals only")

#...sort of

#examining relationship between top 
#sediment-specific OTU and richness of AUT animal-specific OTUs
adf$o40<-abundances$OTU_40
excl<-ggplot(adf[adf$treatment=="SED",], aes(x=o40, y=aurich)) + geom_point(shape=1, position=position_jitter(width=0, height=0.15)) +
  theme_bw() + labs(x="Relative abundance of OTU_40", y="Number of AUT-specific taxa present") 

ggplot(adf, aes(x=sedprop, y=aurich)) + geom_point() + stat_smooth(method=lm)

