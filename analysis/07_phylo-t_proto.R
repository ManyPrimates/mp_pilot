### R script for phylogenetic t-testing

## Goal is to incorporate this into RMarkdown


library(phytools)
library(dplyr)
library(tidyverse)



mp.tree <- read.nexus('analysis/TreeBlock_10kTrees_Primates_Version3.nex')
# summary(mp.tree[[1]])
# head(mp.tree[[1]])

mp.tree[[1]]$tip.label


plotTree(mp.tree[[1]], type='cladogram')


### Preparing test objects

phylo.df = aggregate(mp_data$correct, by=list(mp_data$species,mp_data$condition),
                     FUN=mean)

phylo.df$sd = aggregate(mp_data$correct, by=list(mp_data$species,mp_data$condition),
                      FUN=sd)[[3]]

colnames(phylo.df) = c('species','condition','mean','sd')

phylo.df %>%  
  filter(species != 'black_faced_spider_monkey') -> phylo.df


## Reorder df to match tree

# Manuel: I noticed that the species names did not match the performance averages in the data frames belwo. I changed the code to fix that.

mp.tree[[1]]$tip.label
phylo.df$species = as.factor(phylo.df$species)

phylo.df <- phylo.df%>%
  mutate(species = case_when(species == "chimpanzee" ~ "Pan_troglodytes_verus",
                                   species == "ring_tailed_lemur" ~ "Lemur_catta",
                                   species == "bonobo"~ "Pan_paniscus",
                                   species == "gorilla"~ "Gorilla_gorilla_gorilla",
                                   species == "black_and_white_ruffed_lemur"~ "Varecia_variegata_variegata",
                                   species == "brown_capuchin_monkey" ~ "Cebus_apella",
                                   species == "long_tailed_macaque" ~ "Macaca_fascicularis",
                                   species == "squirrel_monkey"~ "Saimiri_sciureus",
                                   species == "orangutan"~ "Pongo_pygmaeus",
                                   species == "barbary_macaque"~ "Macaca_sylvanus",
                                   species == "rhesus_macaque"~ "Macaca_mulatta"),
         species = as.factor(species))


# not sure what that code does ... delete?
## DMA - this makes the row order of the testing data match 
## the order of entries in the phylogenetic tree

phylo.df$species = relevel(phylo.df$species,'Macaca_mulatta')
phylo.df$species = relevel(phylo.df$species,'Macaca_fascicularis')
phylo.df$species = relevel(phylo.df$species,'Pongo_pygmaeus')
phylo.df$species = relevel(phylo.df$species,'Pan_troglodytes_verus')
phylo.df$species = relevel(phylo.df$species,'Pan_paniscus')
phylo.df$species = relevel(phylo.df$species,'Gorilla_gorilla_gorilla')
phylo.df$species = relevel(phylo.df$species,'Saimiri_sciureus')
phylo.df$species = relevel(phylo.df$species,'Cebus_apella')
phylo.df$species = relevel(phylo.df$species,'Varecia_variegata_variegata')
phylo.df$species = relevel(phylo.df$species,'Lemur_catta')

phylo.df = phylo.df[order(phylo.df$species),]



## a one sample t-test is just a paired t-test where the values are paired with the 0 vector

phylo.table.short = data.frame(Ha=numeric(11),H0=rep.int(1/3,11),Ha.SE=numeric(11),H0.SE=rep.int(0,11))
row.names(phylo.table.short) = phylo.df%>%filter(condition == "short")%>%pull(species)

phylo.table.short$Ha = phylo.df$mean[phylo.df$condition=='short']
phylo.table.short$Ha.SE = phylo.df$sd[phylo.df$condition=='short']
phylo.table.short = as.matrix(phylo.table.short)

phylo.table.medium = phylo.table.short
phylo.table.medium[,'Ha'] = phylo.df$mean[phylo.df$condition=='medium']
phylo.table.medium[,'Ha.SE'] = phylo.df$sd[phylo.df$condition=='medium']

phylo.table.long = phylo.table.short
phylo.table.long[,'Ha'] = phylo.df$mean[phylo.df$condition=='long']
phylo.table.long[,'Ha.SE'] = phylo.df$sd[phylo.df$condition=='long']



### Estimate lambda starting value in advance

l.short = phylosig(mp.tree[[1]],x=phylo.table.short[,'Ha'],
         method='lambda')

l.medium = phylosig(mp.tree[[1]],x=phylo.table.medium[,'Ha'],
                   method='lambda')

l.long = phylosig(mp.tree[[1]],x=phylo.table.long[,'Ha'],
                   method='lambda')



### T-testing

short.t = phyl.pairedttest(mp.tree[[1]],
                           x1=phylo.table.short[,'Ha'],x2=phylo.table.short[,'H0'],
                           se1=phylo.table.short[,'Ha.SE'],se2=phylo.table.short[,'H0.SE'],
                           lambda = l.short$lambda
                           )
short.t

medium.t = phyl.pairedttest(mp.tree[[1]],
                           x1=phylo.table.medium[,'Ha'],x2=phylo.table.medium[,'H0'],
                           se1=phylo.table.medium[,'Ha.SE'],se2=phylo.table.medium[,'H0.SE'],
                           lambda = l.medium$lambda
)
medium.t



long.t = phyl.pairedttest(mp.tree[[1]],
                            x1=phylo.table.long[,'Ha'],x2=phylo.table.long[,'H0'],
                            se1=phylo.table.long[,'Ha.SE'],se2=phylo.table.long[,'H0.SE'],
                          lambda = l.long$lambda
)
long.t



### Possible LR tests

## Trying to reproduce _ et al. (2010)
# LR.D = 2.579 - 2*0 # second part of this is the question
# Ddf = 1
# 1 - pchisq(LR.D, Ddf, ncp=0, lower.tail=TRUE, log.p=FALSE) ## this matches Table 2, line 2


LR.D = abs(2*short.t$logL)
Ddf = 1
1 - pchisq(LR.D, Ddf, ncp=0, lower.tail=TRUE, log.p=FALSE)



LR.D = abs(2*medium.t$logL)
Ddf = 1
1 - pchisq(LR.D, Ddf, ncp=0, lower.tail=TRUE, log.p=FALSE)



LR.D = abs(2*long.t$logL)
Ddf = 1
1 - pchisq(LR.D, Ddf, ncp=0, lower.tail=TRUE, log.p=FALSE)
