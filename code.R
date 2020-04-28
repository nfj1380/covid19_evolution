#COde to analyse differences in covid19 genomes

#Written by Nick Fountain-Jones (Nick.FountainJones@utas.edu.au) and Xavier Didelot

library(treeio)
library(ggtree)
library(tidyverse)
library(treestructure) 
library(phylodyn)# needs to options(buildtools.check = function(action) TRUE ) before installing from Git
library(skygrowth)
library(coda)
library(ape)

#------------------------------------------------------------------------
##################Visualize MCC Tree##################
#------------------------------------------------------------------------
#import .tre file
beast <- read.beast('covid_ExpClock1_test3.treReduced3.txt')

#can read nexus too
b <-  read.nexus('BEASTtreeApril2019.tree')


#can see node numbers
ph <- ggtree(beast) + geom_tiplab(size=2)+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
ph  

#add branch lengths etc
p1 <- ggtree(beast , mrsd="2020-03-24") + theme_tree2()+ 
  #geom_tiplab(align=FALSE, linetype='dashed', linesize=0.5, size=1) +  
  #geom_range("length_0.95_HPD", color='red', size=2, alpha=.5) + #gets length estimates 
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior)> 0.8, 
                 x=branch), vjust=0)
p1

#can isolate particular sections of the tree
viewClade(ph+geom_tiplab(size=2), node=764)

#------------------------------------------------------------------------
##################Treedater method##################
#------------------------------------------------------------------------

library(treedater)
library(lubridate)
# first we need an unrooted maximum likelihood undated tree. We used PhyML ot generate this tree

tr<-unroot(read.tree('PhyMLtree'))
plot(tr)

#resolved polytomies
trNoPoly <- multi2di(tr, random = TRUE) # getting weird resuls. I guess this wasn't what you were suggesting Erik?

# I guess this could be needed  -really not sure?
#trNoPoly $edge.length <- trNoPoly $edge.length * 1e-3

plot(trNoPoly,no.mar=T, cex = .2 )

RIGHT = function(x,n){
  substring(x,nchar(x)-n+1)
}
# second we need a named vector containing the dates of the genomes in decimal years. assuming this is contained as a suffix of the tip labels, we can extract this using:

dates<-ymd(RIGHT(tr$tip.label, 10))
dateD <- decimal_date(dates)
names(dateD)<-tr$tip.label

# third we need in the alignment length
l<-29442

#test for relaxed clock
rctest <- relaxedClockTest(tr, dateD,l, nreps = 100, overrideTempConstraint = T,
                 ncpu = 1) #looks like uncorrelated clock is best supported?


# finally we can run treedater and print the result. I kept the strict clock as I was skeptical there was enough signal for a relaxed clock.
resNoPoly<-dater( trNoPoly, dateD, l)
plot(resNoPoly, no.mar=T, cex = .2 ) #weird results - substitutino rate too small?
rootToTipRegressionPlot(resNoPoly) #yuck

pb <- parboot(resresNoPoly, ncpu = 1, nreps = 100) #lineage dating way too early.
plot(pb)

resPoly <-dater( tr, dateD, l) #with polytomies
plot(resPoly, no.mar=T, cex = .2 )
rootToTipRegressionPlot(resPoly)

str(resPoly)
writeNexus(resPoly, "test1.nex")

par(mfrow=c(1,1))
goodnessOfFitPlot(resPoly) # not sure how to interpret this. 

print(resPoly)
str(resPoly)

pbwPoly <- parboot(resPoly, ncpu = 1, nreps = 100) #these results make sense. Takes 4 hours to run
plot(pbwPoly, ggplot=TRUE)

str(pbwPoly)
#plot treeDater object
#can see node numbers

#belwo isn't working yet
phTT <- ggtree(resPoly) + geom_tiplab(size=2)+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
phTT  



#add branch lengths etc
p1 <- ggtree(beast , mrsd="2020-03-24") + theme_tree2()+ 
  #geom_tiplab(align=FALSE, linetype='dashed', linesize=0.5, size=1) +  
  #geom_range("length_0.95_HPD", color='red', size=2, alpha=.5) + #gets length estimates 
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior)> 0.8, 
                 x=branch), vjust=0)
p1

#can isolate particular sections of the tree
viewClade(ph+geom_tiplab(size=2), node=764)

#------------------------------------------------------------------------
##################look for non-random tree structure using treestructure (Volz et al)##################
#------------------------------------------------------------------------

#try min 10 for clade size 
treeSt <-  trestruct(resPoly, minCladeSize = 150, minOverlap = -Inf, nsim = 10000, #100000 didn't change anythin
                     level = 0.05, ncpu = 1, verbosity = 1) 
#20 gets 12 groups - too many?
  #50-200 stable 3 groups. >200 1 group (nearly 50% of the sequences)
treeSt_df <- as.data.frame(treeSt)

plot(treeSt, use_ggtree = TRUE)

#----------------------------
#subset tree - went with three groups 
#---------------------------


tokeepc1<-setdiff(b$tip.label,names(treeSt$clustering)[which(treeSt$clustering==2)])
Clade1nex<-drop.tip(b,tokeepc1)
Clade1<-drop.tip(beast,substr(tokeepc1,2,nchar(tokeepc1)-1))


# going forward in the tree to remove weird outliers
ggtree(Clade1) + geom_tiplab(size=2)+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 

pc1 <- ggtree(Clade1, mrsd="2020-03-24") + theme_tree2()+ 
  geom_tiplab(align=FALSE, linetype='dashed', linesize=0.5, size=1) +  
  geom_range("length_0.95_HPD", color='red', size=2, alpha=.5,) + #gets length estimates 
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior)> 0.9, 
                 x=branch), vjust=0)
pc1

Clade1_names <- as.data.frame(get_taxa_name(tree_view = NULL, node = NULL))

#clade 2

tokeep<-setdiff(b$tip.label,names(treeSt$clustering)[which(treeSt$clustering==3)])

Clade2nex<-drop.tip(b,tokeep)
Clade2<-drop.tip(beast,substr(tokeep,2,nchar(tokeep)-1))

#basic tree
ggtree(Clade2nex) + geom_tiplab(size=2)+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
#more annotation added
pc2 <- ggtree(Clade2, mrsd="2020-03-24") + theme_tree2()+ 
  geom_tiplab(align=FALSE, linetype='dashed', linesize=0.5, size=1) +  
  geom_range("length_0.95_HPD", color='red', size=2, alpha=.5, nodes=730) + #gets length estimates 
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior)> 0.9, 
                 x=branch), vjust=0)
pc2

Clade2_names <- as.data.frame(get_taxa_name(tree_view = NULL, node = NULL))
Clade2_names$seq <- as.factor(Clade1_names$seq)

Clade2_names$seq <- gsub(("'"), (""), Clade1_names$seq)


Clade3 <- tree_subset(beast , node=820,levels_back = 0) 
Clade3nex <- tree_subset(b, node=820,levels_back = 0)#useful for analyses below

ggtree(Clade3nex) + geom_tiplab(size=2)+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3) #complete tree

pc3 <- ggtree(Clade3, mrsd="2020-03-24") + theme_tree2()+ 
  geom_tiplab(align=FALSE, linetype='dashed', linesize=0.5, size=1) +  
  geom_range("length_0.95_HPD", color='red', size=2, alpha=.5, nodes=730) + #gets length estimates 
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior)> 0.9, 
                 x=branch), vjust=0)
pc3

Clade3_names <- as.data.frame(get_taxa_name(tree_view = NULL, node = NULL))


#----------------------------
#Effective population size esimates
#---------------------------

#overall - all clades included. lengthout means ne changes every 3 days or so

GlobalBSP <- BNPR(b, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01,
                   beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                   derivative = FALSE, forward = TRUE)

plot_BNPR(GlobalBSP)

#ps helps account for preferential sampling -not needed in this case as sampling has been consistent enough?

BSpsClade1<- BNPR_PS(Clade1nex, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01, #lenghout
                     beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                     derivative = FALSE, forward = TRUE) #with pref sampling
BSpsClade1a<- BNPR(Clade1nex, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01,
                   beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                   derivative = FALSE, forward = TRUE)

plot_BNPR(BSpsClade1)
plot_BNPR(BSpsClade1a)

BSpsClade2a<- BNPR(Clade2nex, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01, #lengthout = number of grid points
                     beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                     derivative = FALSE, forward = TRUE)

plot_BNPR(BSpsClade2a)

BSpsClade3a<- BNPR(Clade3nex, lengthout = 35, prec_alpha = 0.01, prec_beta = 0.01,
                     beta1_prec = 0.001, fns = NULL, log_fns = TRUE, simplify = TRUE,
                     derivative = FALSE, forward = TRUE)

plot_BNPR(BSpsClade3a)


#----------------------------
#Skygrowth models
#---------------------------

#---------------------------
#Global models (ie. using all sequences rather than partitions)
#---------------------------
#fit with mcmc  - similar results were found with the map verison

globalgrowth <- skygrowth.mcmc(b, res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 2e+07, control=list(thin=1e3) ) 

growth.plot( globalgrowth)+theme_bw()

R.plot(globalgrowth, forward=TRUE, gamma=0.90)+theme_bw()

globalMCMC <- as.mcmc(cbind(globalgrowth$growthrate[,1:(ncol(globalgrowth$growthrate)-1)],globalgrowth$ne,globalgrowth$tau))
effectiveSize(globalMCMC)
save(globalgrowth$, file="global_covid")
load("global_covid")

#---------------------------
#Lineage A
#---------------------------

mcmcfit_c1 <- skygrowth.mcmc(Clade1nex, res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 2e+07, control=list(thin=1e3) ) 
growth.plot( mcmcfit_c1 )+theme_bw()

R.plot(mcmcfit_c1 , forward=TRUE, gamma=0.90)+theme_bw()# Maybe something from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0230405

#check mcmc diagnostics -
c1 <- as.mcmc(cbind(mcmcfit_c1$growthrate[,1:(ncol(mcmcfit_c1$growthrate)-1)],mcmcfit_c1$ne,mcmcfit_c1$tau))
effectiveSize(c1)

#---------------------------
#Lineage B
#---------------------------

mcmcfit_c2 <- skygrowth.mcmc(Clade2nex, res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 1e+07, control=list(thin=1e3) ) 
growth.plot( mcmcfit_c2 )+theme_bw()
R.plot(mcmcfit_c2 , forward=TRUE, gamma=0.90)+theme_bw()


c2 <- as.mcmc(cbind(mcmcfit_c2$growthrate[,1:(ncol(mcmcfit_c2$growthrate)-1)],mcmcfit_c2$ne,mcmcfit_c2$tau))
effectiveSize(c2)

#---------------------------
#Lineage C
#---------------------------

mcmcfit_c3 <- skygrowth.mcmc(Clade3nex, res = 35, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 4e+07, control=list(thin=1e5) ) #not sure how to tune this
growth.plot( mcmcfit_c3 )+theme_bw()

#check convergence

c3 <- as.mcmc(cbind(mcmcfit_c3$growthrate[,1:(ncol(mcmcfit_c3$growthrate)-1)],mcmcfit_c3$ne,mcmcfit_c3$tau))
effectiveSize(c3)

save(c3, file="clade3_covid")
load("clade3_covid")
R.plot(mcmcfit_c3 , forward=TRUE, gamma=0.90)+theme_bw()

str(mcmcfit_c3)

#----------------------------
#Extracting sequences from each clade - not used in the manuscript 
#---------------------------
library(ggmsa)
library(Biostrings) 
library(DECIPHER) 
library(tidyverse)
library(ape)

#BiocManager::install("DECIPHER") #etc
dna <- readDNAStringSet("covid_alignment.fasta")

dnadf <- as.data.frame(readDNAStringSet("covid_alignment.fasta"))

dnaNN <- tibble::rownames_to_column(dnadf, "seq")


clade_1seq <- dnaNN %>% 
  filter( dnaNN$seq %in% Clade1_names$seq )


row.names(clade_1seq) <- as.character(clade_1seq$seq)
 SeqNames <- clade_1seq$seq
clade_1seq$seq <- NULL

names(clade_1seq) <- 'dna'
clade_1seqMat <- as.matrix(clade_1seq)


clade1seqB <-DNAStringSet(clade_1seqMat, use.names=TRUE)
clade1seqB@ranges@NAMES <- SeqNames

#make it a string again! problem here with the names
clade_1seqMatNames <- cbind(clade1seqB,  SeqNames)

library(tigger)
writeFasta(clade1seqB, file= "clade1seq.fasta ")
#re-align - not needed

clade_1seqNG <- RemoveGaps(clade1seqB , "all")# why it wants gaps removed?
alignedC1 <- AlignSeqs(clade_1seqNG)

# view the alignment in a browser (optional)
BrowseSeqs(alignedC1, highlight=0) #awful alignment. 

#----------------------------
#Selection using CorMut/Ape
#---------------------------

# dNdS - works but not sure how to sue the reusltant matrix
#convert to right format
clade1dN <- as.DNAbin(clade1seqB ) #may need to add a extra consensus sequence.

C1dnds <-  dnds( clade1dN  , code = 1, codonstart = 1, quiet = FALSE)


library(CorMut) #doesn't work as failing to load. Moving to FUBAR in datamonkey.
file <- system.file("clade1seq.fasta",package="CorMut")

file <- system.file('C:/Users/Nick FJ/Documents/Covid19/clade1seq.fasta')

clade1dN <- seqFormat(file, format = c("fasta"))

result <- kaksCodon(clade1seqB)

fresult <-  filterSites(result)
head(fresult)


#useful code  
ex.dna <- read.dna("covid_alignment.fasta", format = "fasta")
names(ex.dna)
ggmsa(ex.dna) , color = "Chemistry_NT")

#converting dates

x <- c(2019.862)
library(lubridate)
(f <- format(date_decimal(x), "%d-%m-%Y"))
