"""It is advised to look at the IsoformSwitchAnalyzerR main page before proceeding."""


library("IsoformSwitchAnalyzeR")

#Importing Data from Kallisto
kallistoQuant <- importIsoformExpression(parentDir = "/Users/manthospitoulias/Desktop/MSc_Project_2021/kallisto_results/kallisto")
#Make design matrix - 3 conditions: ACP/FETAL/NFPA
sampleID = colnames(kallistoQuant$abundance)[-1]
#3 conditions: ACP/FETAL/NFPA
condition = c("ACP","ACP","ACP","ACP","ACP","ACP","NFPA","NFPA","NFPA","ACP","ACP","ACP","ACP","ACP","FETAL","FETAL","FETAL","ACP","ACP","ACP","ACP","ACP","ACP","ACP")
myDesign_2 <- data.frame(sampleID,condition)

#Create SwitchAnalyzerList
aSwitchList_2 <- importRdata(
  isoformCountMatrix   = kallistoQuant$counts,
  isoformRepExpression = kallistoQuant$abundance,
  designMatrix         = myDesign_2,
  isoformExonAnnoation = "/d/projects/u/pm005/msc_project_2021/reference_gen/gencode.v38.annotation.gtf.gz",
  isoformNtFasta       = "/d/projects/u/pm005/msc_project_2021/reference_gen/gencode.v38.transcripts.fa.gz",
  showProgress = FALSE
)

#Filtering 
aSwitchListFiltered_2 <- preFilter(
  switchAnalyzeRlist = aSwitchList_2,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)
#Testing for Isoform Switches via DEXSeq
aSwitchListAnalyzed_2 <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = aSwitchListFiltered_2,
  reduceToSwitchingGenes=TRUE
)
extractSwitchSummary(aSwitchListAnalyzed_2)

#extract nucleotide sequences, generate fasta files for external analysis
#to add functional consequences
aSwitchListAnalyzed_2<-extractSequence(aSwitchListAnalyzed_2,removeShortAAseq = TRUE,removeLongAAseq = TRUE,pathToOutput = "/path/")

#The following analysis have both command line applications and webserver based. Use whichever. However, web-based might
#require some tweaking to fit the criteria they set.

#cpc2 analysis 
aSwitchListAnalyzed_2 <- analyzeCPC2(
  switchAnalyzeRlist   = aSwitchListAnalyzed_2,
  pathToCPC2resultFile = "/d/projects/u/pm005/msc_project_2021/IsoSwitchAnalyzer/cpc2/results2.txt",
  removeNoncodinORFs   = FALSE   # because ORF was predicted de novo
)
#PFAM analysis
aSwitchListAnalyzed_2<-analyzePFAM(switchAnalyzeRlist = aSwitchListAnalyzed_2,pathToPFAMresultFile = "/d/projects/u/pm005/msc_project_2021/IsoSwitchAnalyzer/pfam/pfam_res.txt")

#Moving to SignalP analysis
aSwitchListAnalyzed_2<-analyzeSignalP(switchAnalyzeRlist = aSwitchListAnalyzed_2,pathToSignalPresultFile = "/d/projects/u/pm005/msc_project_2021/IsoSwitchAnalyzer/signalP/aa_summary.signalp5")

#Moving to IUswitchanalyzed.
#vector of the result files from IU analysis
IU = c("/d/projects/u/pm005/msc_project_2021/IsoSwitchAnalyzer/IUpred2a/aa0.result",
       "/d/projects/u/pm005/msc_project_2021/IsoSwitchAnalyzer/IUpred2a/aa1.result",
       "/d/projects/u/pm005/msc_project_2021/IsoSwitchAnalyzer/IUpred2a/aa2.result",
       "/d/projects/u/pm005/msc_project_2021/IsoSwitchAnalyzer/IUpred2a/aa3.result",
       "/d/projects/u/pm005/msc_project_2021/IsoSwitchAnalyzer/IUpred2a/aa4.result",
       "/d/projects/u/pm005/msc_project_2021/IsoSwitchAnalyzer/IUpred2a/aa5.result",
       "/d/projects/u/pm005/msc_project_2021/IsoSwitchAnalyzer/IUpred2a/aa6.result")
aSwitchListAnalyzed_2<-analyzeIUPred2A(switchAnalyzeRlist = aSwitchListAnalyzed_2,pathToIUPred2AresultFile = IU)


#predicting alternative splicing
aSwitchListAnalyzed_2<-analyzeAlternativeSplicing(switchAnalyzeRlist = aSwitchListAnalyzed_2)

#Predicting Switch Consequences
aSwitchListAnalyzed_2<-extractSequence(aSwitchListAnalyzed_2,removeShortAAseq = TRUE,removeLongAAseq = TRUE,pathToOutput = "/d/projects/u/pm005/msc_project_2021/IsoSwitchAnalyzer/splited_files/attempt_3/",writeToFile = FALSE)

aSwitchListAnalyzed_2<-analyzeSwitchConsequences(aSwitchListAnalyzed_2,dIFcutoff = 0.1)

#Analysis of Individual Isoform Switching

aSwitchListAnalyzedSubset_2<-subsetSwitchAnalyzeRlist(aSwitchListAnalyzed_2,aSwitchListAnalyzed_2$isoformFeatures$condition_1 == "ACP")


#enrichemnt for analysis of splicing events
splicingEnrichment_2 <- extractSplicingEnrichment(
  aSwitchListAnalyzed_2,
  splicingToAnalyze='all',
  returnResult=TRUE,
  returnSummary=TRUE
)

consequenceEnrichment_2 <-extractConsequenceEnrichment(aSwitchListAnalyzed_2)
extract_overlaps_2 <- extractSwitchOverlap(aSwitchListAnalyzed_2,filterForConsequences = TRUE,plotIsoforms = FALSE)



x<-subsetSwitchAnalyzeRlist(aSwitchListAnalyzed_2,aSwitchListAnalyzed_2$isoformFeatures$condition_1 == "ACP")
x
z<-extractTopSwitches(x,n=100,extractGenes = FALSE,sortByQvals = TRUE,inEachComparison = TRUE,filterForConsequences = TRUE)
#example
switchPlot(x, gene = "GRIA4",condition1 = "ACP",condition2 = "NFPA")
#export 100 swit hes
switchPlotTopSwitches(
  switchAnalyzeRlist = aSwitchListAnalyzed_2,
  onlySigIsoforms = TRUE,
  n = 100,                                             # Set to Inf for all
  filterForConsequences = TRUE,
  fileType = "pdf",                                   # alternative is "png"
  pathToOutput = "/d/projects/u/pm005/msc_project_2021/IsoSwitchAnalyzer/combined_ACP"
)


