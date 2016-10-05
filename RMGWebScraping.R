# Belinda Slakman
# DSCS 6020 Term Project - Web Scraping

# This file will scrape the RMG website at http://rmg.coe.neu.edu for chemical data.

# Load required packages
library('RCurl')
library('XML')

# This function scrapes thermodynamic data from the RMG database according to sourceType (either 'libraries' or 'groups') 
# and name (of the library or group).
ScrapeRMGThermo <- function(sourceType, name){
  
  # form the URL by joining the rmg thermo database website with the source type and source name
  url <- paste0("http://rmg.coe.neu.edu/database/thermo/", sourceType, "/", name)
  webpage <- getURL(url, followlocation=TRUE)
  # convert the page into a line-by-line format
  tc <- textConnection(webpage)
  webpage <- readLines(tc) 
  close(tc)
  # get webpage in tree format
  pagetree_parent <- htmlTreeParse(webpage, useInternalNodes = TRUE)
  
  # Figure out how many species there are
  lastSpecLabel <- unlist(xpathApply(pagetree_parent,"//*/div[@id='contents']/table[@class='thermoData']/tr[position()=last()]/td[1]/a", xmlValue))
  lastSpecSplit <- strsplit(lastSpecLabel, ". ")
  numSpecies <- as.numeric(lastSpecSplit[[1]][1])

  # create empty vectors
  label <- rep(NA, numSpecies)
  adj_list <- rep(NA, numSpecies)
  Hf <- rep(NA, numSpecies)
  Sf <- rep(NA, numSpecies)
  Cp_300 <- rep(NA, numSpecies)
  Cp_1000 <- rep(NA, numSpecies)
  
  for(index in 1:numSpecies){
    # create specific URL for species
    url_specific <- paste0(url, "/", index)
    # Make sure it exists
    if (!url.exists(url_specific)) {next}
    webpage <- getURL(url_specific, followlocation=TRUE)
    # convert the page into a line-by-line format
    tc <- textConnection(webpage)
    webpage <- readLines(tc) 
    close(tc)
    
    # get webpage in tree format
    pagetree <- htmlTreeParse(webpage, useInternalNodes = TRUE)
    
    # get the label for the molecule, and add just its name to the list
    whole_label <- unlist(xpathApply(pagetree,"//*/h1",xmlValue))
    split_label <- strsplit(whole_label, ". ", fixed=TRUE)
    label[index] <- split_label[[1]][[2]]
    
    # get the "alt" attribute of the molecule's image, which corresponds to its adjacency list
    adj_list_attr <- unlist(xpathApply(pagetree,"//*/div[@id='contents']/p/img/@alt"))
    adj_list[index] <- adj_list_attr[[1]]
    
    # Check the thermo format. We'll only process if Group additivity
    number_label <- paste0(index, ". ", label[index])
    thermo_type <- unlist(xpathApply(pagetree_parent, paste0("//*/table[@class='thermoData']/tr[td[a='", number_label, "']]/td[3]"), xmlValue)) 
    if (!thermo_type=='Group additivity') {next}
    
    # get enthalpy of formation, entropy of formation, and heat capcity at 300K and 1000K
    Hf_value <- unlist(xpathApply(pagetree, "//*/table[@class='thermoEntryData']/tr[1]/td[@class='value']/span", xmlValue))
    Hf_processing <- strsplit(Hf_value, " ")
    Hf[index] <- paste0(Hf_processing[[1]][[1]], " ", Hf_processing[[1]][[4]])
    
    Sf_value <- unlist(xpathApply(pagetree, "//*/table[@class='thermoEntryData']/tr[2]/td[@class='value']/span", xmlValue))
    Sf_processing <- strsplit(Sf_value, " ")
    Sf[index] <- paste0(Sf_processing[[1]][[1]], " ", Sf_processing[[1]][[4]])
    
    Cp_300_value <- unlist(xpathApply(pagetree, "//*/table[@class='thermoEntryData']/tr[3]/td[@class='value']/span", xmlValue))
    Cp_300_processing <- strsplit(Cp_300_value, " ")
    Cp_300[index] <- paste0(Cp_300_processing[[1]][[1]], " ", Cp_300_processing[[1]][[4]])
    
    Cp_1000_value <- unlist(xpathApply(pagetree, "//*/table[@class='thermoEntryData']/tr[8]/td[@class='value']/span", xmlValue))
    Cp_1000_processing <- strsplit(Cp_1000_value, " ")
    Cp_1000[index] <- paste0(Hf_processing[[1]][[1]], " ", Cp_1000_processing[[1]][[4]])
  } 
  
  # Create data frame of these 10 species
  thermo <- data.frame(label, adj_list, Hf, Sf, Cp_300, Cp_1000, stringsAsFactors = FALSE)
  thermo <- thermo[complete.cases(thermo[,1]),]
  return(thermo)
}

ScrapeRMGKineticsFromLibrary <- function(libraryName){
  
  # form the URL by joining the rmg kinetics library website with the library name
  url <- paste0("http://rmg.coe.neu.edu/database/kinetics/libraries/", libraryName)
  
  webpage <- getURL(url, followlocation=TRUE)
  # convert the page into a line-by-line format
  tc <- textConnection(webpage)
  webpage <- readLines(tc) 
  close(tc)
  # get webpage in tree format
  pagetree_parent <- htmlTreeParse(webpage, useInternalNodes = TRUE)
  
  # Figure out how many species there are
  lastSpecLabel <- unlist(xpathApply(pagetree_parent,"//*/div[@id='contents']/table[@class='kineticsData']/tr[position()=last()]/td[1]/a", xmlValue))
  lastSpecSplit <- strsplit(lastSpecLabel, ". ")
  numSpecies <- as.numeric(lastSpecSplit[[1]][1])
  
  # create empty vectors for reactants and kinetics data
  reactant_1 <- rep(NA, numSpecies)
  reactant_2 <- rep(NA, numSpecies)
  reactant_3 <- rep(NA, numSpecies)
  product_1 <- rep(NA, numSpecies)
  product_2 <- rep(NA, numSpecies)
  product_3 <- rep(NA, numSpecies)
  A <- rep(NA, numSpecies)
  n <- rep(NA, numSpecies)
  E_A <- rep(NA, numSpecies)
  
  for(index in 1:numSpecies){
    # create specific URL for reaction
    url_specific <- paste0(url, "/", index)
    if (!url.exists(url_specific)) {next}
    webpage <- getURL(url_specific, followlocation=TRUE)
    # convert the page into a line-by-line format
    tc <- textConnection(webpage)
    webpage <- readLines(tc) 
    close(tc)
    
    # get webpage in tree format
    pagetree <- htmlTreeParse(webpage, useInternalNodes = TRUE)
    
    reactants <- unlist(xpathApply(pagetree,"//*/div[@id='contents']/table[@class='reaction']/tr/td[1]/a/img/@alt"))
    reactant_1[index] <- reactants[[1]]
    if(length(reactants) > 1) {reactant_2[index] <- reactants[[2]]}
    if(length(reactants) > 2) {reactant_3[index] <- reactants[[3]]}
    
    products <- unlist(xpathApply(pagetree,"//*/div[@id='contents']/table[@class='reaction']/tr/td[3]/a/img/@alt"))
    product_1[index] <- products[[1]]
    if(length(products) > 1) {product_2[index] <- products[[2]]}
    if(length(products) > 2) {product_3[index] <- products[[3]]}
    
    # Check the kinetics format. We'll only process if Arrhenius
    number_label <- paste0(index, ". ")
    kinetics_type <- unlist(xpathApply(pagetree_parent, paste0("//*/table[@class='kineticsData']/tr[td[a='", number_label, "']]/td[5]"), xmlValue)) 
    if (!kinetics_type=='Arrhenius') {next}
    
    kinetics <- unlist(xpathApply(pagetree,"//*/div[@id='contents']/div[@class='math']", xmlValue))
    kinetics_split1 <- strsplit(kinetics, "= ")
    if(length(kinetics_split1[[1]]) < 2) {next}
    kinetics_split2 <- strsplit(kinetics_split1[[1]][[2]], "T")
    A[index] <- kinetics_split2[[1]][[1]]
    if(length(kinetics_split2[[1]]) < 2) {next}
    kinetics_split3 <- strsplit(kinetics_split2[[1]][[2]], " ")
    if(length(kinetics_split3[[1]]) < 9) {next}
    n[index] <- kinetics_split3[[1]][[2]]
    E_A[index] <- kinetics_split3[[1]][[9]]
  }
  
  kinetics <- data.frame(reactant_1, reactant_2, reactant_3, product_1, product_2, product_3, A, n, E_A, stringsAsFactors = FALSE)
  kinetics <- kinetics[complete.cases(kinetics[,1]),]
  return(kinetics)
}

ScrapeRMGSolvation <- function(){
  url <- "http://rmg.mit.edu/database/solvation/libraries/solute"
  webpage <- getURL(url, followlocation=TRUE)
  # convert the page into a line-by-line format
  tc <- textConnection(webpage)
  webpage <- readLines(tc) 
  close(tc)
  # get webpage in tree format
  pagetree_parent <- htmlTreeParse(webpage, useInternalNodes = TRUE)
  
  # Figure out how many species there are
  lastSpecLabel <- unlist(xpathApply(pagetree_parent,"//*/div[@id='contents']/table[@class='solvationData']/tr[position()=last()]/td[1]/a", xmlValue))
  lastSpecSplit <- strsplit(lastSpecLabel, ". ")
  numSpecies <- as.numeric(lastSpecSplit[[1]][1])
  
  # create empty lists for 10 solutes
  label <- rep(NA, numSpecies)
  adj_list <- rep(NA, numSpecies)
  S <- rep(NA, numSpecies)
  B <- rep(NA, numSpecies)
  E <- rep(NA, numSpecies)
  L <- rep(NA, numSpecies)
  A <- rep(NA, numSpecies)
  V <- rep(NA, numSpecies)
  
  for(index in 1:numSpecies){
    url_solvation <- paste0("http://rmg.mit.edu/database/solvation/libraries/solute/", index)
    if (!url.exists(url_solvation)) {next}
    webpage <- getURL(url_solvation, followlocation=TRUE)
    # convert the page into a line-by-line format
    tc <- textConnection(webpage)
    webpage <- readLines(tc) 
    close(tc)
    
    # get webpage in tree format
    pagetree <- htmlTreeParse(webpage, useInternalNodes = TRUE)
    
    # get the label for the molecule, and add just its name to the list
    whole_label <- unlist(xpathApply(pagetree,"//*/h1",xmlValue))
    split_label <- strsplit(whole_label, ". ", fixed=TRUE)
    label[index] <- split_label[[1]][[2]]
    
    # get the "alt" attribute of the molecule's image, which corresponds to its adjacency list
    adj_list_attr <- unlist(xpathApply(pagetree,"//*/div[@id='contents']/p/a/img/@alt"))
    adj_list[index] <- adj_list_attr[[1]]
    
    # Retrieve solvation data
    S_full <- unlist(xpathApply(pagetree,"//*/div[@id='contents']/table[@class='solvationEntryData']/table[@class='solvationEntryData']/tr[1]/td[@class='value']/span"), xmlValue)
    splitS <- strsplit(xmlValue(S_full[[1]]), " ")
    S[index] <- splitS[[1]][1]
    
    B_full <- unlist(xpathApply(pagetree,"//*/div[@id='contents']/table[@class='solvationEntryData']/table[@class='solvationEntryData']/tr[2]/td[@class='value']/span"), xmlValue)
    splitB <- strsplit(xmlValue(B_full[[1]]), " ")
    B[index] <- splitB[[1]][1]
    
    E_full <- unlist(xpathApply(pagetree,"//*/div[@id='contents']/table[@class='solvationEntryData']/table[@class='solvationEntryData']/tr[3]/td[@class='value']/span"), xmlValue)
    splitE <- strsplit(xmlValue(E_full[[1]]), " ")
    E[index] <- splitE[[1]][1]
    
    L_full <- unlist(xpathApply(pagetree,"//*/div[@id='contents']/table[@class='solvationEntryData']/table[@class='solvationEntryData']/tr[4]/td[@class='value']/span"), xmlValue)
    splitL <- strsplit(xmlValue(L_full[[1]]), " ")
    L[index] <- splitL[[1]][1]
    
    A_full <- unlist(xpathApply(pagetree,"//*/div[@id='contents']/table[@class='solvationEntryData']/table[@class='solvationEntryData']/tr[5]/td[@class='value']/span"), xmlValue)
    splitA <- strsplit(xmlValue(A_full[[1]]), " ")
    A[index] <- splitA[[1]][1]
    
    V_full <- unlist(xpathApply(pagetree,"//*/div[@id='contents']/table[@class='solvationEntryData']/table[@class='solvationEntryData']/tr[6]/td[@class='value']/span"), xmlValue)
    splitV <- strsplit(xmlValue(V_full[[1]]), " ")
    V[index] <- splitV[[1]][1]
  }
  
  solvation <- data.frame(label, adj_list, S, B, E, L, A, V, stringsAsFactors = FALSE)
  solvation <- solvation[complete.cases(solvation[,1]),]
  return(solvation) 
}

#prim_thermo <- ScrapeRMGThermo('libraries', 'primaryThermoLibrary')
#Glarborg_C3_kinetics <- ScrapeRMGKineticsFromLibrary('Glarborg/C3')
#solvation <- ScrapeRMGSolvation()