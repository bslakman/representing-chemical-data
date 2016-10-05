# Belinda Slakman
# DSCS 6020 Term Project - Relational database creation and querying

# This file will create a database from data scraped from the RMG website using the SQLite package in R, 
# and query it for chemical data of interest.

# First, gather some data into data frames from the web scraping algorithms.
# There is only 1 solvation library; scrape 3 kinetics libraries and 3 thermo libraries.
solvation <- ScrapeRMGSolvation()
sulfur_thermo <- ScrapeRMGThermo('libraries', 'SulfurLibrary')
primary_thermo <- ScrapeRMGThermo('libraries', 'primaryThermoLibrary')
DFTQCI_thermo <- ScrapeRMGThermo('libraries', 'DFT_QCI_thermo')
sulfur_DMS_kinetics <- ScrapeRMGKineticsFromLibrary('Sulfur/DMS')
Glarborg_C3_kinetics <- ScrapeRMGKineticsFromLibrary('Glarborg/C3')
GRI_kinetics <- ScrapeRMGKineticsFromLibrary('GRI-Mech3.0')

# Create tables for the database
# Thermo table
label <- c(sulfur_thermo[1:nrow(sulfur_thermo), 1], primary_thermo[1:nrow(primary_thermo), 1], 
           DFTQCI_thermo[1:nrow(DFTQCI_thermo), 1])
thermoLibraryName <- c(rep("SulfurLibrary",nrow(sulfur_thermo)), rep("primaryThermoLibrary",nrow(primary_thermo)), 
                       rep("DFT_QCI_thermo",nrow(DFTQCI_thermo)))
thermo_data <- rbind(sulfur_thermo[1:nrow(sulfur_thermo), 3:6], primary_thermo[1:nrow(primary_thermo), 3:6], 
                 DFTQCI_thermo[1:nrow(DFTQCI_thermo), 3:6])
thermo_table <- data.frame(label, thermoLibraryName, thermo_data)

# Species table with adj list and solvation data
adjlist <- c(sulfur_thermo[1:nrow(sulfur_thermo), 2], primary_thermo[1:nrow(primary_thermo), 2], 
             DFTQCI_thermo[1:nrow(DFTQCI_thermo), 2])
S <- rep(NA, length(adjlist))
B <- rep(NA, length(adjlist))
E <- rep(NA, length(adjlist))
L <- rep(NA, length(adjlist))
A <- rep(NA, length(adjlist))
V <- rep(NA, length(adjlist))
for (molecule in 1:length(label)){
  for(solute in 1:length(solvation)) {
    if (label[molecule]==solvation[solute,'label']) {
      S[molecule] <- solvation[molecule, 'S']
      B[molecule] <- solvation[molecule, 'B']
      E[molecule] <- solvation[molecule, 'E']
      L[molecule] <- solvation[molecule, 'L']
      A[molecule] <- solvation[molecule, 'A']
      V[molecule] <- solvation[molecule, 'V']
    }
  }
}
solvation_table <- unique(data.frame(label, adjlist, S, B, E, L, A, V))

# Kinetics table
# Primary key for this table includes an index for each reaction in the library as well as the library name
kineticsLibraryName <- c(rep("Sulfur/DMS",nrow(sulfur_DMS_kinetics)), rep("Glarborg/C3",nrow(Glarborg_C3_kinetics)), 
                       rep("GRI-Mech3.0",nrow(GRI_kinetics)))
reaction_index <- rep(NA, nrow(sulfur_DMS_kinetics)+nrow(Glarborg_C3_kinetics)+nrow(GRI_kinetics))
count <- 1
for (entry in 1:nrow(sulfur_DMS_kinetics)) {
  reaction_index[count] <- entry
  count <- count + 1
}
for (entry in 1:nrow(Glarborg_C3_kinetics)) {
  reaction_index[count] <- entry
  count <- count + 1
}
for (entry in 1:nrow(GRI_kinetics)) {
  reaction_index[count] <- entry
  count <- count + 1
}
kinetics_data <- rbind(sulfur_DMS_kinetics[1:nrow(sulfur_DMS_kinetics),], Glarborg_C3_kinetics[1:nrow(Glarborg_C3_kinetics),], 
                     GRI_kinetics[1:nrow(GRI_kinetics),])
kinetics_table <- data.frame(reaction_index, kineticsLibraryName, kinetics_data)

# Create SQL database
db <- dbConnect(SQLite(), dbname="RMG.sqlite")
dbWriteTable(conn=db, name="Thermo_Data", value=thermo_table, row.names=FALSE, overwrite=TRUE)
dbWriteTable(conn=db, name="Molecule_Data", value=solvation_table, row.names=FALSE, overwrite=TRUE)
dbWriteTable(conn=db, name="Kinetics_Data", value=kinetics_table, row.names=FALSE, overwrite=TRUE)

# Query the database
# Which thermo libraries does OH appear in and what are the thermodynamic parameters for each?
OH_thermo <- dbGetQuery(db, "SELECT thermoLibraryName, Hf, Sf, Cp_300, Cp_1000 FROM Thermo_Data 
                             WHERE Thermo_Data.label = 'OH'")

# Which species appear in more than 1 thermo library?
duplicate_thermo <- dbGetQuery(db, "SELECT label, thermoLibraryName, Hf, Sf, Cp_300, Cp_1000 FROM Thermo_Data 
                               WHERE label IN (SELECT label FROM Thermo_Data 
                               GROUP BY label HAVING COUNT (label) > 1) ORDER BY label")

# How many reactions does H2 appear in?
H2_reactions <- dbGetQuery(db, "SELECT reaction_index, kineticsLibraryName, reactant_1, reactant_2, reactant_3, product_1, product_2, product_3 
                            FROM Kinetics_Data 
                            WHERE reactant_1 = 'H2' OR reactant_2 = 'H2' OR reactant_3 = 'H2' OR product_1 = 'H2' OR product_2 = 'H2' OR product_3 = 'H2'")

# How many reactions in each library are Arrhenius? (Meaning it has values for A, Ea)
numArrhenius <- dbGetQuery(db, "SELECT kineticsLibraryName, COUNT(*) FROM Kinetics_Data 
                           WHERE A != 'NA' AND E_A != 'NA' GROUP BY kineticsLibraryName")

# Which molecules from primaryThermoLibrary have solvation data?
primThermoWithSolvation <- dbGetQuery(db, "SELECT Molecule_Data.label, adjlist, S, B, E, L, A 
                                      FROM Thermo_Data INNER JOIN Molecule_Data 
                                      WHERE Thermo_Data.thermoLibraryName = 'PrimaryThermoLibrary' AND Thermo_Data.label = Molecule_Data.label AND Molecule_Data.S != 'NA'")