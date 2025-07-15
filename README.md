# The Exo-proteomic Framework for Nitrogen Acquisition from Proteinaceous Organic Matter in the Model Marine Heterotroph *Ruegeria Pomeroyi* DSS-3 

This repository contains information, metadata, and analysis steps for the Master Thesis / Publication by Stemmer, F. (2025), "The Exo-proteomic Framework for Nitrogen Acquisition from Proteinaceous Organic Matter in the Model Marine Heterotroph *Ruegeria Pomeroyi* DSS-3". 

*Last updated*: 07/15/2025

*Corresponding Author*: Fadime R. Stemmer (fstemmer@mit.edu)

## Abstract
Nitrogen-rich proteinaceous compounds are a key component of surface and upper mesopelagic particulate organic matter, and their rapid degradation by heterotrophic bacteria is an important process in releasing bioavailable nitrogen when other nitrogen containing substrates are deficient. Heterotrophs encode a wide variety of hydrolytic enzymes and nutrient transporters, which provide them with the ability to sustain nutrient stresses and inhabit protein-rich niche habitats such as sinking or suspended particles. Despite their important contribution to the oceanic nitrogen cycle, and marine biogeochemical cycles, the exact protein-level response of heterotrophic bacteria to nitrogen limitation in these environments remains unknown. In this study, *Ruegeria pomeroyi* strain DSS-3, a coastal marine heterotrophic bacterium, served as model organism to investigate the shift in protein production and secretion upon nitrogen limitation in the presence of extracellular protein. We tie these observations to extracellular protease rates determined from a fluorescence-based proteolysis assay. Our analyses reveal an increase in extracellular protease activity due to increased transport of these enzymes through the type 1 secretion system into the periplasmic and extracellular space, in response to nitrogen deficiency in DSS-3. Increased amino acid and oligopeptide uptake capabilities as well as higher abundance of amino acid dehydrogenases in the intracellular proteome further highlight this alternative process heterotrophic bacteria employ to challenge nitrogen limitation. Our study contributes to a better understanding of adaptation strategies of marine heterotrophic bacteria to protein-rich and nitrogen-limited niches, such as sinking and suspended particulate organic matter.

## Acknowledgements and Funding
Funding for this project was provided from XXX, to Primary Investigator Saito, Mak. 
Cultures of *R. pomeroyi* were donated by Dr. Mary Ann Moran (UGA) and this project contributes to the mission of the NSF Center for Chemical Currencies of a Microbial Plants [C-CoMP](https://ccomp-stc.org/) 

---
# This Repository
## Directory Structure

- *data*: Folder containing folders with input files for the scripts used for data analysis
    - *growth_OD600*: Folder containing raw data files from OD600 measurements of *R.pomeroyi* DSS-3 cells grown in six different media conditions. Each file is named based on the media that the cells were grown in. ACD = AC Difco marine broth; ProMM = ProMM minimal medium; +Pr: addition of 87µM Casein to the medium; -N: No inorganic nitrogen (NH4Cl) added; -C: No carbon sources added.
        - ACD.txt
        - ProMM.txt
        - ProMM-N-C+Pr.txt
        - ProMM-N+C+Pr.txt
        - ProMM+N-C+Pr.txt
        - ProMM+Pr.txy
    - *growth_FCM*: Folder containing raw data files from flow cytometry cell counts
    - *protease_kinetics*: Folder containing files of fluorescence measurements of FTC-Casein fluorescence upon degradation of cellular proteases
    - *proteomics*: Folder containing input files of normalized peak areas from after spectral processing in Fragpipe and file fomatting
- *scripts*: Folder containing scripts used for the analysis of this dataset
- *figures*: Folder containing figures resulting from the analysis of this dataset using the workflow provided
 
## Raw data availability
Raw data is available upon request to fstemmer@mit.edu. Pre-processing steps are described below under *Material and Methods > Spectral Processing* and the resulting files shared in .txt format in this repository

---

# Material and Methods
This section will mainly focus on the computational analysis conducted in this study. Details about lab methods can be found in the original publication. A summary of the wetlab steps is summarized on following flow chart

<img width="1243" height="698" alt="Screenshot 2025-07-15 at 14 49 42" src="https://github.com/user-attachments/assets/db16bbb9-d598-4981-b4d9-465f840ca3f9" />

## Growth Rates


## Protease Acitivity

## Spectral Processing
Tryptic peptides were analyzed on a Thermo Scientific Orbitrap Astral mass spectrometer using Data Independent Acquisition (DIA). 

### MSFragger - Peptide Identification, Filtering, and FDR Estimation
The generated mass spectra from DIA analysis of *R. pomeroyi* DSS-3 were searched against the translated genomes obtained from NCBI (Accession number: [PRJNA281](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA281/)) using **MSFragger** within [FragPipe (v22.0)](https://github.com/Nesvilab/FragPipe). 

**Peptide Identification Parameters** 
* A precursor mass tolerance: -10-10 ppm
* Fragment mass tolerance: 10 ppm
* Parent tolerance: 10.0 ppm
* Fixed modification: +57 carbamidomethyl cysteine
* Variable modifications: +16 oxidation modification on methionine; +42 acetyl modification on the N-terminus

**Peptide Filtering** 
* Peptide length: 6 - 50 aa
* Precursor charge states +2, +3, +4
* Max missed cleavage: 1

**False discovery rate (FDR) Estimation**
* Decoys generated from a reverse database
* Filter peptides to maximum FDR of 0.01 (1%)
* Combination of individual search results
* Calculation of posterior error probabilities for all peptide spectral matches

### Percolator - Machine learning based FDR filtering
* Reevaluation of Peptide Spectrum Matches (PSMs)
* Filtering of peptides to 1% FDR using PEPs

### DIA-NN - Peptide Quantification and Protein Grouping
* Selection of top-5 fragment ions per peptide
* Calculation of peptide abundances from number of spectra per peptide
* Protein groups were created based on the parsimony principle and grouped if peptides corresponded well to more than one protein. 

## Proteomes - Data Visualization and Statistical Analysis

### Abundance Comparisons 

### Heatmaps

### NMDS

### Volcano Plots

### Normalized peak area - Bar plots

