# GEOmetric

## Step 0
Select relevant GSEs based on key-word searching and reading description.

## Step 1
Download GSEs of interest and make a structured annotation table describing samples and other experiment data. Significance:
1. The description of a GSE gives only a conceptual guide to whether it is good for our purposes. This process helps further determine if a GSE is suitable and what kind of contrasts can be generated given the samples and experimental design. 
2. The annotations make it easier to feed the datasets to the PCA and limma in a standard and reproducible manner. 
3. These annotations are ultimately required for the NURSA summary which is submitted with each result.

A review file contains these annotation columns. **Bold** indicates what should be filled out for the PCA code to work. *Italics* indicates the additional annotations not immediately required for Step 2, but required for the contrast and limma code.  

*xpType*: 1C (1Channel), 2CDirect (2ChannelDirect), 2CCommonRef (2ChannelCommonReference), Htseq (HighthroughputSequencing)
Note: Past files will also have 2CDirectDyeSwap, which is not a necessary distinction when it comes to the code doing the statistical tests, but it's maybe nice for a human to know?

*RefDye*: Reference dye (Cy3 or Cy5). Not necessary for 1C or Htseq experiments.

*Group*: Experimental group. Labels can be automatically assigned based on Node, NodeFunction, BSM and BSMDCD.

*Ignore*: In rare cases, further examination indicates that some samples should not be analyzed. 0 = F, 1 = T.

*XP*: Experiment number. In rare cases, GSEs can actually contain two experiments in one GSE dataset. They would be analyzed separately, so we would want to know which samples belong to which experiment.

isCTRL: Control group. Not totally necessary or always useful, but it's nice to have.

**Batch**: Can use numbers or a text label to indicate batch. 

**Node**: Gene

**NodeFunction**: KD (Knock-down), KI (Knock-in), OE (Overexpression)
Note: Previous used labels include: OE (Overexpression), LOF (Loss of Function)

**BSM**: Biologically active Small Molecule (ligand)

**BSMDCD**: Time points	

BSMDCD2: Dose, concentrations. In practice, group assignments are not based on BSMDCD2 because experiments usually vary time points and
not concentrations for a ligand.  

## Step 2 
Statistical checks of dataset.
1. PCA -- batch effects, mislabeling?
2. Boxplots of samples -- arrays have been normalized, no RNA degradation? 
3. Data are on logged data -- almost always yes (since limma expects log data)

## Step 3
Create contrasts

## Step 4
Limma -> results

## Step 5
Map biosample ids and create a dataset summary file with descriptions for each submission.



