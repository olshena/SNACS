# Overview of SNACS Analysis Workflow

##### This document provides an overview of the various R functions that comprise SNACS.
##### This overview will reproduce this analysis and visualizations from Experiment 5 in the associated SNACS manuscript, in which 2 patients were multiplexed together. 

First we load SNACS and any dependencies. 

```{r}
library(SNACS)
library(heatmap4)
library(rhdf5)
```
SNACS can be initiated directly from the .hdf5 file that outputs from the open-source alignment and cell calling pipline utilized in this manuscript. Alternatively, if FASTQ files from multiplexed SC DNAseq experiments are aligned using a different pipeline (e.g. Mission Bio Tapestri), genotype, hash antibody signal, total depth, and alternate depth matrices will needto be extracted and then inputted into the SNACSList function.  

```{r}
file_path <- "SNACS5.genotypes.hdf5" # This .hdf5 file is deposited in GEO and available for public use. 

barcodes <- h5read(file_path, "/CELL_BARCODES")
variants <- h5read(file_path, "/VARIANTS")
genotypes <- h5read(file_path, "/GT")
hashes <- h5read(file_path, "/ABS/clr")
hashnames <- h5read(file_path, "/AB_DESCRIPTIONS")
totaldepth <- h5read(file_path, "/DP")
altdepth <- h5read(file_path, "/AD")

mutMat <- matrix(as.numeric(genotypes), nrow = nrow(genotypes), ncol = ncol(genotypes))
rownames(mutMat) <- variants
colnames(mutMat) <- barcodes
mutMat[mutMat == 3] <- NA

colnames(hashes) <- barcodes
rownames(hashes) <- hashnames
hashMat <- hashes

depthTotal <- matrix(as.numeric(totaldepth), nrow = nrow(totaldepth), ncol = ncol(totaldepth))
rownames(depthTotal) <- variants
colnames(depthTotal) <- barcodes

depthAlt <- matrix(as.numeric(altdepth), nrow = nrow(altdepth), ncol = ncol(altdepth))
rownames(depthAlt) <- variants
colnames(depthAlt) <- barcodes
```

Next, we create a SNACSlist object. A SNACSlist object is a simple list-based object that SNACS uses to store data. In addition to creating the SNACSlist object, the **SNACSList function** confirms the input data is in the correct format. The SNACSList object outputs the name of the experiment, the total number of SNPs, the total number of cells, and the number of hashes.

```{r}
exptName <- "Experiment5"
hashColors <- c("red", "blue")

snacsObj <- SNACSList(mut=mutMat,hashes=hashMat,exptName=exptName,hashColors=hashColors,
                      depthTotal=depthTotal,depthAlt=depthAlt)
snacsObj

An object of class "SNACSList"
Experiment name: Experiment5
No. of SNPs: 70476
No. of cells: 2651
Hashes: TS.1, TS.2
SNACSList attributes: mut, hashes, exptName, depthTotal, depthAlt, annHash, annCell, annSNP, processLevels
Output files are saved in "../output" folder
```

The SNACSlist object is then filtered using the **filterData function** to remove very low or very high frequency SNPs as these will be less informative for demultiplexing. The default values for this function remove any single cell with <40% of SNPs genotyped, any single SNP with <40% of single cells genotyped, all SNPs mutated in > 95% all single cells, and all SNPs mutated in <5% of all single cells. These parameters are optional arguments in the filterData function and can be adjusted. In the example below, the total number of SNPs decreased to 77, reflecting that many SNPs were unmutated (i.e., 0) in all single cells. The number of single cells remained the same, at 2651.

```{r}
snacsObj <- filterData(snacsObj=snacsObj)

An object of class "SNACSList"
Experiment name: snacsExpt
No. of SNPs: 77
No. of cells: 2651
Hashes: TS.1, TS.2
SNACSList attributes: mut, hashes, exptName, depthTotal, depthAlt, annHash, annCell, annSNP, processLevels
Output files are saved in "../output" folder
```
The **getBestSNPs function** groups single cells into preliminary groups bashed in hash antibody data, and selects the best SNPs that separate the single cells into these groups. In Experiment 5, this step identified 6 distinguishing SNPs, which are stored in the annSNP dataframe in the SNACslist object.

```{r}
outputFormat <- "pdf"
snacsObj <- getBestSNPs(snacsObj,outputFormat=outputFormat)
```

The **generateAntibodyDensityPlot function** generates plots of the hash antibody distrubtion. The actual antibody expression is plotted in black. To generate the expected bimodal hash antibody distribution, we fit a Gaussian distribution to the left-most actual distribution (red line) and reflect the data to the left of the mode about the mode. 


```{r}
generateAntibodyDensityPlot(snacsObj)
```
<img src="SNACS5_abdistribution.png" alt="Experiment 5 Ab Distribution" width="400"/>

The **imputeMissingMutations function** imputes missing data prior to hierarchical clustering. 

```{r}
snacsObj <- imputeMissingMutations(snacsObj=snacsObj)
```

The **clusterCellsWithSNPdata function** performs hierarchical clustering where the number of cells is the number of constituent samples using the imputated mutation data from the SNPs selected by the getBestSNPs function. The function also generates an initial heatmap showing the results of this clustering with the hash antibody distributions, providing a simple visualization of the clustering result. 

```{r}
outputFileName <- "Experiment5"
snacsObj <- clusterCellsWithSNPdata(snacsObj)
```

<img src="SNACS5_Heatmap1.png" alt="Experiment 5 Heatmap after Clustering" width="400"/>

The **makeSnacsCall function** assigns single cells to a single sample or multiplet. The cell cluster obtained from running clusterCellsWith SNPdata are split into sub-clusters, traversing down the hierarchical tree. Clusters are split if a significant difference is found when comparing the daughter nodes of any current node.  


```{r}
snacsObj <- makeSnacsCall(snacsObj)
```
