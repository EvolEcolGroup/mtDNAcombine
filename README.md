









---
[![CircleCI](https://img.shields.io/circleci/build/github/EvolEcolGroup/mtDNAcombine/master?label=build%20%28master%29&logo=circleci&style=plastic&token=dbf19e075582baa257024cdbe7a75edf96f93517)](https://circleci.com/gh/EvolEcolGroup/mtDNAcombine/tree/master)




# Using `mtDNAcombine`

This vignette describes the `mtDNAcombine` package, an `R` library designed to support comparative analyses of Bayesian Skyline Plot (BSP) population histories based on mtDNA sequence data from multiple studies. 

`mtDNAcombine` includes functions to retrieve, align, summarise, and maniplulate seuquences downloaded from GenBank, as well as generating basic BEAST2 input files.  There are also accessory functions to analyse and plot the outputs of BEAST2 runs.

Below is a flow diagram of the processes and steps in the `mtDNAcombine` pipeline.

![plot of chunk unnamed-chunk-1](/home/eleanor/Documents/Projects/mtDNAcombine/inst/extdata/mtDNAcombine_flow.jpg)


If you are running a LINUX operating system you will need to have the following system dependancies installed before attempting to install `mtDNAcombine`. In Ubuntu 18.04, run the below code in your bash shell: 

```r
sudo apt-get update 
sudo apt-get install libmagick++-dev libudunits2-dev libgdal-dev 
  build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
```

Both in Linux and Windows, you will need to install the `devtools` 'R' package:

```r
install.packages("devtools")
library(devtools)
```

Then you will need to install `mtDNAcombine`:

```r
devtools::install_github("EvolEcolGroup/mtDNAcombine")
```

And finally load it:

```r
library(mtDNAcombine)
```




## Creating input accession number file


To start a project comparing mitochondrial DNA from multiple individuals, species, and studies, we first need to have a list of unique GenBank accessions to explore.  

These accessions can be acquired in multiple ways. A simple method would be to undertake a broad search of GenBank, e.g. open the NCBI webpage with code such as below:

```r
browseURL( "https://www.ncbi.nlm.nih.gov/nuccore")
```

Set search terms along the lines of "birds"[porgn:__txid8782] in the Nucleotide database with the 'Genetic compartments - Mitochondria' box checked.
  

Then click the 'Send to:' drop down, check the 'Complete Record' radio button, then under 'Choose Destination:' the 'File' radio button, and finally, under 'Format:' select 'Accession List' then press 'Create File'


The output produced by GenBank is likely to be a '.seq' file, or, if using a list of accessions gathered in a different way (e.g. compiled by hand), a list of accessions saved in a .csv format is okay as long as there are no headers or row names. 


For this vignette we will work with a fixed set of 335 accessions in the file "vignette_accessions.csv" from the pacakge. This file can be accessed from the extdata directory through the 'system.file' command (see below for an example of its usage).



## An initial sweep of available information 


Firstly, we build a dataframe that contains information on all the genes / sequences associated with each accession numbers to explore what information is available.

```r
path_to_file <- system.file("extdata","vignette_accessions.csv", 
                            package="mtDNAcombine")
GB_data <- build_genbank_df(accession_file_name = path_to_file)
```



Within GenBank, the same single sequence is often associated to multiple features (e.g. 'source', 'gene', and 'CDS'). This is visible on the website, where the same sequence is found under multiple 'Feature' tabs. This means that the same sequence will also be grabbed multiple times when scraping data from GenBank, hence GB_data has 860 observations when it was given 335 accession numbers to search.  

It should also be noted at this stage that `GB_data` only has data for 333 unique accession numbers, two less that the number included in the `vignette_accessions.csv`.  A key feature of the `build_genbank_df` function is that it will silently remove any RefSeq accessions included in the list of accession numbers it is given, if they are idenitified as duplicates.  The RefSeq collection aims to provide a collated and stable set of standard reference sequences for studies from all disciplines to build on.  Drawn from genomes already available in GenBank and other community databases, they may duplicate existing accessions and, for our purposes, need to be removed.  The `load_accessions_list` function, called within `build_genbank_df` silently removes these accesions having first checked that the 'original' accessions they duplicate are already included in the dataset. 


## Tidying up raw information 


To clean the dataset for later analysis, we must remove duplicated entries caused by a single sequence being associated with multiple 'Feature' tabs. 

We must also control for the bredth of possible names used to describe a single gene as individual studies/groups/projects upload data to GenBank using a range of possible synonyms, abreviations, and misspellings. The first step is to standadise nomenclture across the dataframe by converting gene names to a user defined set of 'standard' nomenclatures. 

By default, the `standardise_gene_names` function loads a file containing alternate abreviations, common misspellings, and other frequent errors for 18 commonly sequenced mitochondrial genes.  The user can upload a custom file by specifying the different file as the second variable in the function: `standardise_gene_names(df_to_update, names_to_replace)`
 

```r
GB_data <- standardise_gene_names(df_to_update = GB_data)
```

Then we remove the duplicates.

```r
GB_data <- remove_duplicates(df_to_update = GB_data)
nrow(GB_data)
```

```
## [1] 500
```

We do not expect the number of accessions we're exploring and number of observations to match at this stage. Indeed, in this example we see that, from the orignal 335 accession numbers in "vignette_accessions.csv", the GB_data data frame now has 500 observations.  This is because the script captures ALL genes associated with each given accession.  Every submission to GenBank recieves a unique accession number but these individual submissions can contain data for anything from a single gene through to whole genome data.  



## Check what information is available


Firstly, for what genes are there data?

```r
GB_genes <- droplevels(as.data.frame(unique(GB_data$gene_name)))
```

These data are still messy.

We can tidy the data a little by removing some gene names that are unlikely to be useable or comparable with other sequences e.g. removing any names that are just numbers, removing names over a certain length, and/or dropping other common unwanted patterns.

```r
GB_genes <- as.data.frame(GB_genes[grep("[[:alpha:]]", 
                                        GB_genes$`unique(GB_data$gene_name)`), ])

colnames(GB_genes) <- "gene_name"

GB_genes <- as.data.frame(GB_genes[!nchar(as.character(GB_genes$gene_name)) > 20, ])

patt <- c("RNA|unknown|trn|ATP|similar|central|duplicated|tandem|repeated|other|myoglobin")

colnames(GB_genes) <- "gene_name"

GB_genes <- droplevels(as.data.frame(GB_genes[!grepl(patt, x = GB_genes$gene_name), ]))
```


Secondly; For what species are there data? 

```r
GB_species <- droplevels(as.data.frame(unique(GB_data$sci_nam)))
```

At the moment:

```
## [1] "GB_genes has 13  unique gene names while GB_species has 5 unique species names"
```


However, looking at these data in more detail shows us that there are still some spurious entries being included. For example, GB_species includes both *Motacilla alba* and *Motacilla alba alboides*, a recognised subspecies but, in this instance, data that we want to group with *Motacilla alba* more broadly.

```
##             Unique Names
##           Motacilla alba
##     Picoides tridactylus
##        Calidris maritima
##      Pinicola enucleator
##  Motacilla alba alboides
```



## Clean up the species names


As stated previously, the amount of freedom in the formating of descriptive information associated with GenBank submissions means that individual submissions can vary the chosen species name, inlcuding using different levels of detail for taxonomic rank. For example, some studies may use subspecies names where others choose not to. Subspecies recognition is frequently debated and sometimes we may want to group together samples with names that aren't an exact match. A 3 word name, not a 2 word scientific name, is a simple pattern to recognise subspecies and we exploit that here. 


The `check_poss_synyms` function returns a list of scientific names that are longer than 2 words and these names will be outputted as a .csv file; "poss_synyms.csv"

```r
poss_synyms <- check_poss_synyms(data = GB_data)
```

If, after investigation, any of these species names need updating or altering then they can be edited within the .csv file.  As long as the edited file is saved in the same format, then the `standardise_spp_names` function will reload and integrate any updated names.

For this example, the edited file has been called "poss_sysnyms_updated.csv".

```r
GB_data <- standardise_spp_names(data = GB_data, 
              new_names_file = system.file("extdata", "poss_synyms_updated.csv",
                                           package="mtDNAcombine"))
```

After updating the species names samples for *Motacilla alba alboides* are now labelled as *Motacilla alba* and, therefore, group together for downstream analysis. 

```
##          Unique Names
##        Motacilla alba
##  Picoides tridactylus
##     Calidris maritima
##   Pinicola enucleator
```



## Filter dataframe 


Different ways of creating the original list of accession numbers result in different types of noise being introduced to the dataframe.  To retain only sequence data from relevant genes, the dataframe needs to be filtered using the `gene_of_interest` function.

Here we want to look at the *ND2* gene so we subset the dataframe to include information on the gene of interest

```r
GB_by_gene <- gene_of_interest(gene = "ND2", data = GB_data)
```



## Extract the available raw sequence data


By simply using the `get_GB_sequence_data` function and the curated accession list, we can now download raw sequence data associated with the specific gene of interest.  

```r
GB_with_SeqDat <- get_GB_sequence_data(accessions_of_interest = GB_by_gene, 
                      gene = "ND2", new_names_file = "poss_synyms_updated.csv")
```



For this vignette, the GB_with_SeqDat file should now be 333 observations with 8 variables. 

```r
nrow(GB_with_SeqDat)
## [1] 333
ncol(GB_with_SeqDat)
## [1] 8
```



## Store out key data files


At this stage it might be helpful to store summary details on the data as well as keeping all the raw, unaligned, sequence data for each species / gene combination. The `export_details` function writes summary details to .csv files while the `export_sequences` function writes out individual .fasta files for each dataset.

```r
export_details(data = GB_with_SeqDat)

export_sequences(data = GB_with_SeqDat)
```



## Aligning the sequence data 


We have now managed to generate a set of sequences from multiple species covering one gene. However, for each species, these sequences likely come from multiple independent studies and frequently differ in the gene region they analyse.  

Previously, the  `export_sequences` function wrote out a .fasta file of raw, unaligned, sequence data for each species / gene combination, starting the file name with the regular expression 'FOR_ALIGNMENT'. We exploit this pattern to capture the list of file names to explore. 

```r
alignment_files <- list.files(pattern="FOR_ALIGNMENT")
```


For each species, the sequence data needs to be aligned so that we can capture comparable regions of the genome common to each sample.  This is done within the `align_and_summarise` function using the ClustalW algorithm, removing any columns with blanks or ambiguos calls.

```r
align_and_summarise(alignment_files = alignment_files, 
                    max_haps_found_together = 2, minbp = 200)
```


## Diagnostic plots


### Histogram

Depending on the quality/consistency of the raw sequence data, this step can result in a dramatic reduction in the number of base pairs left in the DNA string.  For example, where one or two sequences are very short, or the section of the genome sequenced is different, the overlap between data from separate studies can be very small. In some instances, removal of one or two sequences before alignment could result in a more informative data set.

The impact of the alignment/trimming process is summarised in a diagnostic histogram plot, offering a visual way to identify cases where it would be advantageous to look at the raw data in more detail. The histogram bars show frequency and sequence length of raw, unaligned data and the red line shows the length of the aligned sequences after cropping to the longest section common to all samples. 

![plot of chunk unnamed-chunk-9](/home/eleanor/Documents/Projects/mtDNAcombine/inst/extdata/hist_Pinicola_enucleator.png)

Here we see that, in the *Pinicola enucleator* dataset, the majority of samples have been heavily cropped due to the inclusion of one, shorter, sequence. In this instance, it may be worth reviewing the decision to include the single, much shorter, 450 base pair sample. 

![plot of chunk unnamed-chunk-10](/home/eleanor/Documents/Projects/mtDNAcombine/inst/extdata/hist_Calidris_maritima.png)

Alternatively, the *Calidris maritima* histogram shows that, whilst a few longer sequences have been trimmed by a couple hundered base bairs, the majority of the sequences are being used at nearly full length. This alignment and crop seems good. 



### Network diagram

The `align_and_summarise` function also produces a haplotype network diagram which helps visuliase the level of structure in a population/sample set.  


![plot of chunk unnamed-chunk-11](/home/eleanor/Documents/Projects/mtDNAcombine/inst/extdata/Net_Picoides_tridactylus.png)
Here is an example of a network diagram for data from *Picoides tridactylus*. Plots like these help to quickly flag if there are any extreme outliers in the dataset or if the population is heavily structured. 


### Haplotype frequency 

In datasets that include samples from studies which have uploaded a single representative haplotype, instead of creating a new accession for every sample, an aligned sequence output file is not generated.  Instead, the papers associated with the unique haplotypes are listed in a .csv file along with the species name; the file is is written out as "More_info_df.csv". 

For each of the populations recorded in "More_info_df.csv", the `align_and_summarise` function also creates a file containing each accession number and a frequency column.  These files all have the regular pattern "MAGNIFY"" in the file name.  Orignal published papers must be tracked down to confirm details of sampling frequency, new values can then be recorded in the "freq" column.  Once updated the file needs to be saved in the same format. In cases where sampling frequency data is not available samples must be excluded.



In this vignette-dataset accessions for *Calidris maritima* and *Motacilla alba* are flagged as needing further investigation. Exploration of the original published papers suggest that data for *Motacilla alba* have indeed been uploaded at sampled frequency and this was essentially a "false alarm". Therefore, we don't want to alter this data so the "MAGNIFY_Motacilla_alba.csv"" file can be left as it is, with a default "freq" column value of "1".

However, exploration of the paper *'A review of the subspecies status of the Icelandic Purple Sandpiper Calidris maritima littoralis'* shows that only unique haplotype sequences were uploaded, rather than a new accession being created for every sample.  Therefore, this dataset needs to be manipulated to get to the original sampled frequency. 

For the puruposes of this vignette we have created an updated "MAGNIFY_Calidris_maritima.csv" file (found in ../extdata/) which already contains the values for the number of times each haplotype was sampled in the population.

```r
magnify_file_list <- list.files(pattern="MAGNIFY")

mag_df <- magnify_to_sampled_freq(magnify_file_list = magnify_file_list)
```

After updating the frequency information the sequences are processed as before - haplotype networks are drawn and .fasta files of the aligned sequence data written out.


To keep accurate summary information of the datasets avalible, we now need to combine the original `info_df` and the newly created `mag_df`.  This will give an updated .csv file that contains information on all the data sets we are working with.

```r
info_df <- updating_info_df(original_df = "Info_df.csv", new_df = mag_df)
```



## Filtering by rules


Data inclusion criteria will vary between studies; there will never be a "one-size-fits-all" set up.  Factors such as species life history, species population history, data availablity, data quality, and even broadly the project aims will influence what data are informative.

The following filtering steps are based on a series of rules built around avain mtDNA.

Firstly, we want to drop populations with insufficent sequence data. This includes data with insufficent number of bases, low numbers of haplotypes, low sample size.

```r
info_df <- drop_low_sample_size(info_df = info_df, min_sample = 7)

info_df <- drop_low_haplo_number(info_df = info_df, min_haps = 6 )

info_df <- drop_low_sequence_length(info_df = info_df, min_length = 600)
```


After applying these filters we are left with curated datasets from two species.  We then want to remove any extreme outliers, considered here to be single samples separated from the nearest haplotype with >30 mutations on a branch.  The function `outliers_dropped` writes out an updated version of `info_df` but doesn't return it.  Therefore we need to read in the new version from the working directory.

```r
what_gets_dropped <- outliers_dropped(max_mutations = 30, info_df = info_df)

info_df <- read.csv("Info_df.csv")
```

```r
what_gets_dropped
```

```
##      outlier_accession spp_name            
## [1,] "EU166960.1"      "new_Motacilla_alba"
```




At this point:

* all the orignal accessions from the accession list have been processed, 
* raw sequence data from relevant sections of the genome have been captured,
* sequences have been aligned, 
* low resolution/low quality data have been rejected, 
* outliers have been removed,
* cleaned sequence data has been written out for use by additional tools.

We want to use these processed data to set up BEAST runs.  In order to speed up the process, and reduce the opportunity for human error, we limit the amount of manual set up required by creating basic BEAUti files with the following code. These files will still require some degree of editing in the BEAUti GUI.


The dataset files have been given the prefix "ALIGNED_", making them easy to find. After editing (e.g. dropping outliers), any new versions of aligned data have been given the tag "new_ALIGNED" and should be used in preference to the original files.

```r
aligned_files <- list.files(pattern="ALIGNED")


superseeded <- NULL
for(n in 1:length(aligned_files)){
  if (substr(aligned_files[n],1,3)=="new"){
    superseeded <- rbind(gsub("new_",'',aligned_files[n]), superseeded)
  }
}
aligned_files <- aligned_files[!aligned_files%in%superseeded]
```


## Build basic xml files 

```r
setup_basic_xml(gene_name = "ND2", aligned_files = aligned_files)
```


Once created, the .xml files will need to be manually edited.  For example, setting up the use of bModelTest - at the time of writing not yet an available option in the `babette` package.



# Exploring outputs


Once BEAST runs are completed we need to explore convergence, ESS values, and other metrics.  

Here we present a pipeline for handling outputs from the software package Tracer.

An example BEAST .xml input file can be found in ./extdata/ND2_Carpodacus_erythrinus_BEASTinput.xml . After running this file in BEAST v2 4.6 we used Tracerv1 to format output data for export.  The resulting file is stored as ./extdata/ND2_Carpodacus_erythrinus_TracerOut.txt

We will use this output file to explore simple plotting/visulisation.


## Plotting


A quick look at the structure of the output from Tracer taking only complete rows (NAs can occour at the end of the file but cause issues in later processing).

```r
data <-read.table(system.file("extdata","ND2_Carpodacus_erythrinus_TracerOut.txt", 
                                    package = "mtDNAcombine"), skip=1, header=T)
data<-data[complete.cases(data),]
head(data)
```

```
##       Time    Mean  Median    Upper   Lower
## 1   0.0000 7952376 7075067 16749959 3933864
## 2 145.5466 7955168 7076337 16713627 3953354
## 3 291.0932 7958742 7079056 16695231 3971612
## 4 436.6397 7958827 7080898 16669858 3975339
## 5 582.1863 7950289 7079365 16601677 3975339
## 6 727.7329 7950633 7082537 16579575 3982642
```


A simple coloured plot can now be created using the code below.

Use the file name as the title

```r
file_name <- "ND2_Carpodacus_erythrinus1_TracerOut.txt"
plot_title <- "Common rosefinch"
```

We want to plot the median (log scale) as well as plotting the HPD interval as a coloured polygon.  
If analysing data from multiple genes it can be helpful to differentiate the plots on the basis of gene type.  Here this is done by colouring the HPD according to the gene.
 

```r
plot(log10(data[,3])~data[,1],type="n",ylim=c(3.5,7.5),xlim=c(0,60000), yaxt="n", 
     yaxs="i", ylab = expression("Pop. Size (Log'[10]*')"), xaxs="i", 
     xlab = "Time since present day (yrs)")
axis(2, at=c(4,5,6,7), labels = c("1.E4","1.E5","1.E6","1.E7"), las=2, adj=1,
     cex.axis=0.82)


gene <- substr(file_name, start = 1, stop=4)

for(i in 1:length(data[,1])-1) {
  x<-c(data[i,1],data[i+1,1],data[i+1,1],data[i,1])
  y<-log10(c(data[i,4],data[i+1,4],data[i+1,5],data[i,5]))
        
      
      if(gene=="cytb"){
         polygon(x,y,col="plum3", border="plum3")
       }else{
        polygon(x,y,col="darkolivegreen3", border="darkolivegreen3")
       }
     }
         
#median value as a dashed line 
points(log10(data[,3])~data[,1],type="l", lty=2, lwd=2)
#edge the HPD interval by plotting the upper and lower 95% HPD
points(log10(data[,4])~data[,1],type="l")
points(log10(data[,5])~data[,1],type="l")
        

title(main = plot_title)
```

<img src="./inst/extdata/plotting_median-1.png" title="plot of chunk plotting median" alt="plot of chunk plotting median" style="display: block; margin: auto;" />
