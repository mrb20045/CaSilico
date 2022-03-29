# CaSilico R Package

# Introduction

The efficiency of CRISPR-Cas system is highly depends on well-designed CRISPR RNA (crRNA). To facilitate the use of various types of CRISPR-Cas systems by a wide range of researchers, there is a need for development of computational tools to design crRNAs which cover different CRISPR-Cas systems with off-target analysis capability. Numerous crRNA design tools have been developed but nearly all of them are dedicated to design crRNA for genome editing. 


Hence, we developed a tool matching the needs of both beginners and experts, named CaSilico, which was inspired by the limitations of the current gRNA design tools for designing crRNAs for Cas12, Cas13, and Cas14 CRISPR-Cas systems. Using a list of important features such as mismatch tolerance rules, self-complementarity, GC content, frequency of cleaving base around the target site, target accessibility and PFS (protospacer flanking site) or PAM (protospacer adjacent motif) requirement, CaSilico searches all potential crRNAs in a user-input sequence. Considering these features help users to rank all crRNAs for a sequence and make an informed decision about whether a crRNA is suited for an experiment or not. Our tool is sufficiently flexible to tune some key parameters governing the design of crRNA and identification of off-targets, which can be led to increases the chances of successful CRISPR-Cas experiments.

CaSilico outperforms previous crRNA design tools in the following respects: 1) supporting any reference genome/transcriptome for which a FASTA file is available; 2) designing crRNAs that simultaneously target multiple sequences through conserved region detection among a set of sequences; 3) considering new CRISPR-Cas subtypes; 4) reporting a list of different features for each candidate crRNA, which can help the user to select the best one. Given these capabilities, CaSilico addresses end-user concerns arising from the use of sophisticated bioinformatics algorithms and has a wide range of potential research applications in different areas especially design of crRNA for pathogen diagnosis. 

CaSilico was successfully applied to design crRNAs for different genes in SARS-CoV-2 genome, as some of the crRNAs have been experimentally tested in the previous studies.
![s](https://user-images.githubusercontent.com/9910942/158050606-9c592f1c-c0a4-4a4e-8cf8-ded8e0c7e7b6.png)


CaSilico workflow. (A) CaSilico accepts a single or a set of DNA or RNA sequences to be scanned for crRNA designing. (B) When more than one sequences is given as input, the conserved regions among them automatically detect considering conservation threshold and one of the two different approaches for identifying conserved regions. (C) A sliding window (stride of 1 nt) is employed across the single sequences or conserved region of multiple sequences to specify potentail target sites. (D) CaSilico applies multiple criteria for crRNA designing, performs off-target analysis and returns outputs in an interactive graphical interface and some files such as MSA and secondary structure (E & F).





# Installation

```
library("devtools")

install_github("mrb20045/CaSilico")
```



# A quick example to use CaSilico package
```
library("CaSilico")



#Run with TargetFasta

data <- paste0(system.file(package = "CaSilico"),"/data/3D.fasta")

CaSilico(ResultsFolder="Example",
         TargetFasta=data,
         TargetAccession=NULL,
         CrisprTypes=c("casVI_A"),
         ConservationMethod = 1,
         ConservationThreshold=0.98,
         OffTarget = F,
         OffAsk = F)




#Run with TargetFasta and off-target analysis (local off-target analysis)

data <- paste0(system.file(package = "CaSilico"),"/data/3D.fasta")

genome_dir <- paste0(system.file(package = "CaSilico"),"/data/genomes_off_target")

Local_Fasta <- file.path(genome_dir,list.files(genome_dir))

Local_Name <- c("Genome1","Genome2","Genome3")

CaSilico(ResultsFolder="Example2",
         TargetFasta=data,
         TargetAccession=NULL,
         CrisprTypes=c("casVI_A","casVI_B","casVI_D","casV_A","casV_B","casV_F1"),
         ConservationMethod = 1,
         ConservationThreshold=0.98,
         OffTarget = T,
         OffAsk = F,
         Organism=NULL,
         LocalOff=T,
         LocalFasta=Local_Fasta,
         LocalName=Local_Name)



#Run with TargetFasta and off-target analysis (online off-target analysis)

data <- paste0(system.file(package = "CaSilico"),"/data/3D.fasta")

CaSilico(ResultsFolder="Example2",
         TargetFasta=data ,
         TargetAccession=NULL,
         CrisprTypes=c("casVI_A"),
         ConservationMethod = 1,
         ConservationThreshold=0.98,
         OffTarget = T,
         OffAsk = F,
         Organism="Homo sapiens",
         LocalOff=F)
         
         
 #Run with accession numbers

 CaSilico(ResultsFolder="Example3",
          TargetAccession=c("U15717", "U15718"),
          CrisprTypes=c("casVI_A","casVI_B"),
          ConservationMethod = 1,
          ConservationThreshold=0.98,
          OffTarget = F,
          OffAsk = F)




 #Run with sequence coordinate
 
 CaSilico(ResultsFolder="Example4",
          TargetCoordinate=list(chromosome=c("1","1"),
                             start=c("2000","2000"),
                             end=c("2150","2150"),
                             strand=c("-","-"),
                             species = c("bos_taurus","bos_taurus")),
          CrisprTypes=c("casVI_B","casVI_A"),
          ConservationMethod = 1,
          ConservationThreshold=0.98,
          OffTarget = F,
          OffAsk = F)
          
          
      



```

# License
GPL-3
