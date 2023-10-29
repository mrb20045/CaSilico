#' In silico crRNA design for CRISPR-Cas12 (A, B and F1 subtypes) and CRISPR-Cas13 (A, B and D subtypes).
#'
#'
#' @param ResultsFolder   Name of the folder to save the final results; Default: ReaultsFolder="CaSilico_output".
#' @param TargetFasta   Full path to the target sequence(s) to design crRNA for; Example: TargetFasta="/usr/mrb/casilico/genome.fasta"
#' @param TargetAccession   Genome/gene/transcript accession numbers(s)(RefSeq or ENSEMBL ID)for the target sequence(s) to design crRNA for; Example: TargetAccession=c("U15717", "U15718").
#' @param TargetCoordinate   Genomic coordinate(s)for the target sequence(s) to design crRNA for. CaSilico automatically retrieves the sequences based on the coordinates and their species; Example: TargetCoordinate=list(chromosome=c("1","5"), start=c("2000","12000"), end=c("2150","15500"), strand=c("-","-"), species = c("bos_taurus","bos_taurus")).
#' @param CrisprTypes   A vector of CRSPR subtypes to design crRNAs based on their characteristics. One or more subtypes can be selected; Example: CrisprTypes=c("casVI_A","casVI_B","casVI_D","casV_A","casV_B","casV_F1").
#' @param ConservationMethod   Desired method (1 or 2) to detect the conserved regions when more than one target sequence is submitted. See detailes for more information; Default: 1.
#' @param ConservationThreshold   Desired threshold (between 0 and 1) to detect the conserved regions when more than one target sequence is submitted; Default: 0.95.
#' @param OffTarget   If TRUE, CaSilico search for off-targets in the sequences specified by the Organism or LocalOff argument. This analysis may take a long time based on the number of subtypes in CrisprTypes argument and number of the selected sequences to be searched (Organism or LocalFasta arguments); Default: FALSE.
#' @param OffAsk   If TRUE, in online mode before starting the off target analysis for each subtype, CaSilico ask from the user if he/she wants to continue or not; Default: FALSE.
#' @param Organism   A vector of desired organisms for online off target analysis. CaSilico automatically search the identified crRNAs on the sequences(s)(NR or Refseq database, based on the CRISPR-Cas subtype) of the specified organismes; Example: Organism=c("Homo sapiens","Bos taurus",...).
#' @param LocalOff   If TRUE, CaSilico search for off targets in local mode and Organism argument dose not works. LocalFile have to be set. Default: FALSE.
#' @param LocalFasta   Full path to the fasta file sequence(s) to be used for local off target analysis. One or more fasta file can be addressed; Example: LocalFasta=c("/yourpath1/filename1.fasta","/yourpath2/filename2.fasta").
#' @param LocalName   A vector including desired names related to the LocalFasta argument to determine the name of the results. If it is'nt specified by the user, Casilico considers sequence1, sequence2 and etc based on the number of fasta files in LocalFasta argument; Example: LocalName=c("Name1", "Name2").
#' @param Threads   Number of threads to use, default value is 4.; Example: Threads=8..
#' @return Casilico presents the results in HTML and PDF files that can be found in ReaultsFolder.
#' @section References:
#' ...
#' @section Further details:
#' Casilico accepts three kinds of input formats as target sequence to be scanned for all potential target sites: 1) direct sequence(s) (TargetFasta argument); 2) Accession nember (TargetAccession argument) and 3) genomic coordinate (TargetCoordinate argument). One of these arguments have to be determined by the user as TargetFasta is recommended. Also, there are two different approaches to identify conserved regions (using ConservationMethod argument): 1) Most or all positions in the fragments (with length of w) have to be conserved and ≤2 positions can be polymorphic (based on the defined conservation score), as these positions can be occurred in the specific locations of the window, 2) Most or all positions in the fragments (with length of w) have to be conserved and ≤2 positions can be polymorphic (based on the defined conservation score) in everywhere of the fragment (with length of w).
#' @export
#' @import  seqinr DT RRNA ape htmltools httr tibble dplyr readr
#' @examples
#
#'
#'############# (1) #############
#'# Run with TargetFasta and without off-target analysis
#'
#'# Full path of target sequences
#'data <- paste0(system.file(package = "CaSilico"),"/data/3D.fasta")
#'
#'# Run CaSilico for subtype VI_A with ConservationThreshold = 0.98 and ConservationMethod = 1
#'CaSilico(ResultsFolder="Example",
#'        TargetFasta=data,
#'         TargetAccession=NULL,
#'         CrisprTypes=c("casVI_A"),
#'         ConservationMethod = 1,
#'         ConservationThreshold=0.98,
#'         OffTarget = F,
#'         OffAsk = F)
#'##############################
#'#
#'#
#'#
#'############# (2) #############
#'
#'# Run with TargetFasta and offline off-target analysis
#'
#'# Full path of target sequences
#'data <- paste0(system.file(package = "CaSilico"),"/data/3D.fasta")
#'
#'# Full path to fasta file sequences for local off target analysis
#'genome_dir <- paste0(system.file(package = "CaSilico"),"/data/genomes_off_target")
#'Local_Fasta <- file.path(genome_dir,list.files(genome_dir))
#'
#'
#'# Names related to the LocalFasta argument to determine the name of the results. (optional)
#'Local_Name <- c("Genome1","Genome2","Genome3")
#'
#'# Run CaSilico for subtypes V_A and V_B with ConservationThreshold = 0.98 and ConservationMethod = 1
#'CaSilico(ResultsFolder="Example1",
#'         TargetFasta=data,
#'         TargetAccession=NULL,
#'         CrisprTypes=c("casV_A","casV_B"),
#'         ConservationMethod = 1,
#'         ConservationThreshold=0.98,
#'         OffTarget = T,
#'         OffAsk = F,
#'         Organism=NULL,
#'         LocalOff=T,
#'         LocalFasta=Local_Fasta,
#'         LocalName=Local_Name)
#'##############################
#'#
#'#
#'#
#'############# (3) #############
#'# Run with TargetFasta and online off-target analysis
#'
#'# Full path of target sequences
#'data <- paste0(system.file(package = "CaSilico"),"/data/3D.fasta")
#'
#'# Run CaSilico for subtype VI_A with ConservationThreshold = 0.98 and ConservationMethod = 1
#'# and online off-target analysis for Homo sapiens
#'CaSilico(ResultsFolder="Example2",
#'         TargetFasta=data ,
#'         TargetAccession=NULL,
#'         CrisprTypes=c("casVI_A"),
#'         ConservationMethod = 1,
#'         ConservationThreshold=0.98,
#'         OffTarget = T,
#'         OffAsk = F,
#'         Organism="Homo sapiens",
#'         LocalOff=F)
#'##############################
#'#
#'#
#'#
#'############# (4) #############
#'# Run with accession numbers
#'
#'# Run CaSilico for subtypes VI_A and VI_B with ConservationThreshold = 0.98 and ConservationMethod = 1
#'# The accession numbers used are: "U15717" and "U15718"
#' CaSilico(ResultsFolder="Example3",
#'          TargetAccession=c("U15717", "U15718"),
#'          CrisprTypes=c("casVI_A","casVI_B"),
#'          ConservationMethod = 1,
#'          ConservationThreshold=0.98,
#'          OffTarget = F,
#'          OffAsk = F)
#'##############################
#'#
#'#
#'#
#'############# (5) #############
#' # Run with sequence coordinate
#' 
#' # Run CaSilico for subtypes VI_A and VI_B with ConservationThreshold = 0.98 and ConservationMethod = 1
#' # The desired coordinates are: 
#' # (1) chromosome = 1, start = 2000, end = 2150, strand = "-", species = "bos_taurus"
#' # (2) chromosome = 1, start = 2000, end = 2150, strand = "-", species = "bos_taurus"
#' CaSilico(ResultsFolder="Example4",
#'          TargetCoordinate=list(chromosome=c("1","1"),
#'                             start=c("2000","2000"),
#'                             end=c("2150","2150"),
#'                             strand=c("-","-"),
#'                             species = c("bos_taurus","bos_taurus")),
#'          CrisprTypes=c("casVI_B","casVI_A"),
#'          ConservationMethod = 1,
#'          ConservationThreshold=0.98,
#'          OffTarget = F,
#'          OffAsk = F)
#'##############################
#'



CaSilico=function(ResultsFolder="CaSilico_output",
                  TargetFasta=NULL,
                  TargetAccession=NULL,
                  TargetCoordinate=list(chromosome=NA, start=NA, end=NA, strand=NA, species = NA ),
                  CrisprTypes,
                  ConservationMethod = 1,
                  ConservationThreshold=0.95,
                  OffTarget = F,
                  OffAsk=F,
                  Organism=NULL,
                  LocalOff=F,
                  LocalFasta=c("Path1","Path2"),
                  LocalName=c("Organism1","Organism2"),
                  Threads=4){




  library("seqinr")
  library("DT")
  library("RRNA")
  library("ape")
  library("htmltools")
  library("tibble")
  library("dplyr")
  library("readr")

  defaultW <- getOption("warn")
  options(warn = -1)

  host_system=Sys.info()[['sysname']]

  if (host_system %in% c("Windows","Linux")==F) {
    stop("Currently CaSilico only supports Windows and Linux systems")
  }
  first_main_dir=getwd()
  package_location<-system.file(package = "CaSilico")

  if (host_system=="Windows") {
    MAFFT_location=paste0(package_location,"/","dependency_win/mafft-win")
    RNAfold_path = paste0(package_location,"/","dependency_win/Blastn/RNAfold.exe")
    Blast_location=paste0(package_location,"/","dependency_win/Blastn/blastn.exe")
    makeblastdb_location=paste0(package_location,"/","dependency_win/Blastn/makeblastdb.exe")
    RNAplot_location=paste0(package_location,"/","dependency_win/Blastn/RNAplot.exe")
  }


  ############## change mod sotware ###########


  if (host_system=="Linux") {

    MAFFT_location=paste0(package_location,"/dependency_linux/mafft-linux64")
    RNAfold_path = paste0(package_location,"/dependency_linux/Blastn/RNAfold")
    Blast_location=paste0(package_location,"/dependency_linux/Blastn/blastn")
    makeblastdb_location=paste0(package_location,"/dependency_linux/Blastn/makeblastdb")
    RNAplot_location=paste0(package_location,"/dependency_linux/Blastn/RNAfold")

    path_software=c(RNAfold_path,
                    Blast_location,
                    makeblastdb_location,
                    RNAplot_location,
                    paste0(MAFFT_location,"/mafft.bat",collapse = ""))
    for (path in path_software) {
      chmod_command=paste0("chmod  777  ","'",path,"'")
      lapply(chmod_command, system)
      }
    lapply(paste0("chmod -R 777  ",MAFFT_location), system)
     }
  #############################################
  is_empty=function(X){
    length(X)==0
  }

  # input seq is path
  if (is_empty(TargetFasta)==FALSE) {
    fasta_file=read.fasta(TargetFasta)
    if (exists("fasta_file")==F | length(fasta_file)==0){
      stop("There is a problem in your file or path of file. please check it.")
    }
    TargetAccession=NULL
    TargetCoordinate=NULL
  }


  # input seq is accession number
  if (length(TargetAccession)>=1) {
    input_genes_fasta=read.GenBank(TargetAccession)
    if (is_empty(input_genes_fasta)==T) {
      stop("The access number entered is incorrect. please check it")
    }

    write.dna(input_genes_fasta,
              file ="Input_Genes.fasta",
              format = "fasta",
              append = FALSE,
              nbcol = 6,
              colsep = "",
              colw = 10)

    fasta_file=read.fasta("Input_Genes.fasta")
    file.remove("Input_Genes.fasta")
    TargetFasta=NULL
    TargetCoordinate=NULL
  }


  # input seq is coordinate
  if (all(is.na(TargetCoordinate))==F) {
    library("httr")
    getseq=function(chromosome, start, end, strand, species = "Homo_sapiens"){
      url = paste0("https://useast.ensembl.org/", species, "/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;output=fasta;r=",
                   chromosome, ":", start, "-", end, ";strand=", strand,
                   ";utr5=yes;cdna=yes;intron=yes;utr3=yes;peptide=yes;coding=yes;genomic=unmasked;exon=yes;_format=Text")
      r=GET(url)
      w=content(r, "text")
      s=strsplit(w,split = "\r\n")[[1]]
      paste0(s[-c(1,length(s))],collapse = "")
    }

    seqs_of_cordinate=c()
    for (no_cor in 1:length(TargetCoordinate$chromosome)) {
      seqs_of_cordinate=c(seqs_of_cordinate,
                          paste0(TargetCoordinate$chromosome[no_cor],"@",
                                 TargetCoordinate$start[no_cor],"@",
                                 TargetCoordinate$end[no_cor],"@",
                                 TargetCoordinate$strand[no_cor],"@",
                                 TargetCoordinate$species[no_cor],"@"))
    }

    input_coordinate_seq=c()

    for (cordinate_seq in 1:length(seqs_of_cordinate)) {
      temp_seq=seqs_of_cordinate[cordinate_seq]
      temp_seq=strsplit(temp_seq,split = "@")[[1]]
      temp_sequence="ddd"
      iter=0
      while(all(iter<=150 , temp_sequence=="ddd")){
        try({
          temp_sequence=getseq(temp_seq[1],temp_seq[2],temp_seq[3],temp_seq[4],temp_seq[5])

        },silent = T)
        iter=iter+1
      }
      input_coordinate_seq=c(input_coordinate_seq,temp_sequence)


    }

    if (any(input_coordinate_seq=="ddd")) {
      stop("CaSilico is unable to download your coordinate sequences due to a problem on the Ensembl site.
            Please try again later.")
    }




    write_line_in=c()
    for (i in 1:length(input_coordinate_seq)) {
      write_line_in=c(write_line_in,paste0(">",i),input_coordinate_seq[i])
    }

    file_Conaction<-file("Input_coordinate.fasta")
    writeLines(write_line_in, file_Conaction)

    close(file_Conaction)

    fasta_file=read.fasta("Input_coordinate.fasta")
    file.remove("Input_coordinate.fasta")
    TargetFasta=NULL
    TargetCoordinate=NULL
  }




  final_file_list=list()
  # calculated total nucleotide for mafft limitation (mafft limitation > 1500 sequences * 30kb)
  total_Nucleotide=0
  for (h in 1:length(fasta_file)) {
    total_Nucleotide=total_Nucleotide+length(fasta_file[[h]])
    empty_pos=which(fasta_file[[h]]==" ")
    if (length(empty_pos)>0) {
      final_file_list[[h]]=fasta_file[[h]][-empty_pos]
    }
    else {final_file_list[[h]]=fasta_file[[h]]}

  }

  if ( total_Nucleotide > 40000000 ) {
    print("Restrictions on Multiple Alignment")
  }

  if (total_Nucleotide <= 40000000) {
    if (length(final_file_list)==1) {
      write.fasta(list(final_file_list[[1]],final_file_list[[1]]), c(names(fasta_file),names(fasta_file)), "final_fasta_file.fasta") 
      }else{
      write.fasta(final_file_list, names(fasta_file), "final_fasta_file.fasta")
      }
    setwd(MAFFT_location)
    if (file.exists("final_fasta_file.fasta")) {
      file.remove("final_fasta_file.fasta")
    }
    if (file.exists("final_fasta_file.aln")) {
      file.remove("final_fasta_file.aln")
    }
    setwd(first_main_dir)
    file.copy("final_fasta_file.fasta",MAFFT_location)
    file.remove("final_fasta_file.fasta")
  }





  # Run mafft
  setwd(MAFFT_location)
  if (host_system=="Linux") {
    MAFFT_list<- paste0("./mafft.bat  --auto ","   --thread ",Threads, " --inputorder final_fasta_file.fasta > final_fasta_file.aln")

  }
  if (host_system=="Windows") {
    MAFFT_list<- paste0("mafft  --auto ","   --thread ",Threads, " --inputorder final_fasta_file.fasta > final_fasta_file.aln")
  }
  lapply(MAFFT_list, system)   ###run MAFFT for all input .txt files
  FMDValn <- read.fasta(file = "final_fasta_file.aln", as.string = T)
  file.copy("final_fasta_file.aln",first_main_dir)
  file.remove("final_fasta_file.fasta")
  file.remove("final_fasta_file.aln")
  setwd(first_main_dir)

  topoTypes_sequences <- FMDValn
  topoTypes_names <- names(FMDValn)

  ################# calculate frequency matrix for each base in sequences #########
  score_matrix=matrix(0,5,length(strsplit(topoTypes_sequences[[1]],split = "")[[1]]))
  rownames(score_matrix)=c("g","t","c","a","-")
  for (i in 1:length(topoTypes_sequences)) {
    seq=strsplit(topoTypes_sequences[[i]],split = "")[[1]]
    for (j in 1:length(seq)) {
      # one Nucleotide symbols
      if(seq[j]=="g" | seq[j]=="G"){
        score_matrix[1,j]= score_matrix[1,j]+1}
      else if(seq[j]=="t" | seq[j]=="T" | seq[j]=="U" | seq[j]=="u"){
        score_matrix[2,j]= score_matrix[2,j]+1}
      else if(seq[j]=="c" | seq[j]=="C"){
        score_matrix[3,j]= score_matrix[3,j]+1}
      else if(seq[j]=="a" |seq[j]=="A"){
        score_matrix[4,j]= score_matrix[4,j]+1}
      else if(seq[j]=="-"){
        score_matrix[5,j]= score_matrix[5,j]+1}

      # four Nucleotide symbols
      else if(seq[j]=="n" | seq[j]=="N"  ){
        score_matrix[1:4,j]= score_matrix[1:4,j]+0.25}

      # two Nucleotide symbols
      else if(seq[j]=="r" | seq[j]=="R"  ){
        score_matrix[c(1,4),j]= score_matrix[c(1,4),j]+0.5}
      else if(seq[j]=="y" | seq[j]=="Y"  ){
        score_matrix[c(2,3),j]= score_matrix[c(2,3),j]+0.5}
      else if(seq[j]=="k" | seq[j]=="K"  ){
        score_matrix[1:2,j]= score_matrix[1:2,j]+0.5}
      else if(seq[j]=="m" | seq[j]=="M"  ){
        score_matrix[3:4,j]= score_matrix[3:4,j]+0.5}
      else if(seq[j]=="s" | seq[j]=="S"  ){
        score_matrix[c(1,3),j]= score_matrix[c(1,3),j]+0.5}
      else if(seq[j]=="w" | seq[j]=="W"  ){
        score_matrix[c(2,4),j]= score_matrix[c(2,4),j]+0.5}

      # three Nucleotide symbols
      else if(seq[j]=="b" | seq[j]=="B"  ){
        score_matrix[1:2,j]= score_matrix[1:2,j]+(0.33)
        score_matrix[3,j]= score_matrix[3,j]+(0.34)}
      else if(seq[j]=="d" | seq[j]=="D"  ){
        score_matrix[1:2,j]= score_matrix[1:2,j]+(0.33)
        score_matrix[4,j]= score_matrix[4,j]+(0.34)}
      else if(seq[j]=="h" | seq[j]=="H"  ){
        score_matrix[2:3,j]= score_matrix[2:3,j]+(0.33)
        score_matrix[4,j]= score_matrix[4,j]+(0.34)}
      else if(seq[j]=="v" | seq[j]=="V"  ){
        score_matrix[3:4,j]= score_matrix[3:4,j]+(0.33)
        score_matrix[1,j]= score_matrix[1,j]+(0.34)}
    }
  }
  #################################################################################




  ################### score matrix normalization ####################
  normal_score_matrix=score_matrix/length(topoTypes_sequences)
  normal_score_matrix1=normal_score_matrix
  row_name_matrix=c("g","t","c","a","Del")
  rownames(normal_score_matrix1)=row_name_matrix
  ####################################################################





  ######################### creATE msa sequence #####################
  extract_higher_frequency=function(x){
    for (g in 1:nrow(normal_score_matrix)) {
      if(normal_score_matrix[g,x]==max(normal_score_matrix[,x]))
        break
    }
    normal_score_matrix[g,x]
  }

  final_score=c()
  for (k in 1:ncol(normal_score_matrix)) {
    final_score=c(final_score,extract_higher_frequency(k))
  }

  final_str=toString(names(final_score),)

  scores=as.data.frame(final_score)[[1]]

  str_frags=strsplit(final_str,split = ", ")[[1]]

  final_str_frags=c()
  final_score1=c()
  pso_del=c()
  for (i in 1:length(str_frags)) {
    if (str_frags[i]!="-") {
      final_str_frags=c(final_str_frags,str_frags[i])
      final_score1=c(final_score1,scores[i])
    }
    else if (str_frags[i]=="-") {
      pso_del=c(pso_del,i)
    }
  }

  if (length(pso_del)==0) {
    original_pos_in_matrix=1:length(str_frags)
  }
  if (length(pso_del)!=0) {
    original_pos_in_matrix=setdiff(1:length(str_frags),pso_del)
  }

  ###############################################################








  ########################### Functions ########################

  {

    ########### mismatch delete function####


    if (ConservationMethod==1) {
      mismatch_del=function(score,start,end,seed_start,seed_end){
        del=c()
        mismatch_locat=c()
        for (i in 1:length(score)) {
          if (score[i]<ConservationThreshold) {
            mismatch_locat=c(mismatch_locat,i)
          }}
        if (length(mismatch_locat)>2) {
          del=TRUE
        }
        if (length(mismatch_locat)==2) {
          if (any(mismatch_locat%in%start:end)==TRUE) {
            del=TRUE}
        }
        if (length(mismatch_locat)==2) {
          if (any(mismatch_locat%in%start:end)!=TRUE) {
            del=FALSE}
        }
        if (length(mismatch_locat)==1) {
          if (mismatch_locat%in%seed_start:seed_end==TRUE) {
            del=TRUE}
        }
        if (length(mismatch_locat)==1) {
          if (mismatch_locat%in%seed_start:seed_end!=TRUE) {
            del=FALSE}
        }
        if (length(mismatch_locat)==0) {
          del=FALSE
        }
        del
      }
    }



    if (ConservationMethod==2) {
      mismatch_del=function(score,start,end,seed_start,seed_end){
        del=c()

        if (start<=end & seed_start<=seed_end)
        { if (length(which(score<ConservationThreshold))<3) {
          del=FALSE
        }
          else if (length(which(score<ConservationThreshold))>=3) {
            del=TRUE
          }}
        del
      }
    }
    ###############################






    ############# revers compliment function########
    revers_compliment=function (seq)
    {
      compliment = function(n, seq_type) {
        if (seq_type == "D") {
          if (n == "a") {
            return("t")
          }
          if (n == "A") {
            return("T")
          }
        }
        if (seq_type == "R") {
          if (n == "a") {
            return("u")
          }
          if (n == "A") {
            return("U")
          }
        }
        if (n == "t") {
          return("a")
        }
        if (n == "c") {
          return("g")
        }
        if (n == "g") {
          return("c")
        }
        if (n == "C") {
          return("G")
        }
        if (n == "G") {
          return("C")
        }
        if (n == "T") {
          return("A")
        }
        if (n == "u") {
          return("a")
        }
        if (n == "U") {
          return("A")
        }
        if (n != "U" & n != "A" & n != "T" &
            n != "G" & n != "C" & n != "u" &
            n != "a" & n != "t" & n != "g" &
            n != "c") {
          return(n)
        }
      }
      Rcompliment_sequence = c()
      if (length(seq) > 1) {
        seq = seq
      }
      if (length(seq) == 1) {
        seq = strsplit(seq, split = "")[[1]]
      }
      type = "D"
      if (any(seq == "u") | any(seq == "U")) {
        type = "R"
      }

      for (i in length(seq):1) {
        Rcompliment_sequence = c(Rcompliment_sequence, compliment(seq[i],
                                                                  seq_type = type))
      }
      Rcompliment_sequence
    }





    # compliment  function
    compliment_seq=function (seq)
    {
      compliment = function(n, seq_type) {
        if (seq_type == "D") {
          if (n == "a") {
            return("t")
          }
          if (n == "A") {
            return("T")
          }
        }
        if (seq_type == "R") {
          if (n == "a") {
            return("u")
          }
          if (n == "A") {
            return("U")
          }
        }
        if (n == "t") {
          return("a")
        }
        if (n == "c") {
          return("g")
        }
        if (n == "g") {
          return("c")
        }
        if (n == "C") {
          return("G")
        }
        if (n == "G") {
          return("C")
        }
        if (n == "T") {
          return("A")
        }
        if (n == "u") {
          return("a")
        }
        if (n == "U") {
          return("A")
        }
        if (n != "U" & n != "A" & n != "T" &
            n != "G" & n != "C" & n != "u" &
            n != "a" & n != "t" & n != "g" &
            n != "c") {
          return(n)
        }
      }
      compliment_sequence = c()
      if (length(seq) > 1) {
        seq = seq
      }
      if (length(seq) == 1) {
        seq = strsplit(seq, split = "")[[1]]
      }
      type = "D"
      if (any(seq == "u") | any(seq == "U")) {
        type = "R"
      }

      for (i in 1:length(seq)) {
        compliment_sequence = c(compliment_sequence, compliment(seq[i],
                                                                seq_type = type))
      }
      compliment_sequence
    }





    # revers function
    revers_function=function(v){
      revers_v=c()
      for (i in length(v):1) {
        revers_v=c(revers_v,v[i])
      }
      revers_v
    }

    ##########################










    ################## check internet function ###########
    if (host_system=="Windows") {
      havingIP <- function() {
        if (.Platform$OS.type == "windows") {
          ipmessage <- system("ipconfig", intern = TRUE)
        } else {
          ipmessage <- system("ifconfig", intern = TRUE)
        }
        validIP <- "((25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)[.]){3}(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)"
        any(grep(validIP, ipmessage))
      }
    }
   

    if (host_system=="Linux") {
      havingIP <- function() {
        if (.Platform$OS.type == "unix") {
          ipmessage <- system("ip r", intern = TRUE)
        } else {
          ipmessage <- system("ip r", intern = TRUE)
        }
        validIP <- "((25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)[.]){3}(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)"
        any(grep(validIP, ipmessage))
      }
    }
    
 ##################################







    ################ crRNA plot function for RNA sequence ##############
    crRNA_plot=function(sequence,crstart,start,pos_mutatiopns,
                        type_SNP,main_plot,W,u_rich_number,
                        acc_flank_number,eficieny_number,cr_acc) {
      sequence=toupper(sequence)
      plot(NULL, xlim =c(1,length(sequence)), ylim =c(1,12),type = "n",
           axes =F,main = main_plot ,ann = T,xlab = "sequence",ylab = "")

      axis(side=1, at = 1:length(sequence), labels= sequence,cex.axis = 0.7,
           font = 0.6,col="darkred",col.axis="darkred")

      axis(side=1, at = crstart:(crstart+W-1), labels= sequence[crstart:(crstart+W-1)] ,
           cex.axis = 0.7,col="darkblue", font = 1,col.axis="darkblue")

      end=start+W
      x=c(crstart,crstart,crstart+W-1,crstart+W-1)
      y=c(2,4,4,2)
      polygon(x,y,col = "pink3",border = "pink3")

      text((x[4]+x[1])/2,(y[2]+y[1])/2,labels = " CRISPR RNA",cex = 1.4)
      text(x[1],1,labels = start,cex = 0.8)
      text(x[4],1,labels = end,cex = 0.8)


      if (is.numeric(pos_mutatiopns)) {
        for (mut in 1:length(pos_mutatiopns)) {
          pos_mutatiopns1=pos_mutatiopns+crstart-1
          M=pos_mutatiopns1[mut]

          if (length(pos_mutatiopns)==1) {
            segments(M,4,M,6)
            points(M,6,pch = 16, cex = 1.7, col = "red")
            text((M),6.5,type_SNP[mut],cex = 1.2)
          }
          if (length(pos_mutatiopns)==2) {
            if (abs(pos_mutatiopns[1]-pos_mutatiopns[2])>1) {
              segments(M,4,M,6)
              points(M,6,pch = 16, cex = 1.7, col = "red")
              text((M),6.5,type_SNP[mut],cex = 1.2)
            }
            else if (abs(pos_mutatiopns[1]-pos_mutatiopns[2])==1) {
              segments(M,4,M,6)
              points(M,6,pch = 16, cex = 1.7, col = "red")
              text((M-0.1),6.5,type_SNP[mut],cex = 1.2,srt=90)
            }
          }
        }

      }


      segments(crstart,9,crstart+W-1,9,lwd=1.5,col = "darkblue")
      text(crstart+W-1,9+0.03,">",col = "darkblue")
      text(crstart,9+0.03,"<",col = "darkblue")
      text((x[1]+x[4])/2,9+2,paste(cr_acc,"Accessible",collapse = ""),cex = 0.9,col = "darkblue")



      if (crstart>10) {

        segments(1,(y[2]+y[1])/2,crstart,(y[2]+y[1])/2,col="darkred")
        text(1,(y[2]+y[1])/2+0.03,"<",col="darkred")
        text(x[1]/2,(y[2]+y[1])/2+1,u_rich_number,cex = 0.9,col="darkred")
        text(x[1]/2,(y[2]+y[1])/2+2,paste(acc_flank_number,"Accessible",collapse = ""),cex = 0.9,col="darkred")
      }


      if (length(sequence)-(crstart+W)>10){
        segments(crstart+W-1,(y[2]+y[1])/2,length(sequence),(y[2]+y[1])/2,col="darkred")
        text(length(sequence),(y[2]+y[1])/2+0.03,">",col="darkred")
        text((length(sequence)+x[4])/2,(y[2]+y[1])/2+1,u_rich_number,cex = 0.9,col="darkred")
        text((length(sequence)+x[4])/2,(y[2]+y[1])/2+2,paste(acc_flank_number,"Accessible",collapse = ""),cex = 0.9,col="darkred")
      }


    }
    ###################################################






    ################# crRNA plot function for DNA sequence  ##################
    crRNA_plot_DNA=function(sequence,crstart,start,pos_mutatiopns,
                            type_SNP,main_plot,W,eficieny_number,
                            sign_strand,length_pam) {

      sequence=toupper(sequence)
      end=start+W

      plot(NULL, xlim =c(1,length(sequence)), ylim =c(1,6),type = "n",
           axes =F,main = main_plot ,ann = T,xlab = "sequence",ylab = "")

      axis(side=1, at = 1:length(sequence), labels= sequence,
           cex.axis = 0.7, font = 0.6,col="darkred",col.axis="darkred")

      x=c(crstart,crstart,crstart+W-1,crstart+W-1)
      y=c(2,3,3,2)
      polygon(x,y,col = "pink3",border = "pink3")




      text((x[4]+x[1])/2,(y[2]+y[1])/2,labels = " CRISPR RNA",cex = 1.4)
      segments(1,1.3,length(sequence),1.3,col="darkred")
      text(1:length(sequence),1.38,"--",srt = 90,col="darkred")
      text(1:length(sequence),1.8,toupper(compliment_seq(tolower(sequence))),col="darkred",cex = 0.7)



      if (sign_strand=="+") {
        axis(side=1, at = crstart:(crstart+W-1), labels= sequence[crstart:(crstart+W-1)] ,
             cex.axis = 0.7,col="darkblue", font = 1,col.axis="darkblue")


        segments((crstart-length_pam),2.1,(crstart-1),2.1,col="black")
        text(((crstart-length_pam)+(crstart-1))/2,2.2,"<",srt = 90,col="black")
        text(((crstart-length_pam)+(crstart-1))/2,2.45,"------",srt = 90,col="black")
        text(((crstart-length_pam)+(crstart-1))/2,3,"PAM",col="black")
        text(x[1],1.05,labels = start,cex = 1)
        text(x[4],1.05,labels = end,cex = 1)



        segments((crstart-length_pam),1.3,(crstart-1),1.3,col="orange")
        text((crstart-length_pam):(crstart-1),1.38,"--",srt = 90,col="orange")
        text((crstart-length_pam):(crstart-1),1.8,
             toupper(compliment_seq(tolower(sequence)))[(crstart-length_pam):(crstart-1)],
             col="orange",cex = 0.7)

      }

      if (sign_strand=="-") {
        segments(crstart,1.3,(crstart+W-1),1.3,col="darkblue")
        text(crstart:(crstart+W-1),1.38,"--",srt = 90,col="darkblue")
        text(crstart:(crstart+W-1),1.8,
             toupper(compliment_seq(tolower(sequence)))[crstart:(crstart+W-1)],
             col="darkblue",cex = 0.7)


        axis(side=1, at = (crstart+W):(crstart+W+length_pam-1), labels= sequence[(crstart+W):(crstart+W+length_pam-1)] ,
             cex.axis = 0.7,col="orange", font = 1,col.axis="orange")


        segments((crstart+W),2.1,(crstart+W+length_pam-1),2.1,col="black")
        text(((crstart+W)+(crstart+W+length_pam-1))/2,2.2,"<",srt = 90,col="black")
        text(((crstart+W)+(crstart+W+length_pam-1))/2,2.45,"------",srt = 90,col="black")
        text(((crstart+W)+(crstart+W+length_pam-1))/2,3,"PAM",cex = 1.1,col="black")
        text(x[1],1.05,labels =end ,cex = 1)
        text(x[4],1.05,labels = start,cex = 1)

      }



      if (is.numeric(pos_mutatiopns)) {
        for (mut in 1:length(pos_mutatiopns)) {
          pos_mutatiopns1=pos_mutatiopns+crstart-1
          M=pos_mutatiopns1[mut]
          segments(M,3,M,4)
          points(M,4,pch = 16, cex = 1.7, col = "red")

          text((M-0.1),4.5,type_SNP[mut],cex = 1.1,srt = 90)
        }

      }


    }

    ####################################################






    ############### my Dot plot function ################
    my_Dot_Plot=function(sequence1,sequence2, wsize = 1 , wstep = 1 , nmatch = 1 ,
                         table_col = "goldenrod1",points_col = "red", xlab = "sequence1", ylab = "sequence2")
    {
      #Split the Elements of strin
      if (length(sequence1)==1){
        x <- strsplit(sequence1,split = "")[[1]]}
      if (length(sequence2)==1){
        y <- strsplit(sequence2,split = "")[[1]]}
      if (length(sequence1)!=1){
        x <-sequence1}
      if (length(sequence2)!=1){
        y <- sequence2}





      if (nchar(x[1]) > 1)
        stop("seq1 should be provided as a vector of single chars")
      if (nchar(y[1]) > 1)
        stop("seq2 should be provided as a vector of single chars")
      if (wsize < 1)
        stop("non allowed value for wsize")
      if (wstep < 1)
        stop("non allowed value for wstep")
      if (nmatch < 1)
        stop("non allowed value for nmatch")
      if (nmatch > wsize)
        stop("nmatch > wsize is not allowed")
      mkwin <- function(seq, wsize, wstep) {
        sapply(seq(from = 1, to = length(seq) - wsize + 1, by = wstep),
               function(i) c2s(seq[i:(i + wsize - 1)]))
      }
      wseq1 <- mkwin(x, wsize, wstep)
      wseq2 <- mkwin(y, wsize, wstep)
      if (nmatch == wsize) {
        xy <- outer(wseq1, wseq2, "==")
      }
      else {
        "%==%" <- function(x, y) colSums(sapply(x, s2c) == sapply(y,
                                                                  s2c)) >= nmatch
        xy <- outer(wseq1, wseq2, "%==%")
      }

      x=toupper(x)
      y=toupper(y)
      a=matrix(0,length(x),length(y))
      # Naming the rows and columns of the formed matrix
      row.names(a)=x
      colnames(a)=y


      for (i in 1:nrow(xy)) {
        for (j in 1:ncol(xy)){
          if (xy[i,j]==T) {
            a[i,j]=1
          }
        }
      }



      h_line=seq(0+((seq(0,max(length(x),length(y))-1,
                         l=length(x))[2]-seq(0,max(length(x),length(y))-1,
                                             l=length(x))[1])/2),max(length(x),length(y)),
                 by=seq(0,max(length(x),length(y))-1,
                        l=length(x))[2]-seq(0,max(length(x),length(y))-1,
                                            l=length(x))[1])
      v_line=seq(1+((seq(1,max(length(x),length(y)),
                         l=length(y))[2]-seq(1,max(length(x),length(y)),
                                             l=length(y))[1])/2),max(length(x),length(y)),
                 by=seq(1,max(length(x),length(y)),
                        l=length(y))[2]-seq(1,max(length(x),length(y)),
                                            l=length(y))[1])

      aaaa=0.04*max(length(x),length(y))


      plot(0:max(length(x),length(y)),0:max(length(x),length(y)),type = "n",
           axes =F,frame.plot=TRUE,main ="" ,
           ann = T,xlab = xlab , ylab = ylab)



      aa=c(-aaaa,-aaaa,max(length(x),length(y))+aaaa,max(length(x),length(y))+aaaa)
      bb=c(h_line[length(h_line)],
           max(length(x),length(y))+aaaa,
           max(length(x),length(y))+aaaa,
           h_line[length(h_line)])

      polygon(aa,bb,col = table_col)


      aa1=c(-aaaa,-aaaa,max(length(x),length(y))+aaaa,max(length(x),length(y))+aaaa)
      bb1=c(max(length(x),length(y))-aaaa,
            max(length(x),length(y))+0.5,
            max(length(x),length(y))+0.5,
            max(length(x),length(y))-aaaa)

      polygon(bb1-max(length(x),length(y)),aa1,col = table_col)



      # y axes
      axis(side=2, at = seq(0,max(length(x),length(y))-1,l=length(x)),
           labels= length(x):1,cex.axis = 0.7, font = 2)

      abline(h=h_line,col="black",lwd=1, lty=1)



      text((aa1[2]+h_line[1])/2,
           seq(0,max(length(x),length(y))-1,l=length(x)),
           rev(x),col = "black",cex = 1,pch = 16,font = 2)


      #x axes
      axis(side=3, at = seq(1,max(length(x),length(y)),
                            l=length(y)), labels= 1:length(y),cex.axis = 0.7, font = 2)


      abline(v=v_line, col="black",lwd=1, lty=1)


      abline(v=0.5,col="black",lwd=1, lty=1)

      text(seq(1,max(length(x),length(y)),l=length(y))
           ,(bb[2]+bb[1])/2,
           y,col = "black",cex = 1,pch = 16,font = 2)



      for (i in 1:length(x)){
        for (j in 1:length(y))
          if (a[i,j]==1)
            points(seq(1,max(length(x),length(y)),l=length(y))[j],seq(0,max(length(x),length(y))-1,l=length(x))[length(x)-i+1],cex=1.5,col=points_col,pch = 16)
      }
    }
    #####################################################







    ############### prediction structure of sequence ###############
    pridict_structure=function(seq,RNAfold_location)
    {
      seq_split=strsplit(seq,split = "")[[1]]
      write.fasta(seq_split,"seq_split","seq_split.fasta")
      structure_of_seq=system(paste(RNAfold_location,"  ","seq_split.fasta" ),intern=T)
      file.remove("seq_split.fasta")
      file.remove("seq_split_ss.ps")
      structure_of_seq_final=strsplit(structure_of_seq[3],split = "")[[1]][1:length(seq_split)]
      free_anergy_seq=as.numeric(paste(strsplit(structure_of_seq[3],split = "")[[1]][(length(seq_split)+3):(length(strsplit(structure_of_seq[3],split = "")[[1]])-1)],collapse = ""))
      list( "structure"=structure_of_seq_final, "free_anergy_seq"= free_anergy_seq)
    }
    ################################################################



    ############### secondary structure plot of backbone and spacer ######################
    sec_structure_plot=function(seq,sructure_seq,first_backbone,backbone_name,length_backbone,main,free_energy)
    {
      length_seq=length(strsplit(seq,split = "")[[1]])
      ct=makeCt(paste(sructure_seq,collapse = ""),paste(seq,collapse = ""))
      coord=ct2coord(ct)
      if (first_backbone==T) {
        RNAPlot(coord,labTF=T,pointSize=1.8,lineWd=0.5,modspec=TRUE,ranges=data.frame(min=c(1,length_backbone+1),
                                                                                      max=c(length_backbone,length_seq),
                                                                                      col=c(2,5),desc=c(backbone_name,"Spacer")),main=main)
        RNAPlot(coord,nt=TRUE,add=TRUE,dp=0.7,tsize=1)
      }
      else if (first_backbone==F) {
        RNAPlot(coord,labTF=T,pointSize=1.8,lineWd=0.5,modspec=TRUE,ranges=data.frame(min=c(1,length_seq-length_backbone+1),
                                                                                      max=c(length_seq-length_backbone,length_seq),
                                                                                      col=c(5,2),desc=c("Spacer",backbone_name)),main=main)
        RNAPlot(coord,nt=TRUE,add=TRUE,dp=0.7,tsize=1)
      }
      else if (first_backbone=="N") {
        RNAPlot(coord,labTF=T,pointSize=1.8,lineWd=0.5,modspec=TRUE,ranges=data.frame(min=1,
                                                                                      max=length_seq,
                                                                                      col=5,desc= paste0("Spacer ","(Free Energy = ",free_energy,")")),main=main)
        RNAPlot(coord,nt=TRUE,add=TRUE,dp=0.7,tsize=1)
      }
    }
    ######################################################################################


  }

  ##############################################################







  for (CRISPR_Type in CrisprTypes) {

    if (CRISPR_Type=="casVI_A") {
      W=28
      Nucleotide_rich="t"
      backbone="GATTTAGACTACCCCAAAAACGAAGGGGACTAAAAC"
      local_rich_type="Local_U_Rich"
      Normal_rich_type="Normalized_U_Rich"
      start=8
      end=27
      seed_start=13
      seed_end=24
      blast_database="refseq_rna"
      length_pam_cas=0
      first_backbone=T
      blast_strand="plus"
      ID_code="VI_A"
      bacteria="Leptotrichia wadei"
      cas_protein_type="Cas13a/C2c2"
      max_free_energy=-57
      direct_repeat="Direct Repeat"
    }

    if (CRISPR_Type=="casVI_D") {
      W=22
      Nucleotide_rich="t"
      backbone="TACCCCTACCAACTGGTCGGGGTTTGAAAC"
      local_rich_type="Local_U_Rich"
      Normal_rich_type="Normalized_U_Rich"
      start=3
      end=20
      seed_start=2
      seed_end=8
      blast_database="refseq_rna"
      length_pam_cas=0
      first_backbone=T
      blast_strand="plus"
      ID_code="VI_D"
      bacteria="Ruminococcus flavefaciens"
      cas_protein_type="Cas13d"
      max_free_energy=-45
      direct_repeat="Direct Repeat"
    }

    if (CRISPR_Type=="casVI_B") {
      W=30
      Nucleotide_rich="a"
      backbone="GTTGTGGAAGGTCCAGTTTTGAGGGGCTATTACAAC"
      local_rich_type="Local_A_Rich"
      Normal_rich_type="Normalized_A_Rich"
      start=12
      end=29
      seed_start=12
      seed_end=26
      blast_database="refseq_rna"
      length_pam_cas=2
      first_backbone=F
      blast_strand="plus"
      ID_code="VI_B"
      bacteria="Prevotella sp."
      cas_protein_type="Cas13b/C2c6"
      max_free_energy=-59
      direct_repeat="Direct Repeat"
    }

    if (CRISPR_Type=="casV_B") {
      W=20
      backbone="GTCTAGAGGACAGAATTTTTCAACGGGTGTGCCAATGGCCACTTTCCAGGTGGCAAAGCCCGTTGAGCTTCTCAAATCTGAGAAGTGGCAC"
      start=0
      end=0
      seed_start=0
      seed_end=0
      blast_database="nr"
      length_pam_cas=3
      first_backbone=T
      blast_strand="both"
      ID_code="V_B"
      bacteria="Alicyclobacillus acidiphilus"
      cas_protein_type="Cas12b/C2c1"
      max_free_energy=-102
      direct_repeat="Scaffold"
    }


    if (CRISPR_Type=="casV_A") {
      W=20
      backbone="TAATTTCTACTAAGTGTAGAT"
      start=1
      end=6
      seed_start=1
      seed_end=6
      blast_database="nr"
      length_pam_cas=4
      first_backbone=T
      blast_strand="both"
      ID_code="V_A"
      bacteria="Lachnospiraceae bacterium"
      cas_protein_type="Cas12a/Cpf1"
      max_free_energy=-34
      direct_repeat="Direct Repeat"
    }


    if (CRISPR_Type=="casV_F1") {
      W=20
      backbone="TGCAGAACCCGAATAGACGAATGAAGGAATGCAAC"
      start=9
      end=16
      seed_start=9
      seed_end=16
      blast_database="nr"
      length_pam_cas=0
      first_backbone=T
      blast_strand="both"
      ID_code="V_F1"
      bacteria="Uncultured archaeon"
      cas_protein_type="Cas14a/Cas12f1"
      max_free_energy=-34
      direct_repeat="Direct Repeat"
    }

    if (length(final_str_frags)<W+105) {
      message(paste("Your gene is too short, so casilico can not design  guides RNA for",CRISPR_Type))
    }
    if (length(final_str_frags)>=W+105) {

      backbone_structure_and_anergy=pridict_structure(backbone,RNAfold_location = RNAfold_path)
      backbone_structure=backbone_structure_and_anergy$structure
      backbone_free_enargy=backbone_structure_and_anergy$free_anergy_seq



      ############### create dir for final result ############

      setwd(first_main_dir)
      name_of_fasta_file_genes=ResultsFolder
      result_folder_name=paste0(CRISPR_Type,"-","ConservationMethod_",ConservationMethod,"-","threshold_",ConservationThreshold*100,"-result_",name_of_fasta_file_genes)
      plots_folder_name="Plots"
      sec_plots_folder_name="Secondary Structure"
      MSA_folder="Multiple Sequenc Alignment"
      first_dir=paste0(first_main_dir,"/",result_folder_name)

      if (dir.exists(result_folder_name)) {
        unlink(paste0("./",result_folder_name), recursive = TRUE)
      }

      dir.create(result_folder_name)

      file.copy("final_fasta_file.aln",first_dir)
      #file.remove("final_fasta_file.aln")
      setwd(first_dir)
      dir.create(plots_folder_name)
      if (CRISPR_Type=="casVI_D"|CRISPR_Type=="casVI_A"|CRISPR_Type=="casVI_B"| CRISPR_Type=="casV_F1") {
        dir.create(sec_plots_folder_name)
      }

      dir.create(MSA_folder)
      file.copy("final_fasta_file.aln",MSA_folder)
      file.remove("final_fasta_file.aln")

      setwd(MSA_folder)
      file.rename("final_fasta_file.aln","Multiple Sequence Alignment.aln")
      write.fasta(final_str_frags,"Consensus Sequence","consensus sequence.fasta")
      setwd(first_dir)
      #####################################################








      #################### first empty value ##############

      crRNA_seq=c()
      vv=list()
      normalized_U_rich=list()
      local_U_rich=list()
      start_pos=list()
      mismatch_pos=list()
      mismatch_num=list()
      efficiency=list()
      type_if_SNP=list()
      accessible_side=list()
      accessible_flank_side=list()
      rank=list()
      prtospaser=list()
      strand=list()
      crRNA_free_energy=list()
      Guide_RNA_ID=list()
      link_plots=list()
      GC_content=list()
      score_Gaide=list()

      #####################################################









      ###################### secondary structure ################
      if (CRISPR_Type=="casVI_D"|
          CRISPR_Type=="casVI_A"|
          CRISPR_Type=="casVI_B"|
          CRISPR_Type=="casV_F1"){
        if (CRISPR_Type=="casVI_D"|
            CRISPR_Type=="casVI_A"|
            CRISPR_Type=="casVI_B") {
          u_rich_in_whole=sum(final_score1[which(final_str_frags==Nucleotide_rich)])/sum(final_score1)
        }

        if (length(final_str_frags)<=2000){

          write.fasta(final_str_frags,"final_str_frags","final_str_frags.fasta")
          structure_of_RNAfold=system(paste(RNAfold_path,"  ","final_str_frags.fasta" ),intern=T)
          file.remove("final_str_frags.fasta")
          file.remove("final_str_frags_ss.ps")
          RNAfols_structure_format=strsplit(structure_of_RNAfold[3],split = "")[[1]][1:length(final_str_frags)]

          RNAplot_file_name=paste0(">","final_structure")
          fileConn<-file("sec_stracuter.txt")
          writeLines(c(RNAplot_file_name,
                       paste(final_str_frags,collapse = ""),
                       paste(RNAfols_structure_format,collapse = "")), fileConn)

          close(fileConn)

          if (host_system=="Windows") {
            RNAPlot_command_line=paste(RNAplot_location," -o  svg ", " sec_stracuter.txt ")
            lapply(RNAPlot_command_line, system)
            file.remove("sec_stracuter.txt")
            file.copy("final_structure_ss.svg",plots_folder_name)
            file.remove("final_structure_ss.svg")
          }
          if (host_system=="Linux") {
            RNAPlot_command_line=paste(RNAplot_location, " sec_stracuter.txt ")
            lapply(RNAPlot_command_line, system)
            file.remove("sec_stracuter.txt")
            file.copy("final_structure_ss.ps",plots_folder_name)
            file.remove("final_structure_ss.ps")
          }
        }


        if (length(final_str_frags)>2000) {

          # start
          write.fasta(final_str_frags[1:2000],"final_str_frags","final_str_frags.fasta")
          structure_of_RNAfold=system(paste(RNAfold_path,"  ","final_str_frags.fasta" ),intern=T)
          file.remove("final_str_frags.fasta")
          file.remove("final_str_frags_ss.ps")
          RNAfols_structure_format_start=strsplit(structure_of_RNAfold[3],split = "")[[1]][1:2000]

          RNAplot_file_name=paste0(">","final_structur","_1_to_2000")
          fileConn<-file("sec_stracuter.txt")
          writeLines(c(RNAplot_file_name,
                       paste(final_str_frags[1:2000],collapse = ""),
                       paste(RNAfols_structure_format_start,collapse = "")), fileConn)

          close(fileConn)


          if (host_system=="Windows") {
            RNAPlot_command_line=paste(RNAplot_location," -o  svg ", " sec_stracuter.txt")
            lapply(RNAPlot_command_line, system)
            file.remove("sec_stracuter.txt")
            file.copy("final_structur_1_to_2000_ss.svg",plots_folder_name)
            file.remove("final_structur_1_to_2000_ss.svg")

          }
          if (host_system=="Linux") {
            RNAPlot_command_line=paste(RNAplot_location, " sec_stracuter.txt ")
            lapply(RNAPlot_command_line, system)
            file.remove("sec_stracuter.txt")
            file.copy("final_structur_1_to_2000_ss.ps",plots_folder_name)
            file.remove("final_structur_1_to_2000_ss.ps")
          }



          # end
          write.fasta(final_str_frags[(length(final_str_frags)-2000+1):length(final_str_frags)],"final_str_frags","final_str_frags.fasta")
          structure_of_RNAfold=system(paste(RNAfold_path,"  ","final_str_frags.fasta" ),intern=T)
          file.remove("final_str_frags.fasta")
          file.remove("final_str_frags_ss.ps")
          RNAfols_structure_format_end=strsplit(structure_of_RNAfold[3],split = "")[[1]][1:2000]

          RNAplot_file_name=paste0(">","final_structur","_Last_2000")
          fileConn<-file("sec_stracuter.txt")
          writeLines(c(RNAplot_file_name,
                       paste(final_str_frags[(length(final_str_frags)-2000+1):length(final_str_frags)],collapse = "")
                       ,paste(RNAfols_structure_format_end,collapse = "")), fileConn)

          close(fileConn)
          if (host_system=="Windows") {
            RNAPlot_command_line=paste(RNAplot_location," -o  svg ", "  sec_stracuter.txt ")
            lapply(RNAPlot_command_line, system)
            file.remove("sec_stracuter.txt")
            file.copy("final_structur_Last_2000_ss.svg",plots_folder_name)
            file.remove("final_structur_Last_2000_ss.svg")


          }
          if (host_system=="Linux") {
            RNAPlot_command_line=paste(RNAplot_location, " sec_stracuter.txt ")
            lapply(RNAPlot_command_line, system)
            file.remove("sec_stracuter.txt")
            file.copy("final_structur_Last_2000_ss.ps",plots_folder_name)
            file.remove("final_structur_Last_2000_ss.ps")
          }

        }}
      ###########################################################








      ###################### for casVI_A type & casVI_D type ######################
      if (CRISPR_Type=="casVI_D"|CRISPR_Type=="casVI_A") {

        for (i in 1:(length(final_str_frags)-W+1)) {
          candidate=final_str_frags[i:(i+W-1)]
          candidate_scores=c(final_score1[i:(i+W-1)])
          GC_condidate=round(length(which(candidate=="g"| candidate=="c"))/length(candidate)*100,digits = 2)
          candidate_U_rich=(sum(final_score1[which(final_str_frags[i:(i+W-1)]==Nucleotide_rich)]))
          if (mismatch_del(candidate_scores,start,end,seed_start,seed_end)!=TRUE) {
            select_candidate=paste(revers_compliment(candidate),collapse = "")
            backbone_spacer=paste0(backbone,select_candidate,collapse = "")
            backbone_spacer_structure=pridict_structure(backbone_spacer,RNAfold_location = RNAfold_path)
            spaser_structure=pridict_structure(select_candidate,RNAfold_path)

            if (paste(backbone_spacer_structure$structure[1:length(backbone_structure)],collapse = "")==
                paste(backbone_structure,collapse = "")){
              final_free_anergy=backbone_spacer_structure$free_anergy_seq-backbone_free_enargy
              crRNA_free_energy[[select_candidate]]=paste(round(final_free_anergy,digits = 4),"Natural structure of DR")
              back_name=paste("Direct Repeat (Natural structure) ")
            }

            if(paste(backbone_spacer_structure$structure[1:length(backbone_structure)],collapse = "")!=
               paste(backbone_structure,collapse = "")){
              final_free_anergy=backbone_spacer_structure$free_anergy_seq
              crRNA_free_energy[[select_candidate]]=paste(round(final_free_anergy,digits = 4),"Unnatural structure of DR")
              back_name=paste("Direct Repeat (Unnatural structure) ")
            }

            prtospaser[[select_candidate]]=candidate

            crRNA_seq=c(crRNA_seq,select_candidate)
            vv[[select_candidate]]=round(revers_function(candidate_scores),digits = 4)
            start_pos[[select_candidate]]=i
            GC_content[[select_candidate]]=paste0(GC_condidate,"%")
            mismatch_pos[[select_candidate]]=which(revers_function(candidate_scores)<ConservationThreshold)
            mismatch_num[[select_candidate]]=length(which(candidate_scores<ConservationThreshold))
            rank[[select_candidate]]=length(which(candidate_scores<ConservationThreshold))+1
            efficiency[[select_candidate]]=mean(candidate_scores)
            Guide_RNA_ID[[select_candidate]]=paste0(ID_code,"_",ConservationMethod,"_",which(crRNA_seq==select_candidate))
            pos_spacer_in_sec_stracture=i
            if (length(final_str_frags)>2000){
              if (i<=1000) {
                RNAfols_structure_format=RNAfols_structure_format_start
                pos_spacer_in_sec_stracture=i
              }
              else if (i>1000 & i<length(final_str_frags)-1000+1) {

                write.fasta(final_str_frags[(i-1000):(i+1000-1)],"final_str_frags","final_str_frags.fasta")
                structure_of_RNAfold=system(paste(RNAfold_path,"  ","final_str_frags.fasta" ),intern=T)
                file.remove("final_str_frags.fasta")
                file.remove("final_str_frags_ss.ps")
                RNAfols_structure_format=strsplit(structure_of_RNAfold[3],split = "")[[1]][1:2000]

                RNAplot_file_name=paste0(">","final_structur_",(i-1000),"_to_",(i+1000-1))
                fileConn<-file("sec_stracuter.txt")
                writeLines(c(RNAplot_file_name,
                             paste(final_str_frags[(i-1000):(i+1000-1)],collapse = ""),
                             paste(RNAfols_structure_format,collapse = "")), fileConn)

                close(fileConn)
                if (host_system=="Windows") {
                  RNAPlot_command_line=paste(RNAplot_location," -o  svg ", " sec_stracuter.txt")
                  lapply(RNAPlot_command_line, system)
                  file.remove("sec_stracuter.txt")
                  file.copy(paste0("final_structur_",(i-1000),"_to_",(i+1000-1),"_ss.svg"),plots_folder_name)
                  file.remove(paste0("final_structur_",(i-1000),"_to_",(i+1000-1),"_ss.svg"))


                }
                if (host_system=="Linux") {
                  RNAPlot_command_line=paste(RNAplot_location, " sec_stracuter.txt ")
                  lapply(RNAPlot_command_line, system)
                  file.remove("sec_stracuter.txt")
                  file.copy(paste0("final_structur_",(i-1000),"_to_",(i+1000-1),"_ss.ps"),plots_folder_name)
                  file.remove(paste0("final_structur_",(i-1000),"_to_",(i+1000-1),"_ss.ps"))
                }

                pos_spacer_in_sec_stracture=1000
              }
              else if (i>=length(final_str_frags)-1000) {
                RNAfols_structure_format=RNAfols_structure_format_end
                pos_spacer_in_sec_stracture=2000-(length(final_str_frags)-i)

              }
            }

            if (i<=100) {
              local_U_rich[[select_candidate]]=
                round(((sum(final_score1[(which(final_str_frags[1:(i+W+100-1)]==Nucleotide_rich))]))-candidate_U_rich)/
                        (sum(final_score1[1:(i+W+100-1)])-sum(candidate_scores)),digits = 4)


              normalized_U_rich[[select_candidate]]= round(local_U_rich[[select_candidate]]/u_rich_in_whole,digits = 4)


              a_side=which(RNAfols_structure_format[pos_spacer_in_sec_stracture:(pos_spacer_in_sec_stracture+W-1)]==".")
              a_flank_side=which(RNAfols_structure_format[c((1:pos_spacer_in_sec_stracture-1),
                                                            (pos_spacer_in_sec_stracture+W):(pos_spacer_in_sec_stracture+W+100-1))]==".")
              accessible_side[[select_candidate]]=round(length(a_side)/W, digits=4)
              accessible_flank_side[[select_candidate]]=
                round(length(a_flank_side)/length(c((1:pos_spacer_in_sec_stracture-1),
                                                    (pos_spacer_in_sec_stracture+W):(pos_spacer_in_sec_stracture+W+100-1))), digits=4)

            }
            else if (i>100 & i<length(final_str_frags)-100-W+1) {
              local_U_rich[[select_candidate]]=
                round(((sum(final_score1[(which(final_str_frags[(i-100):(i+W+100-1)]==Nucleotide_rich))+i-100-1]))-candidate_U_rich)/
                        (sum(final_score1[(i-100):(i+W+100-1)])-sum(candidate_scores)),digits = 4)


              normalized_U_rich[[select_candidate]]=
                round(local_U_rich[[select_candidate]]/u_rich_in_whole,digits = 4)

              a_side=which(RNAfols_structure_format[pos_spacer_in_sec_stracture:(pos_spacer_in_sec_stracture+W-1)]==".")
              a_flank_side=which(RNAfols_structure_format[c((pos_spacer_in_sec_stracture-100):(pos_spacer_in_sec_stracture-1),
                                                            (pos_spacer_in_sec_stracture+W):(pos_spacer_in_sec_stracture+W+100-1))]==".")
              accessible_side[[select_candidate]]=round(length(a_side)/W, digits=4)
              accessible_flank_side[[select_candidate]]=
                round(length(a_flank_side)/length(c((pos_spacer_in_sec_stracture-100):(pos_spacer_in_sec_stracture-1),
                                                    (pos_spacer_in_sec_stracture+W):(pos_spacer_in_sec_stracture+W+100-1))), digits=4)
            }
            else if (i>=length(final_str_frags)-100-W) {
              local_U_rich[[select_candidate]]=
                round(((sum(final_score1[(which(final_str_frags[(i-100+1):length(final_str_frags)]==Nucleotide_rich))+i-100]))-candidate_U_rich)/
                        ((sum(final_score1[(i-100+1):length(final_str_frags)]))-sum(candidate_scores)),digits = 4)

              normalized_U_rich[[select_candidate]]= round(local_U_rich[[select_candidate]]/u_rich_in_whole,digits = 4)

              a_side=which(RNAfols_structure_format[pos_spacer_in_sec_stracture:(pos_spacer_in_sec_stracture+W-1)]==".")
              a_flank_side=which(RNAfols_structure_format[c((pos_spacer_in_sec_stracture-100):(pos_spacer_in_sec_stracture-1),
                                                            (pos_spacer_in_sec_stracture+W):2000)]==".")
              accessible_side[[select_candidate]]=round(length(a_side)/W, digits=4)
              accessible_flank_side[[select_candidate]]=
                round(length(a_flank_side)/length(c((pos_spacer_in_sec_stracture-100):(pos_spacer_in_sec_stracture-1),
                                                    (pos_spacer_in_sec_stracture+W):2000)), digits=4)

            }
            if (mismatch_num[[select_candidate]]>0){
              if (mismatch_num[[select_candidate]]==1) {
                sssnp=normal_score_matrix1[,original_pos_in_matrix[i+(W-mismatch_pos[[select_candidate]])]]
                sssnp[which(names(sssnp)==final_str_frags[i+(W-mismatch_pos[[select_candidate]])])]=0
                SNP=sssnp>0
                type_if_SNP[[select_candidate]]=
                  paste(strsplit(toString(names(SNP[SNP==TRUE]),),split = ", " )[[1]],sep = "",collapse = " & ")
              }
              else if (mismatch_num[[select_candidate]]==2) {
                type_if_SNP[[select_candidate]]=c()
                sssnp1=normal_score_matrix1[,original_pos_in_matrix[i+(W-mismatch_pos[[select_candidate]][1])]]
                sssnp1[which(names(sssnp1)==final_str_frags[i+(W-mismatch_pos[[select_candidate]][1])])]=0
                sssnp2=normal_score_matrix1[,original_pos_in_matrix[i+(W-mismatch_pos[[select_candidate]][2])]]
                sssnp2[which(names(sssnp2)==final_str_frags[i+(W-mismatch_pos[[select_candidate]][2])])]=0
                SNP1=sssnp1>0
                SNP2=sssnp2>0
                type_if_SNP[[select_candidate]]=
                  c(paste(strsplit(toString(names(SNP1[SNP1==TRUE]),),split = ", " )[[1]],sep = "",collapse = " & "),
                    paste(strsplit(toString(names(SNP2[SNP2==TRUE]),),split = ", " )[[1]],sep = "",collapse = " & "))
              }}

            if (mismatch_num[[select_candidate]]==0) {
              type_if_SNP[[select_candidate]]="No_SNP"
              mismatch_pos[[select_candidate]]="-"
            }
            ###File name for plot result
            title_for_crRNA_plot=paste0("Information plot for ","(ID=",Guide_RNA_ID[[select_candidate]],") candidate ")
            filename2=paste0("ID-",paste0(ID_code,"_",ConservationMethod,"_",which(crRNA_seq==select_candidate)),".png",collapse = "")

            if (i<=20) {
              setwd(plots_folder_name)
              png(filename2, width = 1300,height = 750)
              layout(matrix(c(2,2,3,3,4,4,2,2,3,3,4,4,1,1,1,1,1,1), 3, 6, byrow = TRUE))
              crRNA_plot(final_str_frags[1:(i+W+20-1)],i,start_pos[[select_candidate]],
                         mismatch_pos[[select_candidate]],type_if_SNP[[select_candidate]],
                         title_for_crRNA_plot,W,paste0(local_U_rich[[select_candidate]]," ",
                                                       local_rich_type," , ",normalized_U_rich[[select_candidate]]," ",Normal_rich_type)
                         ,paste(accessible_flank_side[[select_candidate]]),
                         paste(efficiency[[select_candidate]],"efficiency"),
                         paste(accessible_side[[select_candidate]]))

              my_Dot_Plot(select_candidate,
                          select_candidate,
                          points_col = "cyan3",
                          xlab = "Specer Seq",
                          ylab = "Specer Seq")

              sec_structure_plot(select_candidate,
                                 spaser_structure$structure,
                                 first_backbone = "N",
                                 main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                 free_energy = spaser_structure$free_anergy_seq)

              sec_structure_plot(backbone_spacer,
                                 backbone_spacer_structure$structure,
                                 first_backbone,
                                 back_name,
                                 length(backbone_structure),
                                 paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

              dev.off()
              setwd(first_dir)

              link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                mutate( path = paste0('file:', plots_folder_name, "/", sort(file ))
                        , path = stringr::str_replace_all( path, ' ', '%20')
                        , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
              link_plots[[select_candidate]]=link$link
            }
            if (i>20 & length(final_str_frags)-20-W+1) {
              setwd(plots_folder_name)
              png(filename2, width = 1300,height = 750)
              layout(matrix(c(2,2,3,3,4,4,2,2,3,3,4,4,1,1,1,1,1,1), 3, 6, byrow = TRUE))
              crRNA_plot(final_str_frags[(i-20):(i+W+20-1)],21,
                         start_pos[[select_candidate]],
                         mismatch_pos[[select_candidate]],
                         type_if_SNP[[select_candidate]],title_for_crRNA_plot,W,
                         paste0(local_U_rich[[select_candidate]]," ",local_rich_type," , ",
                                normalized_U_rich[[select_candidate]]," ",Normal_rich_type),
                         paste(accessible_flank_side[[select_candidate]]),
                         paste(efficiency[[select_candidate]],"efficiency"),
                         paste(accessible_side[[select_candidate]]))

              my_Dot_Plot(select_candidate,
                          select_candidate,
                          points_col = "cyan3",
                          xlab = "Specer Seq",
                          ylab = "Specer Seq")

              sec_structure_plot(select_candidate,
                                 spaser_structure$structure,
                                 first_backbone = "N",
                                 main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                 free_energy = spaser_structure$free_anergy_seq)

              sec_structure_plot(backbone_spacer,
                                 backbone_spacer_structure$structure,
                                 first_backbone,
                                 back_name,
                                 length(backbone_structure),
                                 paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

              dev.off()
              setwd(first_dir)
              link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                mutate( path = paste0('file:', plots_folder_name, "/", sort(file ))
                        , path = stringr::str_replace_all( path, ' ', '%20')
                        , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
              link_plots[[select_candidate]]=link$link
            }
            if (i>=length(final_str_frags)-20-W) {
              setwd(plots_folder_name)
              png(filename2, width = 1300,height = 750)
              layout(matrix(c(2,2,3,3,4,4,2,2,3,3,4,4,1,1,1,1,1,1), 3, 6, byrow = TRUE))
              crRNA_plot(final_str_frags[(i-20+1):length(final_str_frags)],21,
                         start_pos[[select_candidate]],
                         mismatch_pos[[select_candidate]],
                         type_if_SNP[[select_candidate]],
                         title_for_crRNA_plot,W,
                         paste0(local_U_rich[[select_candidate]]," ",
                                local_rich_type," , ",
                                normalized_U_rich[[select_candidate]]," ",
                                Normal_rich_type)
                         ,paste(accessible_flank_side[[select_candidate]]),
                         paste(efficiency[[select_candidate]],"efficiency"),
                         paste(accessible_side[[select_candidate]]))


              my_Dot_Plot(select_candidate,
                          select_candidate,
                          points_col = "cyan3",
                          xlab = "Specer Seq",
                          ylab = "Specer Seq")

              sec_structure_plot(select_candidate,
                                 spaser_structure$structure,
                                 first_backbone = "N",
                                 main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                 free_energy = spaser_structure$free_anergy_seq)

              sec_structure_plot(backbone_spacer,
                                 backbone_spacer_structure$structure,
                                 first_backbone,
                                 back_name,
                                 length(backbone_structure),
                                 paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

              dev.off()
              setwd(first_dir)
              link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                mutate( path = paste0('file:', plots_folder_name, "/", sort(file ))
                        , path = stringr::str_replace_all( path, ' ', '%20')
                        , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
              link_plots[[select_candidate]]=link$link
            }
            score_Gaide[[select_candidate]]=sum(c((75-(mismatch_num[[select_candidate]]*5))
                                                  ,(10*2*local_U_rich[[select_candidate]])
                                                  ,(2.5*accessible_side[[select_candidate]])
                                                  ,(2.5*accessible_flank_side[[select_candidate]])
                                                  ,(10*(1-(final_free_anergy/max_free_energy)))))
          }
        }
        blast_query=prtospaser
      }
      #############################################################################









      ############################# for casVI_B type ##########################
      if (CRISPR_Type=="casVI_B") {

        prtospaseraa=list()
        prtospasergg=list()
        prtospasertt=list()
        prtospaserag=list()
        prtospaserga=list()
        prtospaserat=list()
        prtospaserta=list()
        prtospasergt=list()
        prtospasertg=list()

        for (i in 3:(length(final_str_frags)-W+1)) {
          candidate=final_str_frags[i:(i+W-1)]
          candidate_scores=c(final_score1[i:(i+W-1)])
          GC_condidate=round(length(which(candidate=="g"| candidate=="c"))/length(candidate)*100,digits = 2)
          candidate_U_rich=(sum(final_score1[which(final_str_frags[i:(i+W-1)]==Nucleotide_rich)]))
          if (final_str_frags[i-1]  != "c" & final_str_frags[i-2] != "c"){
            if (mismatch_del(candidate_scores,start,end,seed_start,seed_end)!=TRUE) {
              select_candidate=paste(revers_compliment(candidate),collapse = "")

              backbone_spacer=paste0(select_candidate,backbone,collapse = "")
              backbone_spacer_structure=pridict_structure(backbone_spacer,RNAfold_location = RNAfold_path)
              spaser_structure=pridict_structure(select_candidate,RNAfold_path)

              if (paste(backbone_spacer_structure$structure[(W+1):length(backbone_spacer_structure$structure)],collapse = "")==
                  paste(backbone_structure,collapse = "")){
                final_free_anergy=backbone_spacer_structure$free_anergy_seq-backbone_free_enargy
                crRNA_free_energy[[select_candidate]]=paste(round(final_free_anergy,digits = 4),"Natural structure of DR")
                back_name=paste("Direct Repeat (Natural structure) ")
              }

              if(paste(backbone_spacer_structure$structure[(W+1):length(backbone_spacer_structure$structure)],collapse = "")!=
                 paste(backbone_structure,collapse = "")){
                final_free_anergy=backbone_spacer_structure$free_anergy_seq
                crRNA_free_energy[[select_candidate]]=paste(round(final_free_anergy,digits = 4),"Unnatural structure of DR")
                back_name=paste("Direct Repeat (Unnatural structure) ")
              }

              prtospaser[[select_candidate]]=candidate
              prtospaseraa[[select_candidate]]=c("a","a",candidate)
              prtospasergg[[select_candidate]]=c("g","g",candidate)
              prtospasertt[[select_candidate]]=c("t","t",candidate)
              prtospaserag[[select_candidate]]=c("a","g",candidate)
              prtospaserga[[select_candidate]]=c("g","a",candidate)
              prtospaserat[[select_candidate]]=c("a","t",candidate)
              prtospaserta[[select_candidate]]=c("t","a",candidate)
              prtospasergt[[select_candidate]]=c("g","t",candidate)
              prtospasertg[[select_candidate]]=c("t","g",candidate)
              crRNA_seq=c(crRNA_seq,select_candidate)
              vv[[select_candidate]]=round(revers_function(candidate_scores),digits = 4)
              start_pos[[select_candidate]]=i
              GC_content[[select_candidate]]=paste0(GC_condidate,"%")
              mismatch_pos[[select_candidate]]=which(revers_function(candidate_scores)<ConservationThreshold)
              mismatch_num[[select_candidate]]=length(which(candidate_scores<ConservationThreshold))
              rank[[select_candidate]]=length(which(candidate_scores<ConservationThreshold))+1
              efficiency[[select_candidate]]=mean(candidate_scores)
              Guide_RNA_ID[[select_candidate]]=paste0(ID_code,"_",ConservationMethod,"_",which(crRNA_seq==select_candidate))
              pos_spacer_in_sec_stracture=i
              if (length(final_str_frags)>2000){
                if (i<=1000) {
                  RNAfols_structure_format=RNAfols_structure_format_start
                  pos_spacer_in_sec_stracture=i
                }
                else if (i>1000 & i<length(final_str_frags)-1000+1) {
                  write.fasta(final_str_frags[(i-1000):(i+1000-1)],"final_str_frags","final_str_frags.fasta")
                  structure_of_RNAfold=system(paste(RNAfold_path,"  ","final_str_frags.fasta" ),intern=T)
                  file.remove("final_str_frags.fasta")
                  file.remove("final_str_frags_ss.ps")
                  RNAfols_structure_format=strsplit(structure_of_RNAfold[3],split = "")[[1]][1:2000]

                  RNAplot_file_name=paste0(">","final_structur_",(i-1000),"_to_",(i+1000-1))
                  fileConn<-file("sec_stracuter.txt")
                  writeLines(c(RNAplot_file_name,paste(final_str_frags[(i-1000):(i+1000-1)],collapse = ""),
                               paste(RNAfols_structure_format_start,collapse = "")), fileConn)

                  close(fileConn)
                  if (host_system=="Windows") {
                    RNAPlot_command_line=paste(RNAplot_location," -o  svg ", " sec_stracuter.txt")
                    lapply(RNAPlot_command_line, system)
                    file.remove("sec_stracuter.txt")
                    file.copy(paste0("final_structur_",(i-1000),"_to_",(i+1000-1),"_ss.svg"),plots_folder_name)
                    file.remove(paste0("final_structur_",(i-1000),"_to_",(i+1000-1),"_ss.svg"))


                  }
                  if (host_system=="Linux") {
                    RNAPlot_command_line=paste(RNAplot_location, " sec_stracuter.txt ")
                    lapply(RNAPlot_command_line, system)
                    file.remove("sec_stracuter.txt")
                    file.copy(paste0("final_structur_",(i-1000),"_to_",(i+1000-1),"_ss.ps"),plots_folder_name)
                    file.remove(paste0("final_structur_",(i-1000),"_to_",(i+1000-1),"_ss.ps"))
                  }

                  pos_spacer_in_sec_stracture=1000
                }
                else if (i>=length(final_str_frags)-1000) {
                  RNAfols_structure_format=RNAfols_structure_format_end
                  pos_spacer_in_sec_stracture=2000-(length(final_str_frags)-i)
                }
              }

              if (i<=100) {
                local_U_rich[[select_candidate]]=
                  round(((sum(final_score1[(which(final_str_frags[1:(i+W+100-1)]==Nucleotide_rich))]))-candidate_U_rich)/
                          (sum(final_score1[1:(i+W+100-1)])-sum(candidate_scores)),digits = 4)

                normalized_U_rich[[select_candidate]]= round(local_U_rich[[select_candidate]]/u_rich_in_whole,digits = 4)

                a_side=which(RNAfols_structure_format[i:(i+W-1)]==".")
                a_flank_side=which(RNAfols_structure_format[c((1:i-1),(i+W):(i+W+100-1))]==".")
                accessible_side[[select_candidate]]=round(length(a_side)/W, digits=4)

                accessible_flank_side[[select_candidate]]=round(length(a_flank_side)/
                                                                  length(c((1:i-1),(i+W):(i+W+100-1))), digits=4)
              }
              else if (i>100 & i<length(final_str_frags)-100-W+1) {
                local_U_rich[[select_candidate]]=
                  round(((sum(final_score1[(which(final_str_frags[(i-100):(i+W+100-1)]==Nucleotide_rich)+i-100)]))-candidate_U_rich)/
                          (sum(final_score1[(i-100):(i+W+100-1)])-sum(candidate_scores)),digits = 4)

                normalized_U_rich[[select_candidate]]=round(local_U_rich[[select_candidate]]/u_rich_in_whole,digits = 4)

                a_side=which(RNAfols_structure_format[i:(i+W-1)]==".")
                a_flank_side=which(RNAfols_structure_format[c((i-100):(i-1),(i+W):(i+W+100-1))]==".")
                accessible_side[[select_candidate]]=round(length(a_side)/W, digits=4)

                accessible_flank_side[[select_candidate]]=
                  round(length(a_flank_side)/length(c((i-100):(i-1),(i+W):(i+W+100-1))), digits=4)

              }
              else if (i>=length(final_str_frags)-100-W) {
                local_U_rich[[select_candidate]]=
                  round(((sum(final_score1[(which(final_str_frags[(i-100+1):length(final_str_frags)]==Nucleotide_rich)+i-100)]))-candidate_U_rich)/
                          ((sum(final_score1[(i-100+1):length(final_str_frags)]))-sum(candidate_scores)),digits = 4)

                normalized_U_rich[[select_candidate]]= round(local_U_rich[[select_candidate]]/u_rich_in_whole,digits = 4)

                a_side=which(RNAfols_structure_format[i:(i+W-1)]==".")
                a_flank_side=which(RNAfols_structure_format[c((i-100):(i-1),(i+W-1):length(final_str_frags))]==".")
                accessible_side[[select_candidate]]=round(length(a_side)/W, digits=4)

                accessible_flank_side[[select_candidate]]=
                  round(length(a_flank_side)/length(c((i-100):(i-1),(i+W-1):length(final_str_frags))), digits=4)
              }
              if (mismatch_num[[select_candidate]]>0){
                if (mismatch_num[[select_candidate]]==1) {
                  sssnp=normal_score_matrix1[,original_pos_in_matrix[i+(W-mismatch_pos[[select_candidate]])]]
                  sssnp[which(names(sssnp)==final_str_frags[i+(W-mismatch_pos[[select_candidate]])])]=0
                  SNP=sssnp>0

                  type_if_SNP[[select_candidate]]=
                    paste(strsplit(toString(names(SNP[SNP==TRUE]),),split = ", ")[[1]],collapse = " & ",sep="")
                }
                else if (mismatch_num[[select_candidate]]==2) {
                  type_if_SNP[[select_candidate]]=c()
                  sssnp1=normal_score_matrix1[,original_pos_in_matrix[i+(W-mismatch_pos[[select_candidate]][1])]]
                  sssnp1[which(names(sssnp1)==final_str_frags[i+(W-mismatch_pos[[select_candidate]][1])])]=0
                  sssnp2=normal_score_matrix1[,original_pos_in_matrix[i+(W-mismatch_pos[[select_candidate]][2])]]
                  sssnp2[which(names(sssnp2)==final_str_frags[i+(W-mismatch_pos[[select_candidate]][2])])]=0
                  SNP1=sssnp1>0
                  SNP2=sssnp2>0
                  type_if_SNP[[select_candidate]]=
                    c(paste(strsplit(toString(names(SNP1[SNP1==TRUE]),),split = ", "),collapse = " & ",sep=""),
                      paste(strsplit(toString(names(SNP2[SNP2==TRUE]),),split = ", "),collapse = " & ",sep=""))
                }}

              if (mismatch_num[[select_candidate]]==0) {
                type_if_SNP[[select_candidate]]="No_SNP"
                mismatch_pos[[select_candidate]]="-"
              }
              ###File name for plot result
              title_for_crRNA_plot=paste0("Information plot for ","(ID=",Guide_RNA_ID[[select_candidate]],") candidate ")
              filename2=paste0("ID-",paste0(ID_code,"_",ConservationMethod,"_",which(crRNA_seq==select_candidate)),".png",collapse = "")

              if (i<=20) {
                setwd(plots_folder_name)
                png(filename2, width = 1300,height = 750)
                layout(matrix(c(2,2,3,3,4,4,2,2,3,3,4,4,1,1,1,1,1,1), 3, 6, byrow = TRUE))

                crRNA_plot(final_str_frags[1:(i+W+20-1)],i,
                           start_pos[[select_candidate]],
                           mismatch_pos[[select_candidate]],
                           type_if_SNP[[select_candidate]],
                           title_for_crRNA_plot,W,paste0(local_U_rich[[select_candidate]]," ",
                                                         local_rich_type," , ",normalized_U_rich[[select_candidate]]," ",
                                                         Normal_rich_type),paste(accessible_flank_side[[select_candidate]]),
                           paste(efficiency[[select_candidate]],"efficiency"),
                           paste(accessible_side[[select_candidate]]))

                my_Dot_Plot(select_candidate,
                            select_candidate,
                            points_col = "cyan3",
                            xlab = "Specer Seq",
                            ylab = "Specer Seq")

                sec_structure_plot(select_candidate,
                                   spaser_structure$structure,
                                   first_backbone = "N",
                                   main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                   free_energy = spaser_structure$free_anergy_seq)

                sec_structure_plot(backbone_spacer,
                                   backbone_spacer_structure$structure,
                                   first_backbone,
                                   back_name,
                                   length(backbone_structure),
                                   paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

                dev.off()
                setwd(first_dir)
                link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                  mutate( path = paste0('file:', plots_folder_name,"/" , sort(file ))
                          , path = stringr::str_replace_all( path, ' ', '%20')
                          , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
                link_plots[[select_candidate]]=link$link
              }
              if (i>20 & length(final_str_frags)-20-W+1) {
                setwd(plots_folder_name)
                png(filename2, width = 1300,height = 750)
                layout(matrix(c(2,2,3,3,4,4,2,2,3,3,4,4,1,1,1,1,1,1), 3, 6, byrow = TRUE))
                crRNA_plot(final_str_frags[(i-20):(i+W+20-1)],21,
                           start_pos[[select_candidate]],
                           mismatch_pos[[select_candidate]],
                           type_if_SNP[[select_candidate]],
                           title_for_crRNA_plot,W,paste0(local_U_rich[[select_candidate]]," ",
                                                         local_rich_type," , ",normalized_U_rich[[select_candidate]]," ",
                                                         Normal_rich_type),paste(accessible_flank_side[[select_candidate]]),
                           paste(efficiency[[select_candidate]],"efficiency"),
                           paste(accessible_side[[select_candidate]]))

                my_Dot_Plot(select_candidate,
                            select_candidate,
                            points_col = "cyan3",
                            xlab = "Specer Seq",
                            ylab = "Specer Seq")

                sec_structure_plot(select_candidate,
                                   spaser_structure$structure,
                                   first_backbone = "N",
                                   main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                   free_energy = spaser_structure$free_anergy_seq)

                sec_structure_plot(backbone_spacer,
                                   backbone_spacer_structure$structure,
                                   first_backbone,
                                   back_name,
                                   length(backbone_structure),
                                   paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))


                dev.off()
                setwd(first_dir)
                link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                  mutate( path = paste0('file:', plots_folder_name,"/" , sort(file ))
                          , path = stringr::str_replace_all( path, ' ', '%20')
                          , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
                link_plots[[select_candidate]]=link$link
              }
              if (i>=length(final_str_frags)-20-W) {
                setwd(plots_folder_name)
                png(filename2, width = 1300,height = 750)
                layout(matrix(c(2,2,3,3,4,4,2,2,3,3,4,4,1,1,1,1,1,1), 3, 6, byrow = TRUE))
                crRNA_plot(final_str_frags[(i-20+1):length(final_str_frags)],21,
                           start_pos[[select_candidate]],
                           mismatch_pos[[select_candidate]],
                           type_if_SNP[[select_candidate]],
                           title_for_crRNA_plot,W,paste0(local_U_rich[[select_candidate]]," ",
                                                         local_rich_type," , ",normalized_U_rich[[select_candidate]]," ",
                                                         Normal_rich_type),paste(accessible_flank_side[[select_candidate]]),
                           paste(efficiency[[select_candidate]],"efficiency"),
                           paste(accessible_side[[select_candidate]]))


                my_Dot_Plot(select_candidate,
                            select_candidate,
                            points_col = "cyan3",
                            xlab = "Specer Seq",
                            ylab = "Specer Seq")

                sec_structure_plot(select_candidate,
                                   spaser_structure$structure,
                                   first_backbone = "N",
                                   main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                   free_energy = spaser_structure$free_anergy_seq)

                sec_structure_plot(backbone_spacer,
                                   backbone_spacer_structure$structure,
                                   first_backbone,
                                   back_name,
                                   length(backbone_structure),
                                   paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

                dev.off()
                setwd(first_dir)
                link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                  mutate( path = paste0('file:', plots_folder_name,"/" , sort(file ))
                          , path = stringr::str_replace_all( path, ' ', '%20')
                          , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
                link_plots[[select_candidate]]=link$link
              }

              score_Gaide[[select_candidate]]=sum(c((75-(mismatch_num[[select_candidate]]*5))
                                                    ,(10*2*local_U_rich[[select_candidate]])
                                                    ,(2.5*accessible_side[[select_candidate]])
                                                    ,(2.5*accessible_flank_side[[select_candidate]])
                                                    ,(10*(1-(final_free_anergy/max_free_energy)))))
            }}


        }
        blast_query=c(prtospaseraa,
                      prtospasergg,
                      prtospasertt,
                      prtospaserag,
                      prtospaserga,
                      prtospaserat,
                      prtospaserta,
                      prtospasergt,
                      prtospasertg)
      }
      #########################################################################









      ############################# for casV_B type ###########################
      if (CRISPR_Type=="casV_B") {
        final_str_frags_plus_vb=c("l","l","l",final_str_frags,"l","l","l")
        prtospasertta_pos=list()
        prtospaserttg_pos=list()
        prtospaserttt_pos=list()
        prtospaserttc_pos=list()
        prtospasertta_neg=list()
        prtospaserttg_neg=list()
        prtospaserttt_neg=list()
        prtospaserttc_neg=list()
        for (i in 4:(length(final_str_frags_plus_vb)-W+1-6)) {
          candidate=final_str_frags_plus_vb[i:(i+W-1)]
          candidate_scores1=c(final_score1[(i-3):(i+W-1-3)])
          GC_condidate=round(length(which(candidate=="g"| candidate=="c"))/length(candidate)*100,digits = 2)

          if ((final_str_frags_plus_vb[i-2]    == "t" & final_str_frags_plus_vb[i-3]   == "t")|
              (final_str_frags_plus_vb[i+W+1]  == "a" & final_str_frags_plus_vb[i+W+2] == "a")){

            if (mismatch_del(candidate_scores1,start,end,seed_start,seed_end)!=TRUE) {

              if (final_str_frags_plus_vb[i+W+2]  == "a"  &
                  final_str_frags_plus_vb[i+W+1]  == "a") {

                select_candidate=paste(revers_compliment(candidate),collapse = "")
                candidate_scores=revers_function(candidate_scores1)

                backbone_spacer=paste0(backbone,select_candidate,collapse = "")
                backbone_spacer_structure=pridict_structure(backbone_spacer,RNAfold_location = RNAfold_path)
                spaser_structure=pridict_structure(select_candidate,RNAfold_path)

                if (paste(backbone_spacer_structure$structure[1:length(backbone_structure)],collapse = "")==
                    paste(backbone_structure,collapse = "")){
                  final_free_anergy=backbone_spacer_structure$free_anergy_seq-backbone_free_enargy
                  crRNA_free_energy[[select_candidate]]=paste(round(final_free_anergy,digits = 4),"Natural structure of DR")
                  back_name=paste("Direct Repeat (Natural structure) ")
                }

                if(paste(backbone_spacer_structure$structure[1:length(backbone_structure)],collapse = "")!=
                   paste(backbone_structure,collapse = "")){
                  final_free_anergy=backbone_spacer_structure$free_anergy_seq
                  crRNA_free_energy[[select_candidate]]=paste(round(final_free_anergy,digits = 4),"Unnatural structure of DR")
                  back_name=paste("Direct Repeat (Unnatural structure) ")
                }



                prtospaser[[select_candidate]]=revers_compliment(candidate)
                GC_content[[select_candidate]]=paste0(GC_condidate,"%")
                vv[[select_candidate]]=round((candidate_scores),digits = 4)
                start_pos[[select_candidate]]=(length(final_str_frags_plus_vb)-6)-((i-3)+W-1)+1
                mismatch_pos[[select_candidate]]=which((candidate_scores)<ConservationThreshold)
                mismatch_num[[select_candidate]]=length(which(candidate_scores<ConservationThreshold))
                rank[[select_candidate]]=length(which(candidate_scores<ConservationThreshold))+1
                efficiency[[select_candidate]]=mean(candidate_scores)
                strand[[select_candidate]]="-"
                crRNA_seq=c(crRNA_seq,select_candidate)
                Guide_RNA_ID[[select_candidate]]=paste0(ID_code,"_",ConservationMethod,"_",which(crRNA_seq==select_candidate))
                prtospasertta_neg[[select_candidate]]=c("t","t","a",revers_compliment(candidate))
                prtospaserttg_neg[[select_candidate]]=c("t","t","g",revers_compliment(candidate))
                prtospaserttt_neg[[select_candidate]]=c("t","t","t",revers_compliment(candidate))
                prtospaserttc_neg[[select_candidate]]=c("t","t","c",revers_compliment(candidate))
                if (mismatch_num[[select_candidate]]>0){
                  if (mismatch_num[[select_candidate]]==1) {
                    sssnp=normal_score_matrix1[,original_pos_in_matrix[(i-3)+(W-mismatch_pos[[select_candidate]])]]
                    sssnp[which(names(sssnp)==final_str_frags_plus_vb[(i)+(W-mismatch_pos[[select_candidate]])])]=0
                    SNP=sssnp>0

                    type_if_SNP[[select_candidate]]=
                      paste(revers_compliment(strsplit(toString(names(SNP[SNP==TRUE]),),split = ", " )[[1]]),sep = "",collapse = " & ")

                  }
                  else if (mismatch_num[[select_candidate]]==2) {
                    type_if_SNP[[select_candidate]]=c()
                    sssnp1=normal_score_matrix1[,original_pos_in_matrix[(i-3)+(W-mismatch_pos[[select_candidate]][1])]]
                    sssnp1[which(names(sssnp1)==final_str_frags_plus_vb[(i)+(W-mismatch_pos[[select_candidate]][1])])]=0
                    sssnp2=normal_score_matrix1[,original_pos_in_matrix[(i-3)+(W-mismatch_pos[[select_candidate]][2])]]
                    sssnp2[which(names(sssnp2)==final_str_frags_plus_vb[(i)+(W-mismatch_pos[[select_candidate]][2])])]=0

                    SNP1=sssnp1>0
                    SNP2=sssnp2>0

                    type_if_SNP[[select_candidate]]=c(
                      paste(revers_compliment(strsplit(toString(names(SNP1[SNP1==TRUE]),),split = ", " )[[1]]),sep = "",collapse = " & "),
                      paste(revers_compliment(strsplit(toString(names(SNP2[SNP2==TRUE]),),split = ", " )[[1]]),sep = "",collapse = " & "))



                  }}
                if (mismatch_num[[select_candidate]]==0) {
                  type_if_SNP[[select_candidate]]="No_SNP"
                  mismatch_pos[[select_candidate]]="-"
                }
                ###File name for plot result
                title_for_crRNA_plot_DNA=paste0("Information plot for ",Guide_RNA_ID[[select_candidate]], " crRNA candidate")
                title_for_crRNA_Dotplot=paste0("Dotplot for ", Guide_RNA_ID[[select_candidate]]," crRNA candidate")
                filename2=paste0("ID-",paste0(ID_code,"_",ConservationMethod,"_",which(crRNA_seq==select_candidate)),".png",collapse = "")


                if (i<=20) {
                  setwd(plots_folder_name)
                  png(filename2, width = 1500,height = 866)
                  layout(matrix(c(2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  1,1,1,1,1,1,
                                  1,1,1,1,1,1,
                                  1,1,1,1,1,1), 8,6, byrow = TRUE))

                  crRNA_plot_DNA(final_str_frags_plus_vb[4:(i+W+20-1)],i,
                                 start_pos[[select_candidate]],
                                 mismatch_pos[[select_candidate]],
                                 type_if_SNP[[select_candidate]],
                                 title_for_crRNA_plot_DNA,W,
                                 paste(efficiency[[select_candidate]],"Efficiency"),
                                 strand[[select_candidate]],3)

                  my_Dot_Plot(select_candidate,
                              select_candidate,
                              points_col = "cyan3",
                              xlab = "Specer Seq",
                              ylab = "Specer Seq")

                  sec_structure_plot(select_candidate,
                                     spaser_structure$structure,
                                     first_backbone = "N",
                                     main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                     free_energy = spaser_structure$free_anergy_seq)

                  sec_structure_plot(backbone_spacer,
                                     backbone_spacer_structure$structure,
                                     first_backbone,
                                     back_name,
                                     length(backbone_structure),
                                     paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

                  dev.off()
                  setwd(first_dir)



                  link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                    mutate( path = paste0('file:', plots_folder_name,"/" , sort(file ))
                            , path = stringr::str_replace_all( path, ' ', '%20')
                            , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )



                  link_plots[[select_candidate]]=link$link
                }
                if (i>20 & i<length(final_str_frags_plus_vb)-20-W+1-3) {
                  setwd(plots_folder_name)
                  png(filename2, width = 1500,height = 866)
                  layout(matrix(c(2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  1,1,1,1,1,1,
                                  1,1,1,1,1,1,
                                  1,1,1,1,1,1), 8,6, byrow = TRUE))


                  crRNA_plot_DNA(final_str_frags_plus_vb[(i-20):(i+W+20-1)],21,
                                 start_pos[[select_candidate]],
                                 mismatch_pos[[select_candidate]],
                                 type_if_SNP[[select_candidate]],
                                 title_for_crRNA_plot_DNA,W,
                                 paste(efficiency[[select_candidate]],"Efficiency"),
                                 strand[[select_candidate]],3)

                  my_Dot_Plot(select_candidate,
                              select_candidate,
                              points_col = "cyan3",
                              xlab = "Specer Seq",
                              ylab = "Specer Seq")

                  sec_structure_plot(select_candidate,
                                     spaser_structure$structure,
                                     first_backbone = "N",
                                     main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                     free_energy = spaser_structure$free_anergy_seq)

                  sec_structure_plot(backbone_spacer,
                                     backbone_spacer_structure$structure,
                                     first_backbone,
                                     back_name,
                                     length(backbone_structure),
                                     paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

                  dev.off()
                  setwd(first_dir)


                  link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                    mutate( path = paste0('file:', plots_folder_name,"/" , sort(file ))
                            , path = stringr::str_replace_all( path, ' ', '%20')
                            , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )





                  link_plots[[select_candidate]]=link$link
                }
                if (i>=length(final_str_frags_plus_vb)-20-W-3) {
                  setwd(plots_folder_name)
                  png(filename2, width = 1500,height = 866)
                  layout(matrix(c(2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  1,1,1,1,1,1,
                                  1,1,1,1,1,1,
                                  1,1,1,1,1,1), 8,6, byrow = TRUE))

                  crRNA_plot_DNA(final_str_frags_plus_vb[(i-20):(length(final_str_frags_plus_vb)-3)],21,
                                 start_pos[[select_candidate]],mismatch_pos[[select_candidate]],
                                 type_if_SNP[[select_candidate]],title_for_crRNA_plot_DNA,W,
                                 paste(efficiency[[select_candidate]],"Efficiency"),strand[[select_candidate]],3)

                  my_Dot_Plot(select_candidate,
                              select_candidate,
                              points_col = "cyan3",
                              xlab = "Specer Seq",
                              ylab = "Specer Seq")

                  sec_structure_plot(select_candidate,
                                     spaser_structure$structure,
                                     first_backbone = "N",
                                     main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                     free_energy = spaser_structure$free_anergy_seq)

                  sec_structure_plot(backbone_spacer,
                                     backbone_spacer_structure$structure,
                                     first_backbone,
                                     back_name,
                                     length(backbone_structure),
                                     paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

                  dev.off()
                  setwd(first_dir)
                  link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                    mutate( path = paste0('file:', plots_folder_name,"/" , sort(file ))
                            , path = stringr::str_replace_all( path, ' ', '%20')
                            , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
                  link_plots[[select_candidate]]=link$link
                }

                score_Gaide[[select_candidate]]=sum(c((75-(mismatch_num[[select_candidate]]*5))
                                                      ,(10*(1-(final_free_anergy/max_free_energy)))))

              }






              if (final_str_frags_plus_vb[i-2]  == "t" &
                  final_str_frags_plus_vb[i-3]  == "t") {
                select_candidate=paste(candidate,collapse = "")

                backbone_spacer=paste0(backbone,select_candidate,collapse = "")
                backbone_spacer_structure=pridict_structure(backbone_spacer,RNAfold_location = RNAfold_path)
                spaser_structure=pridict_structure(select_candidate,RNAfold_path)

                if (paste(backbone_spacer_structure$structure[1:length(backbone_structure)],collapse = "")==
                    paste(backbone_structure,collapse = "")){
                  final_free_anergy=backbone_spacer_structure$free_anergy_seq-backbone_free_enargy
                  crRNA_free_energy[[select_candidate]]=paste(round(final_free_anergy,digits = 4),"Natural structure of DR")
                  back_name=paste("Direct Repeat (Natural structure) ")
                }

                if(paste(backbone_spacer_structure$structure[1:length(backbone_structure)],collapse = "")!=
                   paste(backbone_structure,collapse = "")){
                  final_free_anergy=backbone_spacer_structure$free_anergy_seq
                  crRNA_free_energy[[select_candidate]]=paste(round(final_free_anergy,digits = 4),"Unnatural structure of DR")
                  back_name=paste("Direct Repeat (Unnatural structure) ")
                }



                candidate_scores=candidate_scores1
                prtospaser[[select_candidate]]=candidate
                GC_content[[select_candidate]]=paste0(GC_condidate,"%")
                vv[[select_candidate]]=round((candidate_scores),digits = 4)
                start_pos[[select_candidate]]=i-3
                mismatch_pos[[select_candidate]]=which((candidate_scores)<ConservationThreshold)
                mismatch_num[[select_candidate]]=length(which(candidate_scores<ConservationThreshold))
                rank[[select_candidate]]=length(which(candidate_scores<ConservationThreshold))+1
                efficiency[[select_candidate]]=mean(candidate_scores)
                strand[[select_candidate]]="+"
                crRNA_seq=c(crRNA_seq,select_candidate)
                Guide_RNA_ID[[select_candidate]]=paste0(ID_code,"_",ConservationMethod,"_",which(crRNA_seq==select_candidate))
                prtospasertta_pos[[select_candidate]]=c("t","t","a",candidate)
                prtospaserttg_pos[[select_candidate]]=c("t","t","g",candidate)
                prtospaserttt_pos[[select_candidate]]=c("t","t","t",candidate)
                prtospaserttc_pos[[select_candidate]]=c("t","t","c",candidate)


                if (mismatch_num[[select_candidate]]>0){
                  if (mismatch_num[[select_candidate]]==1) {
                    sssnp=normal_score_matrix1[,original_pos_in_matrix[(i-3)+(mismatch_pos[[select_candidate]]-1)]]
                    sssnp[which(names(sssnp)==final_str_frags_plus_vb[(i)+(mismatch_pos[[select_candidate]]-1)])]=0
                    SNP=sssnp>0

                    type_if_SNP[[select_candidate]]=
                      paste(strsplit(toString(names(SNP[SNP==TRUE]),),split = ", " )[[1]],sep = "",collapse = " & ")
                  }
                  else if (mismatch_num[[select_candidate]]==2) {
                    type_if_SNP[[select_candidate]]=c()

                    type_if_SNP[[select_candidate]]=c()
                    sssnp1=normal_score_matrix1[,original_pos_in_matrix[(i-3)+(mismatch_pos[[select_candidate]][1]-1)]]
                    sssnp1[which(names(sssnp1)==final_str_frags_plus_vb[(i)+(mismatch_pos[[select_candidate]][1]-1)])]=0
                    sssnp2=normal_score_matrix1[,original_pos_in_matrix[(i-3)+(mismatch_pos[[select_candidate]][2]-1)]]
                    sssnp2[which(names(sssnp2)==final_str_frags_plus_vb[(i)+(mismatch_pos[[select_candidate]][2])-1])]=0
                    SNP1=sssnp1>0
                    SNP2=sssnp2>0

                    type_if_SNP[[select_candidate]]=c(
                      paste(strsplit(toString(names(SNP1[SNP1==TRUE]),),split = ", " )[[1]],sep = "",collapse = " & "),
                      paste(strsplit(toString(names(SNP2[SNP2==TRUE]),),split = ", " )[[1]],sep = "",collapse = " & "))


                  }}
                if (mismatch_num[[select_candidate]]==0) {
                  type_if_SNP[[select_candidate]]="No_SNP"
                  mismatch_pos[[select_candidate]]="-"
                }
                ###File name for plot result
                title_for_crRNA_plot_DNA=paste0("Information plot for ",Guide_RNA_ID[[select_candidate]], " crRNA candidate")
                filename2=paste0("ID-",paste0(ID_code,"_",ConservationMethod,"_",which(crRNA_seq==select_candidate)),".png",collapse = "")

                if (i<=20) {
                  setwd(plots_folder_name)
                  png(filename2, width = 1500,height = 866)
                  layout(matrix(c(2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  1,1,1,1,1,1,
                                  1,1,1,1,1,1,
                                  1,1,1,1,1,1), 8,6, byrow = TRUE))

                  crRNA_plot_DNA(final_str_frags_plus_vb[4:(i+W+20-1)],i,
                                 start_pos[[select_candidate]],
                                 mismatch_pos[[select_candidate]],
                                 type_if_SNP[[select_candidate]],
                                 title_for_crRNA_plot_DNA,W,
                                 paste(efficiency[[select_candidate]],"Efficiency"),
                                 strand[[select_candidate]],3)

                  my_Dot_Plot(select_candidate,
                              select_candidate,
                              points_col = "cyan3",
                              xlab = "Specer Seq",
                              ylab = "Specer Seq")

                  sec_structure_plot(select_candidate,
                                     spaser_structure$structure,
                                     first_backbone = "N",
                                     main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                     free_energy = spaser_structure$free_anergy_seq)

                  sec_structure_plot(backbone_spacer,
                                     backbone_spacer_structure$structure,
                                     first_backbone,
                                     back_name,
                                     length(backbone_structure),
                                     paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

                  dev.off()
                  setwd(first_dir)
                  link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                    mutate( path = paste0('file:', plots_folder_name,"/" , sort(file ))
                            , path = stringr::str_replace_all( path, ' ', '%20')
                            , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
                  link_plots[[select_candidate]]=link$link
                }
                if (i>20 & i<length(final_str_frags_plus_vb)-20-W+1-3) {
                  setwd(plots_folder_name)
                  png(filename2, width = 1500,height = 866)
                  layout(matrix(c(2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  1,1,1,1,1,1,
                                  1,1,1,1,1,1,
                                  1,1,1,1,1,1), 8,6, byrow = TRUE))

                  crRNA_plot_DNA(final_str_frags_plus_vb[(i-20):(i+W+20-1)],21,
                                 start_pos[[select_candidate]],
                                 mismatch_pos[[select_candidate]],
                                 type_if_SNP[[select_candidate]],
                                 title_for_crRNA_plot_DNA,W,
                                 paste(efficiency[[select_candidate]],"Efficiency"),
                                 strand[[select_candidate]],3)

                  my_Dot_Plot(select_candidate,
                              select_candidate,
                              points_col = "cyan3",
                              xlab = "Specer Seq",
                              ylab = "Specer Seq")

                  sec_structure_plot(select_candidate,
                                     spaser_structure$structure,
                                     first_backbone = "N",
                                     main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                     free_energy = spaser_structure$free_anergy_seq)

                  sec_structure_plot(backbone_spacer,
                                     backbone_spacer_structure$structure,
                                     first_backbone,
                                     back_name,
                                     length(backbone_structure),
                                     paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

                  dev.off()
                  setwd(first_dir)
                  link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                    mutate( path = paste0('file:', plots_folder_name,"/" , sort(file ))
                            , path = stringr::str_replace_all( path, ' ', '%20')
                            , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
                  link_plots[[select_candidate]]=link$link
                }
                if (i>=length(final_str_frags_plus_vb)-20-W-3) {
                  setwd(plots_folder_name)
                  png(filename2, width = 1500,height = 866)
                  layout(matrix(c(2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  2,2,3,3,4,4,
                                  1,1,1,1,1,1,
                                  1,1,1,1,1,1,
                                  1,1,1,1,1,1), 8,6, byrow = TRUE))
                  crRNA_plot_DNA(final_str_frags_plus_vb[(i-20):(length(final_str_frags_plus_vb)-3)],21,
                                 start_pos[[select_candidate]],
                                 mismatch_pos[[select_candidate]],
                                 type_if_SNP[[select_candidate]],
                                 title_for_crRNA_plot_DNA,W,
                                 paste(efficiency[[select_candidate]],"Efficiency"),
                                 strand[[select_candidate]],3)

                  my_Dot_Plot(select_candidate,
                              select_candidate,
                              points_col = "cyan3",
                              xlab = "Specer Seq",
                              ylab = "Specer Seq")

                  sec_structure_plot(select_candidate,
                                     spaser_structure$structure,
                                     first_backbone = "N",
                                     main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                     free_energy = spaser_structure$free_anergy_seq)

                  sec_structure_plot(backbone_spacer,
                                     backbone_spacer_structure$structure,
                                     first_backbone,
                                     back_name,
                                     length(backbone_structure),
                                     paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

                  dev.off()
                  setwd(first_dir)
                  link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                    mutate( path = paste0('file:', plots_folder_name,"/" , sort(file ))
                            , path = stringr::str_replace_all( path, ' ', '%20')
                            , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
                  link_plots[[select_candidate]]=link$link
                }
                score_Gaide[[select_candidate]]=sum(c((75-(mismatch_num[[select_candidate]]*5))
                                                      ,(10*(1-(final_free_anergy/max_free_energy)))))
              }


            }
          }

        }

        blast_query=c(prtospasertta_pos,
                      prtospaserttg_pos,
                      prtospaserttt_pos,
                      prtospaserttc_pos,
                      prtospasertta_neg,
                      prtospaserttg_neg,
                      prtospaserttt_neg,
                      prtospaserttc_neg)

      }
      #########################################################################









      ############################# for casV_A type ###########################
      if (CRISPR_Type=="casV_A") {
        final_str_frags_plus_va=c("l","l","l","l",final_str_frags,"l","l","l","l")
        prtospaserttta_pos=list()
        prtospasertttg_pos=list()
        prtospasertttc_pos=list()

        prtospaserttta_neg=list()
        prtospasertttg_neg=list()
        prtospasertttc_neg=list()

        for (i in 5:(length(final_str_frags_plus_va)-W+1-8)) {
          candidate=final_str_frags_plus_va[i:(i+W-1)]
          candidate_scores1=final_score1[(i-4):(i+W-1-4)]
          GC_condidate=round(length(which(candidate=="g"| candidate=="c"))/length(candidate)*100,digits = 2)

          if (final_str_frags_plus_va[i+W]    != "a" &
              final_str_frags_plus_va[i+W+1]  == "a" &
              final_str_frags_plus_va[i+W+2]  == "a" &
              final_str_frags_plus_va[i+W+3]  == "a") {
            select_candidate=paste(revers_compliment(candidate),collapse = "")



            candidate_scores=revers_function(candidate_scores1)

            if(mismatch_del(candidate_scores,start,end,seed_start,seed_end)!=TRUE) {

              backbone_spacer=paste0(backbone,select_candidate,collapse = "")
              backbone_spacer_structure=pridict_structure(backbone_spacer,RNAfold_location = RNAfold_path)
              spaser_structure=pridict_structure(select_candidate,RNAfold_path)

              if (paste(backbone_spacer_structure$structure[1:length(backbone_structure)],collapse = "")==
                  paste(backbone_structure,collapse = "")){
                final_free_anergy=backbone_spacer_structure$free_anergy_seq-backbone_free_enargy
                crRNA_free_energy[[select_candidate]]=paste(round(final_free_anergy,digits = 4),"Natural structure of DR")
                back_name=paste("Direct Repeat (Natural structure) ")
              }

              if(paste(backbone_spacer_structure$structure[1:length(backbone_structure)],collapse = "")!=
                 paste(backbone_structure,collapse = "")){
                final_free_anergy=backbone_spacer_structure$free_anergy_seq
                crRNA_free_energy[[select_candidate]]=paste(round(final_free_anergy,digits = 4),"Unnatural structure of DR")
                back_name=paste("Direct Repeat (Unnatural structure) ")
              }

              prtospaser[[select_candidate]]=revers_compliment(candidate)
              vv[[select_candidate]]=round((candidate_scores),digits = 4)
              GC_content[[select_candidate]]=paste0(GC_condidate,"%")
              start_pos[[select_candidate]]=(length(final_str_frags_plus_va)-8)-((i-4)+W-1)+1
              mismatch_pos[[select_candidate]]=which((candidate_scores)<ConservationThreshold)
              mismatch_num[[select_candidate]]=length(which(candidate_scores<ConservationThreshold))
              rank[[select_candidate]]=length(which(candidate_scores<ConservationThreshold))+1
              efficiency[[select_candidate]]=mean(candidate_scores)
              strand[[select_candidate]]="-"
              crRNA_seq=c(crRNA_seq,select_candidate)
              Guide_RNA_ID[[select_candidate]]=paste0(ID_code,"_",ConservationMethod,"_",which(crRNA_seq==select_candidate))

              prtospaserttta_neg[[select_candidate]]=c("t","t","t","a",revers_compliment(candidate))
              prtospasertttg_neg[[select_candidate]]=c("t","t","t","g",revers_compliment(candidate))
              prtospasertttc_neg[[select_candidate]]=c("t","t","t","c",revers_compliment(candidate))


              if (mismatch_num[[select_candidate]]>0){
                if (mismatch_num[[select_candidate]]==1) {

                  sssnp=normal_score_matrix1[,original_pos_in_matrix[(i-4)+(W-mismatch_pos[[select_candidate]])]]
                  sssnp[which(names(sssnp)==final_str_frags_plus_va[(i)+(W-mismatch_pos[[select_candidate]])])]=0
                  SNP=sssnp>0
                  type_if_SNP[[select_candidate]]=
                    paste(revers_compliment(strsplit(toString(names(SNP[SNP==TRUE]),),split = ", ")[[1]]),collapse = " & ",sep="")

                }
                else if (mismatch_num[[select_candidate]]==2) {
                  type_if_SNP[[select_candidate]]=c()

                  sssnp1=normal_score_matrix1[,original_pos_in_matrix[(i-4)+(W-mismatch_pos[[select_candidate]][1])]]
                  sssnp1[which(names(sssnp1)==final_str_frags_plus_va[(i)+(W-mismatch_pos[[select_candidate]][1])])]=0
                  sssnp2=normal_score_matrix1[,original_pos_in_matrix[(i-4)+(W-mismatch_pos[[select_candidate]][2])]]
                  sssnp2[which(names(sssnp2)==final_str_frags_plus_va[(i)+(W-mismatch_pos[[select_candidate]][2])])]=0

                  SNP1=sssnp1>0
                  SNP2=sssnp2>0

                  type_if_SNP[[select_candidate]]=c(
                    paste(revers_compliment(strsplit(toString(names(SNP1[SNP1==TRUE]),),split = ", ")[[1]]),collapse = " & ",sep=""),
                    paste(revers_compliment(strsplit(toString(names(SNP2[SNP2==TRUE]),),split = ", ")[[1]]),collapse = " & ",sep=""))



                }}

              if (mismatch_num[[select_candidate]]==0) {
                type_if_SNP[[select_candidate]]="No_SNP"
                mismatch_pos[[select_candidate]]="-"
              }
              ###File name for plot result
              title_for_crRNA_plot_DNA=paste0("Information plot for ",Guide_RNA_ID[[select_candidate]], " crRNA candidate")
              title_for_crRNA_Dotplot=paste0("Dotplot for ", Guide_RNA_ID[[select_candidate]]," crRNA candidate")
              filename2=paste0("ID-",paste0(ID_code,"_",ConservationMethod,"_",which(crRNA_seq==select_candidate)),".png",collapse = "")

              if (i<=20) {
                setwd(plots_folder_name)
                png(filename2, width = 1500,height = 866)
                layout(matrix(c(2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                1,1,1,1,1,1,
                                1,1,1,1,1,1,
                                1,1,1,1,1,1), 8,6, byrow = TRUE))

                crRNA_plot_DNA(final_str_frags_plus_va[4:(i+W+20-1)],i,
                               start_pos[[select_candidate]],
                               mismatch_pos[[select_candidate]],
                               type_if_SNP[[select_candidate]],
                               title_for_crRNA_plot_DNA,W,
                               paste(efficiency[[select_candidate]],"Efficiency"),
                               strand[[select_candidate]],4)

                my_Dot_Plot(select_candidate,
                            select_candidate,
                            points_col = "cyan3",
                            xlab = "Specer Seq",
                            ylab = "Specer Seq")

                sec_structure_plot(select_candidate,
                                   spaser_structure$structure,
                                   first_backbone = "N",
                                   main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                   free_energy = spaser_structure$free_anergy_seq)

                sec_structure_plot(backbone_spacer,
                                   backbone_spacer_structure$structure,
                                   first_backbone,
                                   back_name,
                                   length(backbone_structure),
                                   paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

                dev.off()
                setwd(first_dir)
                link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                  mutate( path = paste0('file:', plots_folder_name,"/",  sort(file ))
                          , path = stringr::str_replace_all( path, ' ', '%20')
                          , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
                link_plots[[select_candidate]]=link$link
              }
              if (i>20 & i<length(final_str_frags_plus_va)-20-W+1-4) {
                setwd(plots_folder_name)
                png(filename2, width = 1500,height = 866)
                layout(matrix(c(2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                1,1,1,1,1,1,
                                1,1,1,1,1,1,
                                1,1,1,1,1,1), 8,6, byrow = TRUE))

                crRNA_plot_DNA(final_str_frags_plus_va[(i-20):(i+W+20-1)],21,
                               start_pos[[select_candidate]],
                               mismatch_pos[[select_candidate]],
                               type_if_SNP[[select_candidate]],
                               title_for_crRNA_plot_DNA,W,
                               paste(efficiency[[select_candidate]],"Efficiency"),
                               strand[[select_candidate]],4)

                my_Dot_Plot(select_candidate,
                            select_candidate,
                            points_col = "cyan3",
                            xlab = "Specer Seq",
                            ylab = "Specer Seq")

                sec_structure_plot(select_candidate,
                                   spaser_structure$structure,
                                   first_backbone = "N",
                                   main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                   free_energy = spaser_structure$free_anergy_seq)

                sec_structure_plot(backbone_spacer,
                                   backbone_spacer_structure$structure,
                                   first_backbone,
                                   back_name,
                                   length(backbone_structure),
                                   paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))
                dev.off()

                setwd(first_dir)
                link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                  mutate( path = paste0('file:', plots_folder_name,"/",  sort(file ))
                          , path = stringr::str_replace_all( path, ' ', '%20')
                          , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
                link_plots[[select_candidate]]=link$link
              }
              if (i>=length(final_str_frags_plus_va)-20-W-4) {
                setwd(plots_folder_name)
                png(filename2, width = 1200,height = 800)
                png(filename2, width = 1500,height = 866)
                layout(matrix(c(2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                1,1,1,1,1,1,
                                1,1,1,1,1,1,
                                1,1,1,1,1,1), 8,6, byrow = TRUE))

                crRNA_plot_DNA(final_str_frags_plus_va[(i-20):(length(final_str_frags_plus_va)-4)],21,
                               start_pos[[select_candidate]],
                               mismatch_pos[[select_candidate]],
                               type_if_SNP[[select_candidate]],
                               title_for_crRNA_plot_DNA,W,
                               paste(efficiency[[select_candidate]],"Efficiency"),
                               strand[[select_candidate]],4)

                my_Dot_Plot(select_candidate,
                            select_candidate,
                            points_col = "cyan3",
                            xlab = "Specer Seq",
                            ylab = "Specer Seq")

                sec_structure_plot(select_candidate,
                                   spaser_structure$structure,
                                   first_backbone = "N",
                                   main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                   free_energy = spaser_structure$free_anergy_seq)

                sec_structure_plot(backbone_spacer,
                                   backbone_spacer_structure$structure,
                                   first_backbone,
                                   back_name,
                                   length(backbone_structure),
                                   paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

                dev.off()
                setwd(first_dir)
                link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                  mutate( path = paste0('file:', plots_folder_name,"/",  sort(file ))
                          , path = stringr::str_replace_all( path, ' ', '%20')
                          , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
                link_plots[[select_candidate]]=link$link
              }
              score_Gaide[[select_candidate]]=sum(c((75-(mismatch_num[[select_candidate]]*5))
                                                    ,(10*(1-(final_free_anergy/max_free_energy)))))

            }}






          if (final_str_frags_plus_va[i-4] == "t" &
              final_str_frags_plus_va[i-3] == "t" &
              final_str_frags_plus_va[i-2] == "t" &
              final_str_frags_plus_va[i-1] != "t") {
            select_candidate=paste(candidate,collapse = "")
            candidate_scores=candidate_scores1
            if(mismatch_del(candidate_scores,start,end,seed_start,seed_end)!=TRUE){

              backbone_spacer=paste0(backbone,select_candidate,collapse = "")
              backbone_spacer_structure=pridict_structure(backbone_spacer,RNAfold_location = RNAfold_path)
              spaser_structure=pridict_structure(select_candidate,RNAfold_path)

              if (paste(backbone_spacer_structure$structure[1:length(backbone_structure)],collapse = "")==
                  paste(backbone_structure,collapse = "")){
                final_free_anergy=backbone_spacer_structure$free_anergy_seq-backbone_free_enargy
                crRNA_free_energy[[select_candidate]]=paste(round(final_free_anergy,digits = 4),"Natural structure of DR")
                back_name=paste("Direct Repeat (Natural structure) ")
              }

              if(paste(backbone_spacer_structure$structure[1:length(backbone_structure)],collapse = "")!=
                 paste(backbone_structure,collapse = "")){
                final_free_anergy=backbone_spacer_structure$free_anergy_seq
                crRNA_free_energy[[select_candidate]]=paste(round(final_free_anergy,digits = 4),"Unnatural structure of DR")
                back_name=paste("Direct Repeat (Unnatural structure) ")
              }


              prtospaser[[select_candidate]]=candidate
              vv[[select_candidate]]=round((candidate_scores),digits = 4)
              GC_content[[select_candidate]]=paste0(GC_condidate,"%")
              start_pos[[select_candidate]]=i-4
              mismatch_pos[[select_candidate]]=which((candidate_scores)<ConservationThreshold)
              mismatch_num[[select_candidate]]=length(which(candidate_scores<ConservationThreshold))
              rank[[select_candidate]]=mismatch_num[[select_candidate]]+1
              efficiency[[select_candidate]]=mean(candidate_scores)
              strand[[select_candidate]]="+"
              crRNA_seq=c(crRNA_seq,select_candidate)
              Guide_RNA_ID[[select_candidate]]=paste0(ID_code,"_",ConservationMethod,"_",which(crRNA_seq==select_candidate))

              prtospaserttta_pos[[select_candidate]]=c("t","t","t","a",candidate)
              prtospasertttg_pos[[select_candidate]]=c("t","t","t","g",candidate)
              prtospasertttc_pos[[select_candidate]]=c("t","t","t","c",candidate)



              if (mismatch_num[[select_candidate]]>0){
                if (mismatch_num[[select_candidate]]==1) {

                  sssnp=normal_score_matrix1[,original_pos_in_matrix[(i-4)+(mismatch_pos[[select_candidate]]-1)]]
                  sssnp[which(names(sssnp)==final_str_frags_plus_va[(i)+(mismatch_pos[[select_candidate]]-1)])]=0
                  SNP=sssnp>0
                  type_if_SNP[[select_candidate]]=
                    paste(strsplit(toString(names(SNP[SNP==TRUE]),),split = ", ")[[1]],collapse = " & ",sep="")


                }
                else if (mismatch_num[[select_candidate]]==2) {
                  type_if_SNP[[select_candidate]]=c()

                  sssnp1=normal_score_matrix1[,original_pos_in_matrix[(i-4)+(mismatch_pos[[select_candidate]][1]-1)]]
                  sssnp1[which(names(sssnp1)==final_str_frags_plus_va[(i)+(mismatch_pos[[select_candidate]][1]-1)])]=0
                  sssnp2=normal_score_matrix1[,original_pos_in_matrix[(i-4)+(mismatch_pos[[select_candidate]][2]-1)]]
                  sssnp2[which(names(sssnp2)==final_str_frags_plus_va[(i)+(mismatch_pos[[select_candidate]][2])-1])]=0

                  SNP1=sssnp1>0
                  SNP2=sssnp2>0

                  type_if_SNP[[select_candidate]]=c(
                    paste(strsplit(toString(names(SNP1[SNP1==TRUE]),),split = ", ")[[1]],sep = "" ,collapse = " & "),
                    paste(strsplit(toString(names(SNP2[SNP2==TRUE]),),split = ", ")[[1]],collapse = " & ",sep=""))


                }}

              if (mismatch_num[[select_candidate]]==0) {
                type_if_SNP[[select_candidate]]="No_SNP"
                mismatch_pos[[select_candidate]]="-"
              }
              ###File name for plot result
              title_for_crRNA_plot_DNA=paste0("Information plot for ",Guide_RNA_ID[[select_candidate]], " crRNA candidate")
              title_for_crRNA_Dotplot=paste0("Dotplot for ", Guide_RNA_ID[[select_candidate]]," crRNA candidate")
              filename2=paste0("ID-",paste0(ID_code,"_",ConservationMethod,"_",which(crRNA_seq==select_candidate)),".png",collapse = "")

              if (i<=20) {
                setwd(plots_folder_name)
                png(filename2, width = 1500,height = 866)
                layout(matrix(c(2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                1,1,1,1,1,1,
                                1,1,1,1,1,1,
                                1,1,1,1,1,1), 8,6, byrow = TRUE))

                crRNA_plot_DNA(final_str_frags_plus_va[4:(i+W+20-1)],i,
                               start_pos[[select_candidate]],
                               mismatch_pos[[select_candidate]],
                               type_if_SNP[[select_candidate]],
                               title_for_crRNA_plot_DNA,W,
                               paste(efficiency[[select_candidate]],"Efficiency"),
                               strand[[select_candidate]],4)

                my_Dot_Plot(select_candidate,
                            select_candidate,
                            points_col = "cyan3",
                            xlab = "Specer Seq",
                            ylab = "Specer Seq")

                sec_structure_plot(select_candidate,
                                   spaser_structure$structure,
                                   first_backbone = "N",
                                   main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                   free_energy = spaser_structure$free_anergy_seq)

                sec_structure_plot(backbone_spacer,
                                   backbone_spacer_structure$structure,
                                   first_backbone,
                                   back_name,
                                   length(backbone_structure),
                                   paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

                dev.off()
                setwd(first_dir)
                link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                  mutate( path = paste0('file:', plots_folder_name,"/",  sort(file ))
                          , path = stringr::str_replace_all( path, ' ', '%20')
                          , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
                link_plots[[select_candidate]]=link$link
              }
              if (i>20 & i<length(final_str_frags_plus_va)-20-W+1-4) {
                setwd(plots_folder_name)
                png(filename2, width = 1500,height = 866)
                layout(matrix(c(2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                1,1,1,1,1,1,
                                1,1,1,1,1,1,
                                1,1,1,1,1,1), 8,6, byrow = TRUE))
                crRNA_plot_DNA(final_str_frags_plus_va[(i-20):(i+W+20-1)],21,
                               start_pos[[select_candidate]],
                               mismatch_pos[[select_candidate]],
                               type_if_SNP[[select_candidate]],
                               title_for_crRNA_plot_DNA,W,
                               paste(efficiency[[select_candidate]],"Efficiency"),
                               strand[[select_candidate]],4)

                my_Dot_Plot(select_candidate,
                            select_candidate,
                            points_col = "cyan3",
                            xlab = "Specer Seq",
                            ylab = "Specer Seq")

                sec_structure_plot(select_candidate,
                                   spaser_structure$structure,
                                   first_backbone = "N",
                                   main="Spacer self complementary",
                                   free_energy = spaser_structure$free_anergy_seq)

                sec_structure_plot(backbone_spacer,
                                   backbone_spacer_structure$structure,
                                   first_backbone,
                                   back_name,
                                   length(backbone_structure),"crRNA self complementary")


                dev.off()
                setwd(first_dir)
                link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                  mutate( path =  paste0('file:', plots_folder_name,"/",  sort(file ))
                          , path = stringr::str_replace_all( path, ' ', '%20')
                          , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
                link_plots[[select_candidate]]=link$link
              }
              if (i>=length(final_str_frags_plus_va)-20-W-4) {
                setwd(plots_folder_name)
                png(filename2, width = 1500,height = 866)
                layout(matrix(c(2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                2,2,3,3,4,4,
                                1,1,1,1,1,1,
                                1,1,1,1,1,1,
                                1,1,1,1,1,1), 8,6, byrow = TRUE))

                crRNA_plot_DNA(final_str_frags_plus_va[(i-20):(length(final_str_frags_plus_va)-4)],21,
                               start_pos[[select_candidate]],
                               mismatch_pos[[select_candidate]],
                               type_if_SNP[[select_candidate]],
                               title_for_crRNA_plot_DNA,W,
                               paste(efficiency[[select_candidate]],"Efficiency"),
                               strand[[select_candidate]],4)

                my_Dot_Plot(select_candidate,
                            select_candidate,
                            points_col = "cyan3",
                            xlab = "Specer Seq",
                            ylab = "Specer Seq")

                sec_structure_plot(select_candidate,
                                   spaser_structure$structure,
                                   first_backbone = "N",
                                   main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                   free_energy = spaser_structure$free_anergy_seq)

                sec_structure_plot(backbone_spacer,
                                   backbone_spacer_structure$structure,
                                   first_backbone,
                                   back_name,
                                   length(backbone_structure),
                                   paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

                dev.off()
                setwd(first_dir)
                link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                  mutate( path =  paste0('file:', plots_folder_name,"/",  sort(file ))
                          , path = stringr::str_replace_all( path, ' ', '%20')
                          , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
                link_plots[[select_candidate]]=link$link
              }

              score_Gaide[[select_candidate]]=sum(c((75-(mismatch_num[[select_candidate]]*5))
                                                    ,(10*(1-(final_free_anergy/max_free_energy)))))
            }

          }
        }


        blast_query=c(prtospaserttta_pos,
                      prtospasertttg_pos,
                      prtospasertttc_pos,

                      prtospaserttta_neg,
                      prtospasertttg_neg,
                      prtospasertttc_neg)

      }
      #########################################################################










      ############################# for casV_F1 type ##########################
      if (CRISPR_Type=="casV_F1") {

        for (i in 1:(length(final_str_frags)-W+1)) {
          candidate=final_str_frags[i:(i+W-1)]
          candidate_scores=c(final_score1[i:(i+W-1)])
          GC_condidate=round(length(which(candidate=="g"| candidate=="c"))/length(candidate)*100,digits = 2)

          if (mismatch_del(candidate_scores,start,end,seed_start,seed_end)!=TRUE) {
            select_candidate=paste(revers_compliment(candidate),collapse = "")

            backbone_spacer=paste0(backbone,select_candidate,collapse = "")
            backbone_spacer_structure=pridict_structure(backbone_spacer,RNAfold_location = RNAfold_path)
            spaser_structure=pridict_structure(select_candidate,RNAfold_path)

            if (paste(backbone_spacer_structure$structure[1:length(backbone_structure)],collapse = "")==
                paste(backbone_structure,collapse = "")){
              final_free_anergy=backbone_spacer_structure$free_anergy_seq-backbone_free_enargy
              crRNA_free_energy[[select_candidate]]=paste(round(final_free_anergy,digits = 4),"Natural structure of DR")
              back_name=paste("Direct Repeat (Natural structure) ")
            }

            if(paste(backbone_spacer_structure$structure[1:length(backbone_structure)],collapse = "")!=
               paste(backbone_structure,collapse = "")){
              final_free_anergy=backbone_spacer_structure$free_anergy_seq
              crRNA_free_energy[[select_candidate]]=paste(round(final_free_anergy,digits = 4),"Unnatural structure of DR")
              back_name=paste("Direct Repeat (Unnatural structure) ")
            }



            prtospaser[[select_candidate]]=candidate
            GC_content[[select_candidate]]=paste0(GC_condidate,"%")
            crRNA_seq=c(crRNA_seq,select_candidate)
            vv[[select_candidate]]=round(revers_function(candidate_scores),digits = 4)
            start_pos[[select_candidate]]=i
            mismatch_pos[[select_candidate]]=which(revers_function(candidate_scores)<ConservationThreshold)
            mismatch_num[[select_candidate]]=length(which(candidate_scores<ConservationThreshold))
            rank[[select_candidate]]=length(which(candidate_scores<ConservationThreshold))+1
            efficiency[[select_candidate]]=mean(candidate_scores)
            Guide_RNA_ID[[select_candidate]]=paste0(ID_code,"_",ConservationMethod,"_",which(crRNA_seq==select_candidate))
            pos_spacer_in_sec_stracture=i

            if (length(final_str_frags)>2000){
              if (i<=1000) {
                RNAfols_structure_format=RNAfols_structure_format_start
                pos_spacer_in_sec_stracture=i
              }
              else if (i>1000 & i<length(final_str_frags)-1000+1) {

                write.fasta(final_str_frags[(i-1000):(i+1000-1)],"final_str_frags","final_str_frags.fasta")
                structure_of_RNAfold=system(paste(RNAfold_path,"  ","final_str_frags.fasta" ),intern=T)
                file.remove("final_str_frags.fasta")
                file.remove("final_str_frags_ss.ps")
                RNAfols_structure_format=strsplit(structure_of_RNAfold[3],split = "")[[1]][1:2000]

                RNAplot_file_name=paste0(">","final_structur_",(i-1000),"_to_",(i+1000-1))
                fileConn<-file("sec_stracuter.txt")
                writeLines(c(RNAplot_file_name,paste(final_str_frags[(i-1000):(i+1000-1)],collapse = ""),
                             paste(RNAfols_structure_format_start,collapse = "")), fileConn)

                close(fileConn)

                if (host_system=="Windows") {
                  RNAPlot_command_line=paste(RNAplot_location," -o  svg ", " sec_stracuter.txt")
                  lapply(RNAPlot_command_line, system)
                  file.remove("sec_stracuter.txt")
                  file.copy(paste0("final_structur_",(i-1000),"_to_",(i+1000-1),"_ss.svg"),plots_folder_name)
                  file.remove(paste0("final_structur_",(i-1000),"_to_",(i+1000-1),"_ss.svg"))
                }
                if (host_system=="Linux") {
                  RNAPlot_command_line=paste(RNAplot_location, " sec_stracuter.txt ")
                  lapply(RNAPlot_command_line, system)
                  file.remove("sec_stracuter.txt")
                  file.copy(paste0("final_structur_",(i-1000),"_to_",(i+1000-1),"_ss.ps"),plots_folder_name)
                  file.remove(paste0("final_structur_",(i-1000),"_to_",(i+1000-1),"_ss.ps"))
                }


                pos_spacer_in_sec_stracture=1000
              }
              else if (i>=length(final_str_frags)-1000) {
                RNAfols_structure_format=RNAfols_structure_format_end
                pos_spacer_in_sec_stracture=2000-(length(final_str_frags)-i)
              }
            }



            if (i<=100) {
              a_side=which(RNAfols_structure_format[pos_spacer_in_sec_stracture:(pos_spacer_in_sec_stracture+W-1)]==".")
              a_flank_side=which(RNAfols_structure_format[c((1:pos_spacer_in_sec_stracture-1),
                                                            (pos_spacer_in_sec_stracture+W):(pos_spacer_in_sec_stracture+W+100-1))]==".")
              accessible_side[[select_candidate]]=round(length(a_side)/W, digits=4)
              accessible_flank_side[[select_candidate]]=
                round(length(a_flank_side)/length(c((1:pos_spacer_in_sec_stracture-1),
                                                    (pos_spacer_in_sec_stracture+W):(pos_spacer_in_sec_stracture+W+100-1))), digits=4)

            }
            else if (i>100 & i<length(final_str_frags)-100-W+1) {
              a_side=which(RNAfols_structure_format[pos_spacer_in_sec_stracture:(pos_spacer_in_sec_stracture+W-1)]==".")
              a_flank_side=which(RNAfols_structure_format[c((pos_spacer_in_sec_stracture-100):(pos_spacer_in_sec_stracture-1),
                                                            (pos_spacer_in_sec_stracture+W):(pos_spacer_in_sec_stracture+W+100-1))]==".")
              accessible_side[[select_candidate]]=round(length(a_side)/W, digits=4)
              accessible_flank_side[[select_candidate]]=
                round(length(a_flank_side)/length(c((pos_spacer_in_sec_stracture-100):(pos_spacer_in_sec_stracture-1),
                                                    (pos_spacer_in_sec_stracture+W):(pos_spacer_in_sec_stracture+W+100-1))), digits=4)
            }
            else if (i>=length(final_str_frags)-100-W) {
              a_side=which(RNAfols_structure_format[pos_spacer_in_sec_stracture:(pos_spacer_in_sec_stracture+W-1)]==".")
              a_flank_side=which(RNAfols_structure_format[c((pos_spacer_in_sec_stracture-100):(pos_spacer_in_sec_stracture-1),
                                                            (pos_spacer_in_sec_stracture+W):2000)]==".")
              accessible_side[[select_candidate]]=round(length(a_side)/W, digits=4)
              accessible_flank_side[[select_candidate]]=
                round(length(a_flank_side)/length(c((pos_spacer_in_sec_stracture-100):(pos_spacer_in_sec_stracture-1),
                                                    (pos_spacer_in_sec_stracture+W):2000)), digits=4)

            }


            if (mismatch_num[[select_candidate]]>0){
              if (mismatch_num[[select_candidate]]==1) {
                sssnp=normal_score_matrix1[,original_pos_in_matrix[i+(W-mismatch_pos[[select_candidate]])]]
                sssnp[which(names(sssnp)==final_str_frags[i+(W-mismatch_pos[[select_candidate]])])]=0
                SNP=sssnp>0
                type_if_SNP[[select_candidate]]=
                  paste(strsplit(toString(names(SNP[SNP==TRUE]),),split = ", ")[[1]],collapse = " & ",sep="")
              }
              else if (mismatch_num[[select_candidate]]==2) {
                type_if_SNP[[select_candidate]]=c()
                sssnp1=normal_score_matrix1[,original_pos_in_matrix[i+(W-mismatch_pos[[select_candidate]][1])]]
                sssnp1[which(names(sssnp1)==final_str_frags[i+(W-mismatch_pos[[select_candidate]][1])])]=0
                sssnp2=normal_score_matrix1[,original_pos_in_matrix[i+(W-mismatch_pos[[select_candidate]][2])]]
                sssnp2[which(names(sssnp2)==final_str_frags[i+(W-mismatch_pos[[select_candidate]][2])])]=0

                SNP1=sssnp1>0
                SNP2=sssnp2>0
                type_if_SNP[[select_candidate]]=c(
                  paste(strsplit(toString(names(SNP1[SNP1==TRUE]),),split = ", ")[[1]],collapse = " & ",sep=""),
                  paste(strsplit(toString(names(SNP2[SNP2==TRUE]),),split = ", ")[[1]],collapse = " & ",sep=""))

              }}
            if (mismatch_num[[select_candidate]]==0) {
              type_if_SNP[[select_candidate]]="No_SNP"
              mismatch_pos[[select_candidate]]="-"
            }
            ###File name for plot result
            title_for_crRNA_plot=paste0("Information plot for ",Guide_RNA_ID[[select_candidate]], " crRNA candidate")
            filename2=paste0("ID-",paste0(ID_code,"_",ConservationMethod,"_",which(crRNA_seq==select_candidate)),".png",collapse = "")

            if (i<=20) {
              setwd(plots_folder_name)
              png(filename2, width = 1300,height = 750)
              layout(matrix(c(2,2,3,3,4,4,2,2,3,3,4,4,1,1,1,1,1,1), 3, 6, byrow = TRUE))

              crRNA_plot(final_str_frags[1:(i+W+20-1)],
                         i,
                         start_pos[[select_candidate]],
                         mismatch_pos[[select_candidate]],
                         type_if_SNP[[select_candidate]],
                         title_for_crRNA_plot,W,"",paste(accessible_flank_side[[select_candidate]]),
                         paste(efficiency[[select_candidate]],"efficiency"),
                         paste(accessible_side[[select_candidate]]))

              my_Dot_Plot(select_candidate,
                          select_candidate,
                          points_col = "cyan3",
                          xlab = "Specer Seq",
                          ylab = "Specer Seq")

              sec_structure_plot(select_candidate,
                                 spaser_structure$structure,
                                 first_backbone = "N",
                                 main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                 free_energy = spaser_structure$free_anergy_seq)

              sec_structure_plot(backbone_spacer,
                                 backbone_spacer_structure$structure,
                                 first_backbone,
                                 back_name,
                                 length(backbone_structure),
                                 paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

              dev.off()
              setwd(first_dir)
              link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                mutate( path = paste0('file:',  plots_folder_name,"/",  sort(file ))
                        , path = stringr::str_replace_all( path, ' ', '%20')
                        , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
              link_plots[[select_candidate]]=link$link
            }
            if (i>20 & length(final_str_frags)-20-W+1) {
              setwd(plots_folder_name)
              png(filename2, width = 1300,height = 750)
              layout(matrix(c(2,2,3,3,4,4,2,2,3,3,4,4,1,1,1,1,1,1), 3, 6, byrow = TRUE))
              crRNA_plot(final_str_frags[(i-20):(i+W+20-1)],21,
                         start_pos[[select_candidate]],
                         mismatch_pos[[select_candidate]],
                         type_if_SNP[[select_candidate]],
                         title_for_crRNA_plot,W,"",
                         paste(accessible_flank_side[[select_candidate]]),paste(efficiency[[select_candidate]],"efficiency"),
                         paste(accessible_side[[select_candidate]]))

              my_Dot_Plot(select_candidate,
                          select_candidate,
                          points_col = "cyan3",
                          xlab = "Specer Seq",
                          ylab = "Specer Seq")

              sec_structure_plot(select_candidate,
                                 spaser_structure$structure,
                                 first_backbone = "N",
                                 main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                 free_energy = spaser_structure$free_anergy_seq)

              sec_structure_plot(backbone_spacer,
                                 backbone_spacer_structure$structure,
                                 first_backbone,
                                 back_name,
                                 length(backbone_structure),
                                 paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

              dev.off()
              setwd(first_dir)
              link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                mutate( path = paste0('file:',  plots_folder_name,"/",  sort(file ))
                        , path = stringr::str_replace_all( path, ' ', '%20')
                        , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
              link_plots[[select_candidate]]=link$link
            }
            if (i>=length(final_str_frags)-20-W) {
              setwd(plots_folder_name)
              png(filename2, width = 1300,height = 750)
              layout(matrix(c(2,2,3,3,4,4,2,2,3,3,4,4,1,1,1,1,1,1), 3, 6, byrow = TRUE))

              crRNA_plot(final_str_frags[(i-20+1):length(final_str_frags)],21,
                         start_pos[[select_candidate]],mismatch_pos[[select_candidate]],
                         type_if_SNP[[select_candidate]],title_for_crRNA_plot,W,"",
                         paste(accessible_flank_side[[select_candidate]]),
                         paste(efficiency[[select_candidate]],"efficiency"),
                         paste(accessible_side[[select_candidate]]))

              my_Dot_Plot(select_candidate,
                          select_candidate,
                          points_col = "cyan3",
                          xlab = "Specer Seq",
                          ylab = "Specer Seq")

              sec_structure_plot(select_candidate,
                                 spaser_structure$structure,
                                 first_backbone = "N",
                                 main=paste0("Spacer Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"),
                                 free_energy = spaser_structure$free_anergy_seq)

              sec_structure_plot(backbone_spacer,
                                 backbone_spacer_structure$structure,
                                 first_backbone,
                                 back_name,
                                 length(backbone_structure),
                                 paste0("crRNA Self complementarity ","(ID=",Guide_RNA_ID[[select_candidate]],")"))

              dev.off()
              setwd(first_dir)
              link= tibble( file = dir( file.path( '.', plots_folder_name ), pattern=filename2 ) ) %>%
                mutate( path = paste0('file:',  plots_folder_name,"/",  sort(file ))
                        , path = stringr::str_replace_all( path, ' ', '%20')
                        , link =  paste0('<a  target=_blank href=', path, '>', file,'</a>' ) )
              link_plots[[select_candidate]]=link$link
            }

            score_Gaide[[select_candidate]]=sum(c((75-(mismatch_num[[select_candidate]]*5))
                                                  ,(2.5*accessible_side[[select_candidate]])
                                                  ,(2.5*accessible_flank_side[[select_candidate]])
                                                  ,(10*(1-(final_free_anergy/max_free_energy)))))


          }

        }
        blast_query=prtospaser
      }
      #########################################################################








      ############ replace sec plots to Secondary structure folder ############

      if (CRISPR_Type=="casVI_D"|CRISPR_Type=="casVI_A"|CRISPR_Type=="casVI_B"|CRISPR_Type=="casV_F1") {
        current_folder <- paste0(first_dir,"/Plots")
        new_folder <- paste0(first_dir,"/Secondary Structure")
        if (host_system=="Windows") {
          list_of_files <- list.files(current_folder, ".svg")
        }

        if (host_system=="Linux") {
          list_of_files <- list.files(current_folder, ".ps")
        }

        file.copy(file.path(current_folder,list_of_files), new_folder)
        file.remove(file.path(current_folder,list_of_files))

      }

      #########################################################################







      ################# View data (data.frame construction) ############

      if (length(crRNA_seq)>0) {
        data=as.data.frame(crRNA_seq)
        data$Guide_RNA_ID=Guide_RNA_ID
        names(data)[names(data) == "crRNA_seq"] <- "Guide_RNA"
        if (length(topoTypes_names)>1) {
          data$Major_Allele_Frequency=vv
        }

        data$Start=start_pos

        if (CRISPR_Type=="casV_B"|CRISPR_Type=="casV_A"){
          data$Strand=strand
        }
        data$GC_content=GC_content

        if (length(topoTypes_names)>1) {
          data$Mismatch_Number=mismatch_num
          data$Mismatch_Position=mismatch_pos
          data$Type_Of_SNP=type_if_SNP
          data$Efficiency=efficiency
          names(data)[names(data) == "Efficiency"] <- "Conservation Efficiency"
        }

        if (CRISPR_Type=="casVI_D"|CRISPR_Type=="casVI_A"|CRISPR_Type=="casVI_B"){
          data$local_U_rich=local_U_rich
          names(data)[names(data) == "local_U_rich"] <- local_rich_type
          data$Normalized_U_rich=normalized_U_rich
          names(data)[names(data) == "Normalized_U_rich"] <- Normal_rich_type
        }

        if (CRISPR_Type=="casVI_D"| CRISPR_Type=="casVI_A"|CRISPR_Type=="casVI_B"|CRISPR_Type=="casV_F1"){
          data$Protospacer_Accessibility_Score=accessible_side
          data$Local_Accessibility_Score=accessible_flank_side
        }

        data$Rank=rank
        data$Self_Complementarity_Free_Energy=crRNA_free_energy
        data$Score=unlist(score_Gaide)

        ###################################################################








        ###################### DT options #####################


        my.options <- list(autoWidth = T,
                           searching = T,
                           ordering = T,
                           filter = 'none',
                           lengthChange = T,
                           lengthMenu = c(5,10,50,100,1000),
                           paging = T,
                           info = T,
                           pageLength = 5,
                           searchHighlight = TRUE,
                           searchable=F)



        header.style <- "th { font-family: 'Times New Roman'; font-weight: bold; color: white ; background-color: #000000;}"


        ########################################################







        ###### checking internet connection and  blast with rules (off_target_prediction) ######

        setwd(first_dir)
        if (OffTarget==T) {
          if (LocalOff!=T & havingIP()==TRUE){
            if (is.null(Organism)!=T) {
              if (OffAsk) {
                min_approximate_blast_time=length(blast_query)*0.3*length(Organism)
                max_approximate_blast_time=length(blast_query)*1*length(Organism)
                question <- readline(paste0("off-Targets Prediction for ",ID_code," CRISPR-Cas System ","may be take long time",
                                            " (About ",
                                            min_approximate_blast_time," to ",
                                            max_approximate_blast_time,
                                            " minutes) ",
                                            "would you like to continue or no? (Y/N)"))
              }
              else if (OffAsk==F) {
                question="y"
              }
            }
            else{
              message('If you want to do online mode for off-Target Prediction,
                       please Ensure you introduced the "Organism" parameter.')
              question="n"
            }
          }

          if (LocalOff == T) {
            if (all(!(LocalFasta %in% c("path1", "path2")))) {
              if (length(LocalName) == length(LocalFasta)) {
                Organism = LocalName
              }
              else if (length(LocalName) != length(LocalFasta)) {
                Organism = paste0("User_Sequence_",
                                  1:length(LocalFasta))
              }
              if (OffAsk) {
                question <- readline("Depending on the size of your genomes, this step may be time consuming, Do you want to continue?(Y/N)")
              }
              else if (OffAsk == F) {
                question = "y"
              }
            }
            else {
              message("please Ensure you introduced the \"LocalFasta\" parameter",
                      "\n                 and check the \"LocalOff\" parameter be true.")
              question = "n"
            }
          }

          if(question== 'y'| question== 'Y'){

            Guide_RNA=data$Guide_RNA
            len_list=rep(0,length(Guide_RNA))
            count_mm0=as.list(len_list)
            count_mm1=as.list(len_list)
            count_mm2=as.list(len_list)
            count_mm3=as.list(len_list)
            count_mm4=as.list(len_list)

            data_off_target=as.data.frame(Guide_RNA)
            data_off_target$Guide_RNA_ID=data$Guide_RNA_ID
            for (org in Organism) {

              list_blast_mm_0=list()
              list_blast_mm_1=list()
              list_blast_mm_2=list()

              if (CRISPR_Type=="casV_B"|CRISPR_Type=="casV_A") {
                list_blast_mm_3=list()
                list_blast_mm_4=list()
              }
              write.fasta( blast_query,names(blast_query), file.out = "blast_query.FASTA")


              # online blast
              if (LocalOff!=T){
                entrez_query=paste(org,' [organism]',sep = "")

                entrez_query <- paste0('"',entrez_query,'"')

                blast_command=paste(Blast_location,
                                    '  -task blastn   -db ',  blast_database ,
                                    '    -entrez_query   ' , entrez_query,
                                    '     -word_size  7    -evalue  10    -query   blast_query.FASTA        -out  blast_out.txt    ' ,
                                    '  -strand  ', blast_strand,
                                    '    -outfmt "10 delim=@ qseqid sacc  qseq    sseq  length    qstart  qend sstart   send   gapopen   mismatch "    -num_alignments   5000        -remote',sep = '' ,collapse =  "" )
                lapply(blast_command, system)
              }

              #offline blast
              if (LocalOff==T){
                number_of_genome_file=which(org==Organism)
                file.copy(LocalFasta[number_of_genome_file],first_dir)
                neme_genome=basename(LocalFasta[number_of_genome_file])
                file.rename(neme_genome,"genome_name.fasta")
                makeblastdb_command=paste(  makeblastdb_location   ,
                                            '  -in genome_name.fasta  -out  genome_name_index  -dbtype nucl  -parse_seqids')
                lapply(makeblastdb_command, system)


                blast_command=paste(Blast_location,
                                    '  -task blastn   -db   genome_name_index   ' ,
                                    '  -word_size  7    -evalue  10    -query   blast_query.FASTA        -out  blast_out.txt    ' ,
                                    '  -strand  ', blast_strand,'  -num_threads  ', Threads,
                                    '     -outfmt "10 delim=@ qseqid sacc  qseq    sseq  length    qstart  qend sstart   send   gapopen   mismatch "    -num_alignments   5000 ',
                                    sep = '' ,collapse =  "" )
                lapply(blast_command, system)

                removes_file=list.files(first_dir,pattern = "genome_name*")
                file.remove(removes_file)

              }



              info_blast = file.info("blast_out.txt")


              if (info_blast$size!=0) {
                if (CRISPR_Type=="casVI_D"|CRISPR_Type=="casVI_A"|CRISPR_Type=="casV_F1"){
                  blast_result=read.table("blast_out.txt",sep = "@")
                  blast_result=rbind(blast_result[blast_result$V5==((W+length_pam_cas)-2) & blast_result$V11==0 & blast_result$V10==0,],
                                     blast_result[blast_result$V5==((W+length_pam_cas)-1) & blast_result$V11<=1 & blast_result$V10==0,],
                                     blast_result[blast_result$V5==(W+length_pam_cas) & blast_result$V11<=2 & blast_result$V10==0,])
                }


                else if (CRISPR_Type=="casVI_B"){
                  blast_result=read.table("blast_out.txt",sep = "@")
                  blast_result=rbind(blast_result[blast_result$V5==((W+length_pam_cas)-2) & blast_result$V6==1 & blast_result$V11==0 & blast_result$V10==0,],
                                     blast_result[blast_result$V5==((W+length_pam_cas)-1) & blast_result$V6==1 & blast_result$V11<=1 & blast_result$V10==0,],
                                     blast_result[blast_result$V5==((W+length_pam_cas)) & blast_result$V6==1 & blast_result$V11<=2 & blast_result$V10==0,])
                }

                else if (CRISPR_Type=="casV_B"|CRISPR_Type=="casV_A"){
                  blast_result=read.table("blast_out.txt",sep = "@")
                  blast_result=rbind(blast_result[blast_result$V5==((W+length_pam_cas)-4) & blast_result$V6==1 & blast_result$V11==0 & blast_result$V10==0,],
                                     blast_result[blast_result$V5==((W+length_pam_cas)-3) & blast_result$V6==1 & blast_result$V11<=1 & blast_result$V10==0,],
                                     blast_result[blast_result$V5==((W+length_pam_cas)-2) & blast_result$V6==1 & blast_result$V11<=2 & blast_result$V10==0,],
                                     blast_result[blast_result$V5==((W+length_pam_cas)-1) & blast_result$V6==1 & blast_result$V11<=3 & blast_result$V10==0,],
                                     blast_result[blast_result$V5==(W+length_pam_cas) & blast_result$V6==1 & blast_result$V11<=4 & blast_result$V10==0,])
                }



                if (length(blast_result$V2)!=0) {

                  for (crRNA_predict in 1:length(data$Guide_RNA)) {
                    list_blast_mm_0[[crRNA_predict]]="No_Hit"
                    list_blast_mm_1[[crRNA_predict]]="No_Hit"
                    list_blast_mm_2[[crRNA_predict]]="No_Hit"

                    if (CRISPR_Type=="casV_B"|CRISPR_Type=="casV_A") {
                      list_blast_mm_3[[crRNA_predict]]="No_Hit"
                      list_blast_mm_4[[crRNA_predict]]="No_Hit"
                    }
                    genes_mm0=c()
                    genes_mm1=c()
                    genes_mm2=c()
                    genes_mm3=c()
                    genes_mm4=c()
                    for (blast_crRNA_predict in 1:length(unique(blast_result$V1))){
                      if (data$Guide_RNA[crRNA_predict]==unique(blast_result$V1)[blast_crRNA_predict]) {
                        s_blast_result=blast_result[blast_result$V1==unique(blast_result$V1)[blast_crRNA_predict],]
                        for (gene in 1:length(s_blast_result$V2))  {
                          s_blast_result_gene=s_blast_result[gene,]
                          pos_mismach_blast_crRNA_gene=which(strsplit(s_blast_result_gene$V3,split = "")[[1]]!=strsplit(s_blast_result_gene$V4,split = "")[[1]])
                          #pam mismatch checked
                          if (any(pos_mismach_blast_crRNA_gene%in%0:length_pam_cas)!=TRUE) {

                            if (CRISPR_Type=="casVI_A"|CRISPR_Type=="casVI_D"|CRISPR_Type=="casVI_B"|CRISPR_Type=="casV_F1"){

                              if ((W+length_pam_cas)-s_blast_result_gene$V5==1) {
                                if (s_blast_result_gene$V6==1) {
                                  pos_mismach_blast_crRNA_gene=c(pos_mismach_blast_crRNA_gene,(W+length_pam_cas))
                                }
                                else if (s_blast_result_gene$V6!=1) {
                                  pos_mismach_blast_crRNA_gene=c(1,pos_mismach_blast_crRNA_gene)
                                }
                              }
                              if ((W+length_pam_cas)-s_blast_result_gene$V5==2) {
                                if (s_blast_result_gene$V6==1) {
                                  pos_mismach_blast_crRNA_gene=c(pos_mismach_blast_crRNA_gene,(W+length_pam_cas)-1,(W+length_pam_cas))
                                }
                                else if (s_blast_result_gene$V6==2) {
                                  pos_mismach_blast_crRNA_gene=c(1,pos_mismach_blast_crRNA_gene,(W+length_pam_cas))
                                }
                                else if (s_blast_result_gene$V6==3) {
                                  pos_mismach_blast_crRNA_gene=c(1,2,pos_mismach_blast_crRNA_gene)
                                }
                              }
                              seq_gene=rep(1,(W+length_pam_cas))
                              seq_gene[pos_mismach_blast_crRNA_gene]=0

                              if (mismatch_del(seq_gene[(length_pam_cas+1):length(seq_gene)],start,end,seed_start,seed_end)!=TRUE) {
                                gene_off_target=paste0(s_blast_result_gene$V2,"_pos",":",s_blast_result_gene$V8,"-",s_blast_result_gene$V9)
                                #mm0
                                if (length(pos_mismach_blast_crRNA_gene)==0) {
                                  genes_mm0=c(genes_mm0,gene_off_target)
                                }
                                #mm1
                                else if (length(pos_mismach_blast_crRNA_gene)==1) {
                                  genes_mm1=c(genes_mm1,gene_off_target)
                                }
                                #mm2
                                else if (length(pos_mismach_blast_crRNA_gene)==2) {
                                  genes_mm2=c(genes_mm2,gene_off_target)
                                }
                              }}

                            if (CRISPR_Type=="casV_B"|CRISPR_Type=="casV_A") {
                              number_mismach_blast_crRNA_gene=(((W+length_pam_cas)-s_blast_result_gene$V5)+s_blast_result_gene$V11)
                              gene_off_target=paste0(s_blast_result_gene$V2,"_pos",":",s_blast_result_gene$V8,"-",s_blast_result_gene$V9)
                              #mm0
                              if (number_mismach_blast_crRNA_gene==0) {
                                genes_mm0=c(genes_mm0,gene_off_target)
                              }
                              #mm1
                              else if (number_mismach_blast_crRNA_gene==1) {
                                genes_mm1=c(number_mismach_blast_crRNA_gene)
                              }

                              #mm2
                              else if (number_mismach_blast_crRNA_gene==2) {
                                genes_mm2=c(genes_mm2,gene_off_target)
                              }

                              #mm3
                              else if (number_mismach_blast_crRNA_gene==3) {
                                genes_mm3=c(genes_mm3,gene_off_target)
                              }

                              #mm4
                              else if (number_mismach_blast_crRNA_gene==4) {
                                genes_mm4=c(genes_mm4,gene_off_target)
                              }
                            }
                          }

                        }


                        #mm0
                        if (length(genes_mm0)>0) {
                          list_blast_mm_0[[crRNA_predict]]=genes_mm0
                          count_mm0[[crRNA_predict]]=count_mm0[[crRNA_predict]]+length(genes_mm0)
                        }

                        #mm1
                        if (length(genes_mm1)>0) {
                          list_blast_mm_1[[crRNA_predict]]=genes_mm1
                          count_mm1[[crRNA_predict]]=count_mm1[[crRNA_predict]]+length(genes_mm1)
                        }

                        #mm2
                        if (length(genes_mm2)>0) {
                          list_blast_mm_2[[crRNA_predict]]=genes_mm2
                          count_mm2[[crRNA_predict]]=count_mm2[[crRNA_predict]]+length(genes_mm2)
                        }

                        #mm3
                        if (length(genes_mm3)>0) {
                          list_blast_mm_3[[crRNA_predict]]=genes_mm3
                          count_mm3[[crRNA_predict]]=count_mm3[[crRNA_predict]]+length(genes_mm3)
                        }

                        #mm4
                        if (length(genes_mm4)>0) {
                          list_blast_mm_4[[crRNA_predict]]=genes_mm4
                          count_mm4[[crRNA_predict]]=count_mm4[[crRNA_predict]]+length(genes_mm4)
                        }

                        break
                      }
                    }

                  }


                  data_off_target$blast_mm_0=list_blast_mm_0
                  data_off_target$blast_mm_1=list_blast_mm_1
                  data_off_target$blast_mm_2=list_blast_mm_2
                  names(data_off_target)[names(data_off_target) == "blast_mm_0"] <- paste("No mismatched off targets in",org)
                  names(data_off_target)[names(data_off_target) == "blast_mm_1"] <- paste("1 mismatched off targets in",org)
                  names(data_off_target)[names(data_off_target) == "blast_mm_2"] <- paste("2 mismatched off Targets in",org)
                  if (CRISPR_Type=="casV_B"|CRISPR_Type=="casV_A"){
                    data_off_target$blast_mm_3=list_blast_mm_3
                    data_off_target$blast_mm_4=list_blast_mm_4
                    names(data_off_target)[names(data_off_target) == "blast_mm_3"] <- paste("3 mismatched off targets in",org)
                    names(data_off_target)[names(data_off_target) == "blast_mm_4"] <- paste("4 mismatched off targets in",org)
                  }


                }

                if (length(blast_result$V2)==0) {
                  data_off_target$blast_mm_0="No_Hit"
                  data_off_target$blast_mm_1="No_Hit"
                  data_off_target$blast_mm_2="No_Hit"
                  names(data_off_target)[names(data_off_target) == "blast_mm_0"] <- paste("No mismatched off targets in",org)
                  names(data_off_target)[names(data_off_target) == "blast_mm_1"] <- paste("1 mismatched off targets in",org)
                  names(data_off_target)[names(data_off_target) == "blast_mm_2"] <- paste("2 mismatched off targets in",org)
                  if (CRISPR_Type=="casV_B"|CRISPR_Type=="casV_A"){
                    data_off_target$blast_mm_3="No_Hit"
                    data_off_target$blast_mm_4="No_Hit"
                    names(data_off_target)[names(data_off_target) == "blast_mm_3"] <- paste("3 mismatched off targets in",org)
                    names(data_off_target)[names(data_off_target) == "blast_mm_4"] <- paste("4 mismatched off targets in",org)
                  }
                }
              }

              if (info_blast$size==0) {
                data_off_target$blast_mm_0="No_Hit"
                data_off_target$blast_mm_1="No_Hit"
                data_off_target$blast_mm_2="No_Hit"
                names(data_off_target)[names(data_off_target) == "blast_mm_0"] <- paste("No mismatched off targets in",org)
                names(data_off_target)[names(data_off_target) == "blast_mm_1"] <- paste("1 mismatched off targets in",org)
                names(data_off_target)[names(data_off_target) == "blast_mm_2"] <- paste("2 mismatched off targets in",org)
                if (CRISPR_Type=="casV_B"|CRISPR_Type=="casV_A"){
                  data_off_target$blast_mm_3="No_Hit"
                  data_off_target$blast_mm_4="No_Hit"
                  names(data_off_target)[names(data_off_target) == "blast_mm_3"] <- paste("3 mismatched off targets in",org)
                  names(data_off_target)[names(data_off_target) == "blast_mm_4"] <- paste("4 mismatched off targets in",org)
                }
              }



            }
            data$No_Mismatch_off_Target=count_mm0
            data$One_Mismatch_off_Targets=count_mm1
            data$Two_Mismatch_off_Targets=count_mm2
            col_number=3

            if (CRISPR_Type=="casV_B"|CRISPR_Type=="casV_A") {
              data$Three_Mismatch_off_Targets=count_mm3
              data$Four_Mismatch_off_Targets=count_mm4
              col_number=col_number+2
            }


            setwd(first_dir)
            html_result_folder="Off_Target-result.html"
            dir.create( html_result_folder)

            number_of_link_for_off=0
            for (mm_type in (length(data)-(col_number-1)):length(data)) {
              for ( guide in 1:length(data$Guide_RNA)) {
                if (data[guide,mm_type]>0) {
                  number_of_link_for_off=number_of_link_for_off+1
                  zz=c("NO_mismatch","1_mismatch","2_mismatches","3_mismatches","4_mismatches")
                  type_of_mm=which(mm_type==(length(data)-(col_number-1)):length(data))
                  number_of_organism=length(Organism)
                  if (number_of_organism==1) {
                    show_col=c(1,2,type_of_mm+2)
                  }
                  if (number_of_organism>1) {
                    show_col_1_plus=(type_of_mm+2)+1:(number_of_organism-1)*col_number
                    show_col=c(1,2,type_of_mm+2,show_col_1_plus)
                  }




                  header.names <- c("",names(data_off_target[guide,][,show_col]))
                  sketch_for_each_gene <-  withTags(table(
                    style(type = "text/css", header.style),
                    thead(
                      tr(
                        lapply(header.names, th, style = "text-align: center ; border-right-width: 1px;border-right-style: solid; border-right-color: white; border-bottom-width: 1px; border-bottom-style: solid; border-bottom-color: white")
                      )
                    )
                  ))


                  rr=as.data.frame(lapply(data_off_target[guide,][,show_col], function(x) {gsub(",", "\n", x)}))
                  cc=as.data.frame(lapply(rr, function(x) {gsub("c\\(", "", x)}))
                  bb=as.data.frame(lapply(cc, function(x) {gsub("\\)", "", x)}))
                  data_each_spacer_off_target=as.data.frame(lapply(bb, function(x) {gsub('"', '', x)}))


                  result_of_guide= datatable(data_each_spacer_off_target, filter = 'none',
                                             options =my.options,
                                             container = sketch_for_each_gene,
                                             rownames = TRUE,
                                             caption = paste0('Results Table for ',CRISPR_Type,' (Guide RNA=',data_off_target[guide,1],')'))

                  result_of_guide <- formatStyle(result_of_guide,
                                                 columns = " ",
                                                 backgroundColor =  "#000000",
                                                 borderBottomColor =  "#000000",
                                                 borderBottomStyle = "solid",
                                                 borderBottomWidth = "1px",
                                                 borderCollapse = "collapse",
                                                 borderRightColor =  "#000000",
                                                 borderRightStyle = "solid",
                                                 borderRightWidth = "1px",
                                                 color = "rgb(255, 255, 255)",
                                                 fontFamily = "Times New Roman",
                                                 fontSize = "13px",
                                                 fontWeight = "bold",
                                                 lineHeight = "normal",
                                                 paddingBottom = "2.6px",
                                                 paddingLeft = "5.2px",
                                                 paddingRight = "5.2px",
                                                 paddingTop = "2.6px",
                                                 textAlign = "center",
                                                 verticalAlign = "middle",

                  )
                  result_of_guide <- formatStyle(result_of_guide,
                                                 columns = 1:(length(Organism)+2),
                                                 fontFamily = "Times New Roman",
                                                 fontSize = "13px",
                                                 lineHeight = "normal",
                                                 paddingBottom = "2.6px",
                                                 paddingLeft = "5.2px",
                                                 paddingRight = "5.2px",
                                                 paddingTop = "2.6px",
                                                 textAlign = "center",
                                                 verticalAlign = "middle")




                  name_html_Guide=paste0(zz[type_of_mm],"_",data[[guide,2]],".html",collapse = "")


                  setwd(html_result_folder)
                  htmlwidgets::saveWidget(result_of_guide, name_html_Guide )

                  setwd(first_dir)
                  Guide_html_Link = tibble( file = dir( file.path( '.', html_result_folder ), pattern= name_html_Guide ) ) %>%
                    mutate( path = paste0('file:', html_result_folder,"/",  sort(file ))
                            , path = stringr::str_replace_all( path, ' ', '%20')
                            , link =  paste0('<a  target=_blank href=', path, '>', data[guide,mm_type],'</a>' ) )


                  data[guide,mm_type]=Guide_html_Link$link
                }

              }

            }


            if (number_of_link_for_off==0 ) {
              if (dir.exists(html_result_folder)) {
                unlink(paste0("./",html_result_folder), recursive = TRUE)
              }

            }



            ############# create sketch for DT package ###########

            zz=c("NO_mismatch","1_mismatch","2_mismatches","3_mismatches","4_mismatches")
            name_subcolumn_off_target=zz[1: col_number]
            header.names <- c("",names(data_off_target))
            # The container parameter allows us to design the header of the table
            # using CSS
            sketch_off_Target <-  withTags(table(
              style(type = "text/css", header.style),
              thead(
                tr(
                  lapply(header.names[1:3], th,rowspan = 2, style = "text-align: center ; border-right-width: 1px;border-right-style: solid; border-right-color: white; border-bottom-width: 1px; border-bottom-style: solid; border-bottom-color: white"),
                  lapply(Organism, th,colspan = col_number, style = "text-align: center ; border-right-width: 1px;border-right-style: solid; border-right-color: white; border-bottom-width: 1px; border-bottom-style: solid; border-bottom-color: white")

                ),
                tr(
                  rep(lapply(name_subcolumn_off_target, th),length(Organism))
                )
              )
            ))

            ######################################################




            rr=as.data.frame(lapply(data_off_target, function(x) {gsub(",", "\n", x)}))
            cc=as.data.frame(lapply(rr, function(x) {gsub("c\\(", "", x)}))
            bb=as.data.frame(lapply(cc, function(x) {gsub("\\)", "", x)}))
            data_off_target=as.data.frame(lapply(bb, function(x) {gsub('"', '', x)}))



            result_off_targets= datatable(data_off_target, filter = 'none', options = my.options,rownames = TRUE,container = sketch_off_Target)

            result_off_targets <- formatStyle(result_off_targets,
                                              columns = " ",
                                              backgroundColor =  "#000000",
                                              borderBottomColor =  "#000000",
                                              borderBottomStyle = "solid",
                                              borderBottomWidth = "1px",
                                              borderCollapse = "collapse",
                                              borderRightColor =  "#000000",
                                              borderRightStyle = "solid",
                                              borderRightWidth = "1px",
                                              color = "rgb(255, 255, 255)",
                                              fontFamily = "Times New Roman",
                                              fontSize = "13px",
                                              fontWeight = "bold",
                                              lineHeight = "normal",
                                              paddingBottom = "2.6px",
                                              paddingLeft = "5.2px",
                                              paddingRight = "5.2px",
                                              paddingTop = "2.6px",
                                              textAlign = "center",
                                              verticalAlign = "middle",

            )
            result_off_targets <- formatStyle(result_off_targets,
                                              columns = 1:length(data_off_target),
                                              fontFamily = "Times New Roman",
                                              fontSize = "13px",
                                              lineHeight = "normal",
                                              paddingBottom = "2.6px",
                                              paddingLeft = "5.2px",
                                              paddingRight = "5.2px",
                                              paddingTop = "2.6px",
                                              textAlign = "center",
                                              verticalAlign = "middle")

            htmlwidgets::saveWidget(result_off_targets, "result2.html")


            html_result=read_file("result2.html")
            line_heder=paste0('<header><h1 style="color:blue;font-size:30px;" align="middle">',
                              paste0('Off_Targets Results Table for ',ID_code,' CRISPR-Cas System'),'</h1>')
            line_img_tehran=paste('<img src=   "https://drive.google.com/thumbnail?id=1NrT3Zw15mdetkgB4xLZ7GlZnDj-uAeNw"','width= 300 height= 130',' align="left" >')
            line_img_MRBLAB=paste('<img src=   "https://drive.google.com/thumbnail?id=1BWDyB9RcBdes4dRNRh5ff3lUxMaUwdME" ','width= 100 height= 100',' align="left" >')


            file_Conaction<-file("Off_Targets_Information.Html")
            writeLines(c(line_heder,line_img_tehran,html_result,line_img_MRBLAB), file_Conaction)

            close(file_Conaction)


            file.remove("result2.html")
            file.remove("blast_out.txt")
            file.remove("blast_query.FASTA")
          }



          if(question == 'n'| question == 'N'){
            message("Running continue without Off_Target prediction" )
          }



          if (havingIP()!=TRUE) {
            message("Problem in Internet connection ",
                    "continue without off_Target result")
          }

        }

        #######################################################################################












        ######################### assign link for figure and sort data ##################
        setwd(first_dir)

        data$Plot=link_plots
        data=data[order(-data$Score),]
        data$Score=NULL
        #################################################################################







        ################### create sketch for DT package ################################
        if (OffTarget==T) {
          if (question== 'y'| question== 'Y') {

            name_subcolumn_off_target=zz[1: col_number]
            header.names <- c("",names(data)[1:(length(names(data))-(col_number+1))])
            # The container parameter allows us to design the header of the table
            # using CSS
            sketch_Specer_informatin <-  withTags(table(
              style(type = "text/css", header.style),
              thead(
                tr(
                  lapply(header.names, th,rowspan = 2, style = "text-align: center ; border-right-width: 1px;border-right-style: solid; border-right-color: white; border-bottom-width: 1px; border-bottom-style: solid; border-bottom-color: white"),
                  lapply("Off_Target", th,colspan = col_number, style = "text-align: center ; border-right-width: 1px;border-right-style: solid; border-right-color: white; border-bottom-width: 1px; border-bottom-style: solid; border-bottom-color: white"),
                  th(rowspan = 2,"Plots",style = "text-align: center ; border-right-width: 1px;border-right-style: solid; border-right-color: white; border-bottom-width: 1px; border-bottom-style: solid; border-bottom-color: white")

                ),
                tr(
                  lapply(name_subcolumn_off_target, th)
                )
              )
            ))
          }

          else if (question== 'n'| question== 'N') {

            header.names <- c("",names(data))
            sketch_Specer_informatin <-  withTags(table(
              style(type = "text/css", header.style),
              thead(
                tr(
                  lapply(header.names, th, style = "text-align: center ; border-right-width: 1px;border-right-style: solid; border-right-color: white; border-bottom-width: 1px; border-bottom-style: solid; border-bottom-color: white")
                )
              )
            ))
          }
        }
        if (OffTarget==F) {

          header.names <- c("",names(data))
          sketch_Specer_informatin <-  withTags(table(
            style(type = "text/css", header.style),
            thead(
              tr(
                lapply(header.names, th, style = "text-align: center ; border-right-width: 1px;border-right-style: solid; border-right-color: white; border-bottom-width: 1px; border-bottom-style: solid; border-bottom-color: white")
              )
            )
          ))
        }

        #################################################################################








        ################# final result html report and save ###############


        result_spacers= datatable(data, filter = 'none', options = my.options,rownames = TRUE,container = sketch_Specer_informatin ,
                                  caption = paste0(direct_repeat,' from ',bacteria,' (',backbone,")"))



        result_spacers <- formatStyle(result_spacers,
                                      columns = " ",
                                      backgroundColor =  "#000000",
                                      borderBottomColor =  "#000000",
                                      borderBottomStyle = "solid",
                                      borderBottomWidth = "1px",
                                      borderCollapse = "collapse",
                                      borderRightColor =  "#000000",
                                      borderRightStyle = "solid",
                                      borderRightWidth = "1px",
                                      color = "rgb(255, 255, 255)",
                                      fontFamily = "Times New Roman",
                                      fontSize = "13px",
                                      fontWeight = "bold",
                                      lineHeight = "normal",
                                      paddingBottom = "2.6px",
                                      paddingLeft = "5.2px",
                                      paddingRight = "5.2px",
                                      paddingTop = "2.6px",
                                      textAlign = "center",
                                      verticalAlign = "middle",

        )
        result_spacers <- formatStyle(result_spacers,
                                      columns = 1:length(data),
                                      fontFamily = "Times New Roman",
                                      fontSize = "13px",
                                      lineHeight = "normal",
                                      paddingBottom = "2.6px",
                                      paddingLeft = "5.2px",
                                      paddingRight = "5.2px",
                                      paddingTop = "2.6px",
                                      textAlign = "center",
                                      verticalAlign = "middle")


        setwd(first_dir)
        htmlwidgets::saveWidget(result_spacers, "result1.html")


        #####################################################################



        #edit HTML result
        html_result=read_file("result1.html")
        line_heder=paste0('<header><h1 style="color:blue;font-size:30px;" align="middle">',
                          paste0('Results Table for ',ID_code,' CRISPR-Cas System'),'</h1>')
        line_img_tehran=paste('<img src=   "https://drive.google.com/thumbnail?id=1NrT3Zw15mdetkgB4xLZ7GlZnDj-uAeNw"','width= 20% height= 20%',' align="left" >')
        line_img_MRBLAB=paste('<img src=    "https://drive.google.com/thumbnail?id=1BWDyB9RcBdes4dRNRh5ff3lUxMaUwdME"','width= 12% height= 20%',' align="left" >')
        line_help='<p>more detailed description of columns is available in the <a href="https://sadegh65v.github.io/casilico/Casilico%20Guide.html" rel="nofollow noreferrer">HELP</a> section</p>'

        file_Conaction<-file("Spacers_Information.Html")
        writeLines(c(line_heder,line_img_tehran,html_result,line_help,line_img_MRBLAB), file_Conaction)

        close(file_Conaction)


        file.remove("result1.html")

        setwd(first_main_dir)
      }


      if (length(crRNA_seq)==0  ) {

        message("No candidate result pridicted for ",CRISPR_Type,
                " please decrease mismatch threshold or change the ConservationMethod used ")
        setwd(first_main_dir)
        unlink(paste0("./",result_folder_name), recursive = TRUE)

      }

    }

  }
  file.remove("final_fasta_file.aln")
  return("Design_complte")
  options(warn = defaultW)

}









#' MAFFT (Multiple Sequence Alignment) in R
#' @param TargetFasta  	Full path to the target sequence(s); Example: TargetFasta="/usr/mrb/casilico/genome.fasta"
#' @param mafft_Options  	Options to run mafft; Example: MafftOptions="--auto   --inputorder"
#' @return msa file (output of mafft), always output format is Clustal and name is final_fasta_file.aln
#' @section Further details: mafftR reports output as a cluster in final_fasta_file.aln file.
#' details for Options at https://mafft.cbrc.jp/alignment/software/manual/manual.html
#' @export mafftR
#'
#'
#'
mafftR=function(TargetFasta, MafftOptions = "--auto   --inputorder"){

  fasta_file=read.fasta(TargetFasta)
   host_system=Sys.info()[['sysname']]

  if (host_system %in% c("Windows","Linux")==F) {
    stop("Currently CaSilico only supports Windows and Linux systems")
  }
  package_location<-system.file(package = "CaSilico")

  if (host_system=="Windows") {
    MAFFT_location=paste0(package_location,"/","dependency_win/mafft-win")
  }


  ############## change mod sotware ###########


  if (host_system=="Linux") {

    MAFFT_location=paste0(package_location,"/dependency_linux/mafft-linux64")
    
    lapply(paste0("chmod -R 777  ",MAFFT_location), system)
     }
  final_fasta_file=fasta_file
  write.fasta(final_fasta_file,names(final_fasta_file),"final_fasta_file.fasta")
  file.copy("final_fasta_file.fasta",MAFFT_location)
  file.remove("final_fasta_file.fasta")


  first_dir=getwd()
  setwd(MAFFT_location)
  MAFFT_list<- paste("mafft    --clustalout",MafftOptions,"  final_fasta_file.fasta > final_fasta_file.aln")
  lapply(MAFFT_list, system)   ###run MAFFT for all input .txt files
  file.copy("final_fasta_file.aln",first_dir)
  file.remove("final_fasta_file.fasta")
  file.remove("final_fasta_file.aln")
  setwd(first_dir)
}






#' Dot plot for short sequence
#' @param sequence1  First sequence
#' @param sequence2  Second sequence
#' @param wsize the  Size in chars of the moving window.
#' @param wstep the  Size in chars for the steps of the moving window. Use wstep == wsize for non-overlapping windows.
#' @param nmatch     If the number of match per window is greater than or equal to nmatch, then a dot is produced.
#' @param table_col  Table color by default table_col = "goldenrod1".
#' @param points_col Point color by default points_col = "red".
#' @param xlab       Label for x axis.
#' @param ylab       Label for y axis.
#' @return plots.
#' @section References:
#' ...
#' @section Further details:In order to visualize this function better, short sequences are recommended.
#' ...
#' @export VDotplot
#'
#'

VDotplot=function(sequence1,sequence2, wsize = 1 , wstep = 1 , nmatch = 1 ,
                              table_col = "goldenrod1",points_col = "red", xlab = "sequence1", ylab = "sequence2")
{
  #Split the Elements of strin
  if (length(sequence1)==1){
    x <- strsplit(sequence1,split = "")[[1]]}
  if (length(sequence2)==1){
    y <- strsplit(sequence2,split = "")[[1]]}
  if (length(sequence1)!=1){
    x <-sequence1}
  if (length(sequence2)!=1){
    y <- sequence2}





  if (nchar(x[1]) > 1)
    stop("seq1 should be provided as a vector of single chars")
  if (nchar(y[1]) > 1)
    stop("seq2 should be provided as a vector of single chars")
  if (wsize < 1)
    stop("non allowed value for wsize")
  if (wstep < 1)
    stop("non allowed value for wstep")
  if (nmatch < 1)
    stop("non allowed value for nmatch")
  if (nmatch > wsize)
    stop("nmatch > wsize is not allowed")
  mkwin <- function(seq, wsize, wstep) {
    sapply(seq(from = 1, to = length(seq) - wsize + 1, by = wstep),
           function(i) c2s(seq[i:(i + wsize - 1)]))
  }
  wseq1 <- mkwin(x, wsize, wstep)
  wseq2 <- mkwin(y, wsize, wstep)
  if (nmatch == wsize) {
    xy <- outer(wseq1, wseq2, "==")
  }
  else {
    "%==%" <- function(x, y) colSums(sapply(x, s2c) == sapply(y,
                                                              s2c)) >= nmatch
    xy <- outer(wseq1, wseq2, "%==%")
  }

  x=toupper(x)
  y=toupper(y)
  a=matrix(0,length(x),length(y))
  # Naming the rows and columns of the formed matrix
  row.names(a)=x
  colnames(a)=y


  for (i in 1:nrow(xy)) {
    for (j in 1:ncol(xy)){
      if (xy[i,j]==T) {
        a[i,j]=1
      }
    }
  }



  h_line=seq(0+((seq(0,max(length(x),length(y))-1,
                     l=length(x))[2]-seq(0,max(length(x),length(y))-1,
                                         l=length(x))[1])/2),max(length(x),length(y)),
             by=seq(0,max(length(x),length(y))-1,
                    l=length(x))[2]-seq(0,max(length(x),length(y))-1,
                                        l=length(x))[1])
  v_line=seq(1+((seq(1,max(length(x),length(y)),
                     l=length(y))[2]-seq(1,max(length(x),length(y)),
                                         l=length(y))[1])/2),max(length(x),length(y)),
             by=seq(1,max(length(x),length(y)),
                    l=length(y))[2]-seq(1,max(length(x),length(y)),
                                        l=length(y))[1])

  aaaa=0.04*max(length(x),length(y))


  plot(0:max(length(x),length(y)),0:max(length(x),length(y)),type = "n",
       axes =F,frame.plot=TRUE,main ="" ,
       ann = T,xlab = xlab , ylab = ylab)



  aa=c(-aaaa,-aaaa,max(length(x),length(y))+aaaa,max(length(x),length(y))+aaaa)
  bb=c(h_line[length(h_line)],
       max(length(x),length(y))+aaaa,
       max(length(x),length(y))+aaaa,
       h_line[length(h_line)])

  polygon(aa,bb,col = table_col)


  aa1=c(-aaaa,-aaaa,max(length(x),length(y))+aaaa,max(length(x),length(y))+aaaa)
  bb1=c(max(length(x),length(y))-aaaa,
        max(length(x),length(y))+0.5,
        max(length(x),length(y))+0.5,
        max(length(x),length(y))-aaaa)

  polygon(bb1-max(length(x),length(y)),aa1,col = table_col)



  # y axes
  axis(side=2, at = seq(0,max(length(x),length(y))-1,l=length(x)),
       labels= length(x):1,cex.axis = 0.7, font = 2)

  abline(h=h_line,col="black",lwd=1, lty=1)



  text((aa1[2]+h_line[1])/2,
       seq(0,max(length(x),length(y))-1,l=length(x)),
       rev(x),col = "black",cex = 1,pch = 16,font = 2)


  #x axes
  axis(side=3, at = seq(1,max(length(x),length(y)),
                        l=length(y)), labels= 1:length(y),cex.axis = 0.7, font = 2)


  abline(v=v_line, col="black",lwd=1, lty=1)


  abline(v=0.5,col="black",lwd=1, lty=1)

  text(seq(1,max(length(x),length(y)),l=length(y))
       ,(bb[2]+bb[1])/2,
       y,col = "black",cex = 1,pch = 16,font = 2)



  for (i in 1:length(x)){
    for (j in 1:length(y))
      if (a[i,j]==1)
        points(seq(1,max(length(x),length(y)),l=length(y))[j],seq(0,max(length(x),length(y))-1,l=length(x))[length(x)-i+1],cex=1.5,col=points_col,pch = 16)
  }
}





#' Compute and return revers compliment of a sequence
#' @param seq  Sequence of DNA or RNA.
#' @return Reverse compliment of sequence .
#' @section References:
#' ...
#' @section Further details:
#' ...
#' @export Rcompliment
#'
#'
Rcompliment=function (seq)
{
  compliment = function(n, seq_type) {
    if (seq_type == "D") {
      if (n == "a") {
        return("t")
      }
      if (n == "A") {
        return("T")
      }
    }
    if (seq_type == "R") {
      if (n == "a") {
        return("u")
      }
      if (n == "A") {
        return("U")
      }
    }
    if (n == "t") {
      return("a")
    }
    if (n == "c") {
      return("g")
    }
    if (n == "g") {
      return("c")
    }
    if (n == "C") {
      return("G")
    }
    if (n == "G") {
      return("C")
    }
    if (n == "T") {
      return("A")
    }
    if (n == "u") {
      return("a")
    }
    if (n == "U") {
      return("A")
    }
    if (n != "U" & n != "A" & n != "T" &
        n != "G" & n != "C" & n != "u" &
        n != "a" & n != "t" & n != "g" &
        n != "c") {
      return(n)
    }
  }
  Rcompliment_sequence = c()
  if (length(seq) > 1) {
    seq = seq
  }
  if (length(seq) == 1) {
    seq = strsplit(seq, split = "")[[1]]
  }
  type = "D"
  if (any(seq == "u") | any(seq == "U")) {
    type = "R"
  }

  for (i in length(seq):1) {
    Rcompliment_sequence = c(Rcompliment_sequence, compliment(seq[i],
                                                              seq_type = type))
  }
  Rcompliment_sequence
}









#' Compute and return compliment of a sequence.
#' @param seq  Sequence of DNA or RNA
#' @return Compliment of sequence .
#' @section References:
#' ...
#' @section Further details:
#' ...
#' @export compliment
#'
#'
compliment=function (seq)
  {
    compliment = function(n, seq_type) {
      if (seq_type == "D") {
        if (n == "a") {
          return("t")
        }
        if (n == "A") {
          return("T")
        }
      }
      if (seq_type == "R") {
        if (n == "a") {
          return("u")
        }
        if (n == "A") {
          return("U")
        }
      }
      if (n == "t") {
        return("a")
      }
      if (n == "c") {
        return("g")
      }
      if (n == "g") {
        return("c")
      }
      if (n == "C") {
        return("G")
      }
      if (n == "G") {
        return("C")
      }
      if (n == "T") {
        return("A")
      }
      if (n == "u") {
        return("a")
      }
      if (n == "U") {
        return("A")
      }
      if (n != "U" & n != "A" & n != "T" &
          n != "G" & n != "C" & n != "u" &
          n != "a" & n != "t" & n != "g" &
          n != "c") {
        return(n)
      }
    }
    compliment_sequence = c()
    if (length(seq) > 1) {
      seq = seq
    }
    if (length(seq) == 1) {
      seq = strsplit(seq, split = "")[[1]]
    }
    type = "D"
    if (any(seq == "u") | any(seq == "U")) {
      type = "R"
    }

    for (i in 1:length(seq)) {
      compliment_sequence = c(compliment_sequence, compliment(seq[i],
                                                              seq_type = type))
    }
    compliment_sequence
  }

