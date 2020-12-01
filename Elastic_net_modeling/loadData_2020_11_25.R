makeDataTables <- function(filepaths){
  
  # Load data table CLIN and break it into observations and labs
  # Universal ids are in rownames()
  CLIN <- read.table(filepaths["clin"], sep="\t", row.names = 1, header=TRUE,
                     na.strings =c("", "  ", "NA"))
  if(any(CLIN$Diet.Protein == "Red  Meat" | CLIN$Diet.Protein == "White meat")){
    CLIN$Diet.Protein[CLIN$Diet.Protein == "Red  Meat"] <- "Red Meat"
    CLIN$Diet.Protein[CLIN$Diet.Protein == "White meat"] <- "White Meat"
  }
  CLIN <- CLIN[,sapply(1:ncol(CLIN), function(x) !(all(is.na(CLIN[,x]))))]
  num_fields <- c(
    "Age.At.Diet.Start",
    "Diet.Start.Date",
    "TG.Ave.",
    "TC.Ave.",
    "LDL.Ave.",
    "HDLC.Ave.",
    "ApoA1.Ave.",
    "ApoB.Ave.",
    "BMI",
    "Per_BF",
    "SBP.Avg.",
    "DBP.Avg.",
    "Endopat_AI",
    "Endopat_AI_75Per",
    "Hip..cm.",
    "Waist..cm.",
    "WaistIC..cm.",
    "Wt...Lbs..Avg.",
    "Height..cm.",
    "Diet..Compliance"
  )
  
  # Convert dates to days since 1 Jan 1970
  require(lubridate)
  CLIN$Diet.Start.Date <- as.numeric(mdy(CLIN[,"Diet.Start.Date"]))

  CLIN_num <- CLIN[,num_fields]
  
  # Replace factors with dummy variables
  categ <- names(CLIN)[!(names(CLIN) %in% num_fields)]
  ref_grp <- c(
    "Sat..Fat.Arm-Low",
    "Sex-Female",
    "Diet.Order-NWR",
    "X1st.Diet-N",
    "X2nd.Diet-N",
    "X3rd.Diet-N",
    "Diet.Protein-Baseline"
  )
  
  # Format Dummy variables as matrices and cbind to CLIN
  dummies <- sapply(1:length(categ), function(x){
    my.categ <- CLIN[,categ[x]]
    my.categ <- droplevels(factor(my.categ))
    nu <- as.numeric(my.categ)
    l  <- unique(nu)
    m  <- sapply(l, function(y) as.integer(nu==y))
    m  <- m[,apply(m, 2, function(y) !all(is.na(y)))]
    colnames(m) <- paste(categ[x],levels(my.categ),sep="-")
    mm <- m[,-which(colnames(m)==ref_grp[x]),drop=FALSE]
    (mm)
  })
  dum <- cbind(dummies[[1]],dummies[[2]])
  for( i in 3:length(dummies)){ dum <- cbind(dum,dummies[[i]])}
  
  CLIN <- cbind(CLIN_num, dum)
  
  
  #### Loading raw count tables
  
  # Load count table, change labels to match TMAL
  cutC = loadCtTab(filepaths, gene="cutC", rownames(CLIN), t=FALSE)
  # yea = loadCtTab(filepaths, gene="yea", rownames(CLIN), t=FALSE)
  # cai = loadCtTab(filepaths, gene="cai", rownames(CLIN), t=FALSE)
  # cut = loadCtTab(filepaths, gene="cut", rownames(CLIN), t=FALSE)
  # cnt = loadCtTab(filepaths["cnt"], filepaths["lib"], gene="cnt", t=FALSE),
  # tor = loadCtTab(filepaths, gene="tor", rownames(CLIN), t=FALSE)
  # grd = loadCtTab(filepaths, gene="grd", rownames(CLIN), t=FALSE)

  
  #### Loading Peter's USEARCH Alignment Counts
  
  # Peter.cts.pp <- preprocessPeter(filepaths, rownames(CLIN))
  
  # Fecal mass spec data
  FLABS <- loadFecalData(filepaths["fdata"])
  
  # Serum mass spec data
  SLABS <- loadSerumData(filepaths["smap"], filepaths["pdata"])
  
  SLABS <- SLABS[,c("Plasma.TMAO...然.",
                    "Choline...然.",
                    "Betaine..uM.",
                    "Carnitine..然.",
                    "Butyrobetaine..然.",
                    "Crotonobetaine..然."
  )]
  
  # SLABS$TMAO.high <- ifelse(SLABS$Plasma.TMAO...然. >= 4.5, 1, 0)
  
  # load("2020-02-18/fer.tmao.Rdata")
  load("urine.tmao.Rdata")
  
  
  my.data <- list(Base=CLIN,FLABS=FLABS,SLABS=SLABS,
                  cutC=cutC,
                  Urinetmao = urine.tmao)
 
  (
    list(dat = mergeByRows(my.data),
         cols =  rep(names(my.data), sapply(my.data,ncol)))
  )
  
}

makeDeltaTable <- function(clinPath){
  
  clin <- read.table(clinPath, header=TRUE, sep="\t")
  delta.table <- data.frame(PID = character(), Non_Meat=character(), 
                            Red_Meat=character())
  for(i in unique(clin$PID)){
    r = clin[which(clin$Diet == "Red_Meat" & clin$PID == i),"Sample"]
    if(length(r) == 0) r <- NA
    n = clin[which(clin$Diet == "Non_Meat" & clin$PID == i),"Sample"]
    if(length(n) == 0) n <- NA
    
    delta.table <- rbind(delta.table, data.frame(PID = i, 
                                                 Non_Meat=n, 
                                                 Red_Meat=r
    )
    )
  }
  
  delta.table$Non_Meat <- sapply(delta.table$Non_Meat, function(x){
    substr(x, 0, nchar(as.character(x))-1)
  })
  
  delta.table$Red_Meat <- sapply(delta.table$Red_Meat, function(x){
    substr(x, 0, nchar(as.character(x))-1)
  })
  
  delta.table <- delta.table[complete.cases(delta.table),]
  
  (delta.table)
  
}

mkDataTableDelta <- function(clinPath, dat, datCols){
  ## returns dtab: absolute change, ddtab: percent change, 
  ## predictor set labels
  
  ptab <- makeDeltaTable(clinPath)
  
  baselineCols <- c(
    "Age.At.Diet.Start",
    "Diet.Start.Date",
    "Sex-Male",
    "Diet.Order-NRW",
    "Diet.Order-RNW",
    "Diet.Order-RWN",
    "Diet.Order-WNR",
    "Diet.Order-WRN",
    "X1st.Diet-R",
    "X1st.Diet-W",
    "X2nd.Diet-R",
    "X2nd.Diet-W",
    "X3rd.Diet-R",
    "X3rd.Diet-W"
  )
  
  excludeCols <- c(
    "Diet.Protein-Non-Meat",
    "Diet.Protein-Red Meat",
    "Diet.Protein-White Meat"
  )
  
  datBaseline <- dat[,baselineCols]
  datT        <- dat[,-grep(paste(c(baselineCols,excludeCols),collapse="|"), names(dat))]
  
  dtab <- matrix(nrow=nrow(ptab), ncol=ncol(datT),
                 dimnames=list(ptab$PID, colnames(datT)))
  
  ddtab <- matrix(nrow=nrow(ptab), ncol=ncol(datT),
                  dimnames=list(ptab$PID, colnames(datT)))
  
  my.min <- min(dat[dat>0], na.rm=TRUE) 
  
  for(i in ptab$PID){
    r = unlist(datT[ptab[ptab$PID == i, "Red_Meat"],])
    n = unlist(datT[ptab[ptab$PID == i, "Non_Meat"],])
    n = sapply(n, function(x) ifelse(x != 0, x, my.min ))
    dtab[i,]  <- r - n
    ddtab[i,] <- 100*(r - n) / n
  }
  
  dtab  <- cbind(dtab,  datBaseline[ptab$Non_Meat,])
  ddtab <- cbind(ddtab, datBaseline[ptab$Non_Meat,])
  
  (
    list(
      dtab  = dtab,
      ddtab = ddtab,
      cols  = unlist(
        sapply(colnames(dtab), function(x){  # pred categories
          datCols[grep(paste('^',x,'$',sep=''), colnames(dat))]    # mapped from datCols
        })
      )
    )
  )
  
  
}

#######################################
#
#  makeDataTables Supporting Methods
#    1. loadDiets      - Formatting diet data to construct serum labs IDs
#    2. loadCtTab      - Preprocess and format read count tables
#    3. loadFecalData  - Formatting fecal mass spectrometry data
#    4. loadSerumData  - Formatting serum mass spectrometry data
#
#######################################

loadDiets <- function(f.path){
  METADATA.APP <- read.table(f.path, 
                             sep="\t", header=T, stringsAsFactors = FALSE)
  md <- METADATA.APP[,c("Label", "PID", "Diet.Protein", "Diet..Compliance")]
  
  # recode ambiguous variables
  md$Diet.Protein[md$Diet.Protein == "Red  Meat"] <- "Red Meat"
  md$Diet.Protein[md$Diet.Protein == "White meat"] <- "White Meat"
  
  ( md )
}

loadCtTab <- function(fp, gene,
                      t=FALSE,
                      test=FALSE){
  ## fp is the full path to the raw count table
  
  if(gene == "cutC"){
    ctTab <- read.table(fp[gene],
                        sep = "\t", header=T, row.names=1)
    Len = 2.545
    
  }
  else{
  ctTab <- read.table(grep(paste(gene,"_bgc",sep=""),fp,value=TRUE),
                      sep = "\t", header=T, row.names=1)
  
  # Homolog lengths
  Len <- sapply(colnames(ctTab), function(x){
    spl <- strsplit(x, "\\.")[[1]]; l <- length(spl)
    (as.numeric(spl[l]) - as.numeric(spl[l-1]))
  })
  
  Len <- Len / 1e3  # convert to kB
  
  }
  
  if(t) ctTab <- t(ctTab)
  
  
  # Adjust the rownames
  rownames(ctTab) <- sapply(rownames(ctTab), function(x) substring(x, 1, nchar(x)-1))
  
  ## Filter columns from trailing whitespace
  ctTab <- ctTab[,colSums(ctTab) > 0 & !(apply(ctTab, 2, function(x) all(is.na(x)))), 
                 drop = FALSE]
  
  ## Normalize to seqDepth and average homolog length
  # Load seq depth
  seqDepth <- read.table(fp["lib"],sep="\t", row.names=1)
  seqDepth <- seqDepth[rownames(ctTab),]
  seqDepth = seqDepth / 1e6 # convert to millions of reads
  
  # Divide by library size
  if(nrow(ctTab) != length(seqDepth)) stop("Error")
  ctTab.lib <- apply(ctTab, 2, function(x) x / seqDepth)
  
  # Divide by refence length
  if(ncol(ctTab) != length(Len)) stop("Error")
  ctTab.lib.len <- sweep(ctTab.lib, 1, Len, FUN = "/")
  # ctTab.lib.len <- t(apply(ctTab.lib, 1, function(x) x / Len))
  
  ## Filtering done in pp
  
  (ctTab.lib.len)
  
  
}

loadFecalData <- function(f.path){
  
  load(f.path)
  
  TMAL <- FECAL.APP
  
  TMAL <- TMAL[TMAL[,"Fecal.PMC.ID"] != '' & !is.na(TMAL[,"Fecal.PMC.ID"]),]  # Remove measurements with no sample label
  
  rownames(TMAL) <- sapply(TMAL[,"Fecal.PMC.ID"], function(x) substr(x, 3, nchar(as.character(x)))) # Remove 'AP' from sample label
  
  # Sort by sample IDs
  TMAL <- TMAL[order(rownames(TMAL)),c("d6.TMA..18.Hr..然.",
                                       "d9.TMA..18.Hr..然.",
                                       "TMA..18.Hr..然.",
                                       "d6.Choline..18.Hr..然.",
                                       "d9.Carnitine..18.Hr..然.",
                                       "d9.GBB..18.Hr..然.",
                                       "d9.CTB.18.Hr.Area.Ratio",
                                       "d6.TMA..36.Hr..然.",
                                       "d9.TMA..36.Hr..然.",
                                       "TMA..36.Hr..然.",
                                       "d6.Choline..36.Hr..然.",
                                       "d9.Carnitine..36.Hr..然.",
                                       "d9.GBB..36.Hr..然.",
                                       "d9.CTB.36.Hr.Area.Ratio")]       
  (TMAL)
}

loadSerumData <- function(Smap,f.path){
  
  ## Smap has two columns: Krauss ID and FecalID
  
  load(f.path)
  Smap <- read.table(Smap, sep="\t")
  
  PLASMA.APP <- PLASMA.APP[PLASMA.APP$Krauss.Plasma.ID %in% Smap[,2],]
  rownames(PLASMA.APP) <- sapply(PLASMA.APP$Krauss.Plasma.ID, function(x) Smap[Smap[,2] == x, 1][1])
  
  # Sort by sample IDs
  PLASMA.APP <- PLASMA.APP[order(rownames(PLASMA.APP)),]       
  (PLASMA.APP)
  
}


preprocessPeter <- function(fp, samps){
  Peter.cts <- read.csv(fp["Pcts"],  header=TRUE, row.names=1)
  rownames(Peter.cts) <- sapply(rownames(Peter.cts), function(x)
    substring(x, 1, nchar(x)-1))
  
  Peter.cts <- Peter.cts[samps,]
  Peter.cts <- Peter.cts[apply(Peter.cts,1,function(x) !all(is.na(x))),]
  
  
  Peter.lens   <- read.csv(fp["Plens"], header=TRUE, row.names=1)
  
  Peter.lens$Length <- Peter.lens$Length / 1e3  # Convert to kb
  
  Peter.lens <- Peter.lens[colnames(Peter.cts),]
  
  seqDepth <- read.table(fp["Plib"],sep=",",header=T, row.names=1)
  seqDepth <- seqDepth / 1e6 # Convert to Mb
  rownames(seqDepth) <- sapply(rownames(seqDepth), function(x)
    substring(x, 1, nchar(x)-1))
  seqDepth <- seqDepth[rownames(Peter.cts),]
  
  # Divide by library size
  if(nrow(Peter.cts) != length(seqDepth)) stop("Error 1")
  Peter.cts.lib <- apply(Peter.cts, 2, function(x) x / seqDepth)
  
  # Divide by refence length
  if(ncol(Peter.cts) != length(Peter.lens)) stop("Error 2")
  Peter.cts.lib.len <- t(apply(Peter.cts.lib, 1, function(x) x / Peter.lens))
  
  (Peter.cts.lib.len)
  
}




