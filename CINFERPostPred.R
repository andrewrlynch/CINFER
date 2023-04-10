# function to run posterior predictive NetLogo simulation in CINFER
CINFERPostPred <- function(r, s, t){
   Chrs <- paste("Chr", c(1:22, "X"), c("p", "q"), sep = "")
   trim.brackets <- function(x) gsub("\\[|\\]", "", x)

   netlogopath <- file.path("NetLogo 6.3.0")
   # modelpath <- file.path(netlogopath, ifelse(g == "exponential",
   #                                            "models/2021-06-28-ARL_segmental_aneuploidy.nlogo",
   #                                            ifelse(g == "constant wright-fisher", 
   #                                                   "models/2021-07-02-ARL_constant_wrightfisher_segmental_aneuploidy.nlogo", 
   #                                                   stop())))
   modelpath <- file.path(netlogopath, "models/2021-06-28-ARL_segmental_aneuploidy.nlogo")
   outpath <- file.path(netlogopath)
   
   nl <- nl(nlversion = "6.3.0",
            nlpath = netlogopath,
            modelpath = modelpath,
            jvmmem = 1024)
   
   nl@experiment <- experiment(expname="postpred",
                               outpath=outpath,
                               repetition=5,
                               tickmetrics="false",
                               idsetup="setup",
                               idgo="go",
                               runtime=round(t),
                               evalticks=round(t),
                               metrics=c("count tumorcells",
                                         Chrs),
                               variables = list("s" = list(values = c(0, s))),
                               constants = list("initial-number" = 100,
                                                "end-population" = 100000,
                                                "pmisseg" = round(r / 46, 3),
                                                "initial-ploidy" = 2,
                                                "end-ticks" = round(t),
                                                "outputcount" = 100,
                                                "periodicity" = "\"constant\"",
                                                "mitosisProbability" = 2,
                                                "pbreak" = 0,
                                                "model" = "\"Abundance\""
                                                #"model" = sprintf("\"%s\"", m)
                                                )
                               )

            nl@simdesign <- simdesign_ff(nl=nl,
                                       nseeds=1)
      
      # Evaluate nl object:
      eval_variables_constants(nl)
      
      # Run all simulations (loop over all siminputrows and simseeds)
      chrdata <- withProgressShiny(run_nl_all(nl), message = "Running simulations. This may take a bit.")
      
      # Data post-process
      chrdata <- data.frame(chrdata[,c("[run number]", "pmisseg", "s", "[step]")], apply(chrdata[,Chrs], 2, trim.brackets))
      chrdata <- cSplit(chrdata, Chrs, sep = " ", direction = "long")
      postpredsumstats <- ddply(chrdata, c("X.run.number.", "s"), function(x){
         chrMat <- as.matrix(as.data.frame(x)[,Chrs])
         if(nrow(na.omit(chrMat)) < 3){
            return(data.frame(
               #AvgPloidy = NA,
               AvgAneuploidy = NA,
               #MKV = NA,
               NormMKV = NA,
               #AvgEuc = NA,
               #NormAvgEuc = NA,
               #Sackin = NA,
               #Colless = NA,
               #SackinNorm = NA ,
               CollessNorm = NA ,
               #Cherries = NA,
               NormCherries = NA,
               #Pitchforks = NA,
               #NormPitchforks = NA,
               #AvgLadder = NA,
               #NormAvgLadder = NA,
               #Stairs1 = NA,
               #Stairs2 = NA
            ))
         }
         
         clusMat <- hclust(dist(chrMat))
         
         colless.permute <- vector()
         for(i in 1:100){
            chrMat.permute <- apply(chrMat, 2, permute)
            chrClust.permute <- hclust(dist(chrMat.permute))
            colless.permute[i] <- colless.phylo(as.phylo(chrClust.permute), normalise = T)
         }
         
         data.frame(
            AvgAneuploidy = mean(apply(chrMat, 1, var, na.rm = T)), 
            NormMKV = mean(apply(chrMat, 2, var)) / mean(apply(chrMat, 1, mean)), 
            NormUniqueClones = nrow(unique(chrMat)) / nrow(chrMat),
            SackinNorm = sackin.phylo(as.phylo(clusMat), normalise = TRUE), 
            CollessNorm = colless.phylo(as.phylo(clusMat), normalise = TRUE),
            CollessPermute = mean(colless.permute),
            NormCherries = cherries(as.phylo(clusMat)) / nrow(chrMat) 
         )
      })
      postpredsumstats <- postpredsumstats %>% 
         pivot_longer(cols = -c(X.run.number., s))
      return(postpredsumstats)
}
