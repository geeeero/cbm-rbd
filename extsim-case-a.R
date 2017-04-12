# --------------------------------------------------------------------------- #
# CBM systems paper - extended simulation study for revision -
# Case A: failure times as expected
# --------------------------------------------------------------------------- #

# packages, code files, example system
source("extsim-setup.R")

# simulate 100 25-cycle machines with priors as given
set.seed(2211)
repetitions <- 100
opcycles <- 25
Asimdata <- list()
for (i in 1:repetitions){
  Asimdata[[i]] <- brWeibullData(opcycles, br1beta, br1mttf)
}

# result lists
AsimCBMcpu <- list() # our model
AsimCBMepu <- list() # do not update params during cycle, but at end of cycle 
AsimCBMnpu <- list() # never update params
AsimABMepu <- list() # age-based policy with parameter update
AsimABMnpu <- list() # age-based policy without parameter update
AsimCM <- list()     # corrective policy

# run simulations (six separate loops possible - use plyr::llply or lapply instead?)
for (i in 1:repetitions){
  #cat("Repetition", i, ": CBM-cpu\n")
  AsimCBMcpu[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = Asimdata[[i]], n0y0 = br1n0y0,
                               beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401,
                               onecycle = FALSE)
  #cat("Repetition", i, ": CBM-epu\n")
  AsimCBMepu[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = Asimdata[[i]], n0y0 = br1n0y0,
                               beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401,
                               prior = TRUE, onecycle = FALSE)
  #cat("Repetition", i, ": CBM-npu\n")
  AsimCBMnpu[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = Asimdata[[i]], n0y0 = br1n0y0,
                               beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401,
                               prior = TRUE, cycleupdate = FALSE, onecycle = FALSE)
  #cat("Repetition", i, ": ABM-epu\n")
  AsimABMepu[[i]] <- simNcycleAgebased(sys = br, ctypes = brctypes, compfts = Asimdata[[i]], n0y0 = br1n0y0,
                                       beta = br1beta, tnowstep = 0.1, hor = 4, seqlen = 401)
  #cat("Repetition", i, ": ABM-npu\n")
  AsimABMnpu[[i]] <- simNcycleAgebased(sys = br, ctypes = brctypes, compfts = Asimdata[[i]], n0y0 = br1n0y0,
                                       beta = br1beta, tnowstep = 0.1, hor = 4, seqlen = 401, cycleupdate = FALSE)
  #cat("Repetition", i, ": CM\n")
  AsimCM[[i]] <- simNcycleCorrective(sys = br, ctypes = brctypes, compfts = Asimdata[[i]], tnowstep = 0.1)
}

# save result lists
str1 <- "Asim"
str2 <- c("CBMcpu", "CBMepu", "CBMnpu", "ABMepu", "ABMnpu", "CM")

ststr <- paste("br1sim", 1:3, sep="")
enstr <- c("Tt01", "prTt01", "prprTt01", "a", "apr", "c")
savelist <- as.vector(sapply(ststr, function(brstr) paste(brstr, enstr, sep="")))
save(list = savelist, file = "papersimobjects")
load("papersimobjects")

#