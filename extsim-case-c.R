# --------------------------------------------------------------------------- #
# CBM systems paper - extended simulation study for revision -
# Case C: failure times later than expected
# --------------------------------------------------------------------------- #

# packages, code files, example system
source("extsim-setup.R")

# simulate 100 25-cycle machines
set.seed(2411)
repetitions <- 100
opcycles <- 25
Csimdata <- list()
for (i in 1:repetitions){
  Csimdata[[i]] <- brWeibullData(opcycles, br1beta, 2*br1mttf)
}

# result lists
CsimCBMcpu <- list() # our model
CsimCBMepu <- list() # do not update params during cycle, but at end of cycle 
CsimCBMnpu <- list() # never update params
CsimABMepu <- list() # age-based policy with parameter update
CsimABMnpu <- list() # age-based policy without parameter update
CsimCM <- list()     # corrective policy

# run simulations (six separate loops possible - use plyr::llply or lapply instead?)
for (i in 1:repetitions){
  #cat("Repetition", i, ": CBM-cpu\n")
  CsimCBMcpu[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = Csimdata[[i]], n0y0 = br1n0y0,
                               beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401,
                               onecycle = FALSE)
  #cat("Repetition", i, ": CBM-epu\n")
  CsimCBMepu[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = Csimdata[[i]], n0y0 = br1n0y0,
                               beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401,
                               prior = TRUE, onecycle = FALSE)
  #cat("Repetition", i, ": CBM-npu\n")
  CsimCBMnpu[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = Csimdata[[i]], n0y0 = br1n0y0,
                               beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401,
                               prior = TRUE, cycleupdate = FALSE, onecycle = FALSE)
  #cat("Repetition", i, ": ABM-epu\n")
  CsimABMepu[[i]] <- simNcycleAgebased(sys = br, ctypes = brctypes, compfts = Csimdata[[i]], n0y0 = br1n0y0,
                                       beta = br1beta, tnowstep = 0.1, hor = 4, seqlen = 401)
  #cat("Repetition", i, ": ABM-npu\n")
  CsimABMnpu[[i]] <- simNcycleAgebased(sys = br, ctypes = brctypes, compfts = Csimdata[[i]], n0y0 = br1n0y0,
                                       beta = br1beta, tnowstep = 0.1, hor = 4, seqlen = 401, cycleupdate = FALSE)
  #cat("Repetition", i, ": CM\n")
  CsimCM[[i]] <- simNcycleCorrective(sys = br, ctypes = brctypes, compfts = Csimdata[[i]], tnowstep = 0.1)
}

# save result lists
str1 <- "Csim"
str2 <- c("CBMcpu", "CBMepu", "CBMnpu", "ABMepu", "ABMnpu", "CM")
save(list = paste(str1, str2), file = "extsim-c-objects.RData", compression = TRUE)


#