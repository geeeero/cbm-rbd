# --------------------------------------------------------------------------- #
# CBM systems paper - extended simulation study for revision -
# Case B: failure times earlier than expected
# --------------------------------------------------------------------------- #

# packages, code files, example system
source("extsim-setup.R")

# simulate 100 25-cycle machines
set.seed(2211)
repetitions <- 100
opcycles <- 25
Bsimdata <- list()
for (i in 1:repetitions){
  Bsimdata[[i]] <- brWeibullData(opcycles, br1beta, 0.5*br1mttf)
}

# result lists
BsimCBMcpu <- list() # our model
BsimCBMepu <- list() # do not update params during cycle, but at end of cycle 
BsimCBMnpu <- list() # never update params
BsimABMepu <- list() # age-based policy with parameter update
BsimABMnpu <- list() # age-based policy without parameter update
BsimCM <- list()     # corrective policy

# run simulations (six separate loops possible - use plyr::llply or lapply instead?)
for (i in 1:repetitions){
  #cat("Repetition", i, ": CBM-cpu\n")
  BsimCBMcpu[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = Bsimdata[[i]], n0y0 = br1n0y0,
                               beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401,
                               onecycle = FALSE)
  #cat("Repetition", i, ": CBM-epu\n")
  BsimCBMepu[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = Bsimdata[[i]], n0y0 = br1n0y0,
                               beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401,
                               prior = TRUE, onecycle = FALSE)
  #cat("Repetition", i, ": CBM-npu\n")
  BsimCBMnpu[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = Bsimdata[[i]], n0y0 = br1n0y0,
                               beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401,
                               prior = TRUE, cycleupdate = FALSE, onecycle = FALSE)
  #cat("Repetition", i, ": ABM-epu\n")
  BsimABMepu[[i]] <- simNcycleAgebased(sys = br, ctypes = brctypes, compfts = Bsimdata[[i]], n0y0 = br1n0y0,
                                       beta = br1beta, tnowstep = 0.1, hor = 4, seqlen = 401)
  #cat("Repetition", i, ": ABM-npu\n")
  BsimABMnpu[[i]] <- simNcycleAgebased(sys = br, ctypes = brctypes, compfts = Bsimdata[[i]], n0y0 = br1n0y0,
                                       beta = br1beta, tnowstep = 0.1, hor = 4, seqlen = 401, cycleupdate = FALSE)
  #cat("Repetition", i, ": CM\n")
  BsimCM[[i]] <- simNcycleCorrective(sys = br, ctypes = brctypes, compfts = Bsimdata[[i]], tnowstep = 0.1)
}

# save result lists
str1 <- "Bsim"
str2 <- c("CBMcpu", "CBMepu", "CBMnpu", "ABMepu", "ABMnpu", "CM")
save(list = paste(str1, str2), file = "extsim-b-objects.RData", compression = TRUE)


#