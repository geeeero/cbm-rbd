# --------------------------------------------------------------------------- #
# CBM systems paper - code for example figs - define system with other prior
# --------------------------------------------------------------------------- #

#source("weibull-sys-functions.R")
#source("cbm-sim2.R")
#source("plotdefs.R")

# full system for prior
br <- graph.formula(s -- M -- C1:C2:C3:C4, P1:P2:P3:P4 -- t,
                    C1 -- P1, C2 -- P2, C3 -- P3, C4 -- P4, s -- H -- P3:P4)
brctypes <- list("C"=c("C1", "C2", "C3", "C4"), "H"=c("H"), "M"=c("M"), "P"=c("P1", "P2", "P3", "P4"))
br <- setCompTypes(br, brctypes)
brsysign <- computeSystemSurvivalSignature(br)

#brplotlayout <- matrix(c(0, 1, rep(2,4), rep(3, 4), 4, 1,
#                         0, 0, rep(seq(1.5, -1.5), 2), 0, -1), ncol = 2)
#plot(br, layout = brplotlayout)

# order is C H M P
br1beta <- c(2, 1, 2.5, 1.5)
br1mttf <- c(3, 10, 8, 5)
br1n0 <- c(1, 1, 1, 1) # weak prior
br1n0y0 <- list(C=c(br1n0[1], failuretolambda(br1mttf[1], br1beta[1])),
                H=c(br1n0[2], failuretolambda(br1mttf[2], br1beta[2])),
                M=c(br1n0[3], failuretolambda(br1mttf[3], br1beta[3])),
                P=c(br1n0[4], failuretolambda(br1mttf[4], br1beta[4])))
brfts0 <- list(C=NULL, H=NULL, M=NULL, P=NULL)


#