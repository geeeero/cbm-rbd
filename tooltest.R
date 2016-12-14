# --------------------------------------------------------------------------- #
# code for the excel tool
# --------------------------------------------------------------------------- #

# cbm-sys-other for defining the rbd, setting prior parameters and visualizing them

# failure history (things after tnow are ignored!)
brcompfts <- list(C1 = NA, C2 = 6, C3 = 7, C4 = NA, H = 8,  M = NA, P1 = NA, P2 = 3, P3 = 4, P4 = NA)
# failure history plot
brfstdf <- compfts2df(brcompfts)
plotfts(brcompfts)

# address all arguments!
toolbrt0 <- taustarhist(sys = br, ctypes = brctypes, compfts = brcompfts,
                        n0y0 = br1n0y0, beta = br1beta, tnowvec = 0, hor = 4,
                        seqlen = 401, prior = FALSE, cu = 1, cp = 0.2,
                        onecycle = FALSE, tool = TRUE, verbose = FALSE)
toolbrt0 <- taustarhistprocessor(toolbrt0)
toolplot(toolbrt0)

#
