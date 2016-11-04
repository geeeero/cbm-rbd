# --------------------------------------------------------------------------- #
# cbm systems paper - plots for campi slotevent talk
# --------------------------------------------------------------------------- #

source("weibull-sys-functions.R")
source("plotdefs.R")

# component reliability plots
sysA <- graph.formula(s -- A -- t)
betaA <- 2
compA <- list(survsign = computeSystemSurvivalSignature(setCompTypes(sysA, list("A" = "A"))),
              n0y0 = list(A = c(1, failuretolambda(5, betaA))),
              beta = betaA)
compArel0 <-sysrelnowhor(compA$survsign, compA$n0y0, compA$beta, list(A = NULL), tnow = 0, hor = 10, seqlen=101)
compArel5 <-sysrelnowhor(compA$survsign, compA$n0y0, compA$beta, list(A = 5), tnow = 5, hor = 10, seqlen=101)
pdf("talk-campi/compArel0.pdf", width = minirelsize, height = minirelsize)
ggplot(compArel0, aes(x = t, y = rel)) + geom_line() + minirel
dev.off()
pdf("talk-campi/compArel5.pdf", width = minirelsize, height = minirelsize)
ggplot(compArel5, aes(x = t, y = rel)) + geom_line() + minirel
dev.off()

sysB <- graph.formula(s -- B -- t)
betaB <- 1.2
compB <- list(survsign = computeSystemSurvivalSignature(setCompTypes(sysB, list("B" = "B"))),
              n0y0 = list(B = c(1, failuretolambda(8, betaB))),
              beta = betaB)
compBrel0 <-sysrelnowhor(compB$survsign, compB$n0y0, compB$beta, list(B = NULL), tnow = 0, hor = 10, seqlen=101)
compBrel5 <-sysrelnowhor(compB$survsign, compB$n0y0, compB$beta, list(B = NULL), tnow = 5, hor = 10, seqlen=101)
pdf("talk-campi/compBrel0.pdf", width = minirelsize, height = minirelsize)
ggplot(compBrel0, aes(x = t, y = rel)) + geom_line() + minirel
dev.off()
pdf("talk-campi/compBrel5.pdf", width = minirelsize, height = minirelsize)
ggplot(compBrel5, aes(x = t, y = rel)) + geom_line() + minirel
dev.off()


#