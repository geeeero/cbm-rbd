# component prior predictive reliability functions to visualize prior assumptions

source("cbm-sys-other.R")

comphor <- 10; complen <- 201
onecomp <- graph.formula(s -- 1 -- t)
Ccomp <- setCompTypes(onecomp, list(C = "1"))
Hcomp <- setCompTypes(onecomp, list(H = "1"))
Mcomp <- setCompTypes(onecomp, list(M = "1"))
Pcomp <- setCompTypes(onecomp, list(P = "1"))
Crel1 <- sysrelnowhor(computeSystemSurvivalSignature(Ccomp), br1n0y0[1], br1beta[1], brfts0[1],
                      tnow = 0, hor = comphor, seqlen = complen)
Hrel1 <- sysrelnowhor(computeSystemSurvivalSignature(Hcomp), br1n0y0[2], br1beta[2], brfts0[2],
                      tnow = 0, hor = comphor, seqlen = complen)
Mrel1 <- sysrelnowhor(computeSystemSurvivalSignature(Mcomp), br1n0y0[3], br1beta[3], brfts0[3],
                      tnow = 0, hor = comphor, seqlen = complen)
Prel1 <- sysrelnowhor(computeSystemSurvivalSignature(Pcomp), br1n0y0[4], br1beta[4], brfts0[4],
                      tnow = 0, hor = comphor, seqlen = complen)
compdf1 <- rbind(data.frame(Crel1, comp = "C"), data.frame(Hrel1, comp = "H"),
                 data.frame(Mrel1, comp = "M"), data.frame(Prel1, comp = "P"))
# in different facets
compprior1fig1 <- ggplot(compdf1, aes(x = t, y = rel)) + geom_line(aes(group = comp)) + 
  ylab("Reliability") + scale_x_continuous(breaks=seq(0, 10, by=2), minor_breaks=0:10) +
  facet_wrap(~ comp, nrow = 2, scales = "free_y") + coord_cartesian(ylim = c(0, 1))
pdf("compprior1fig1.pdf", width = 6, height = 4)
print(compprior1fig1)
dev.off()
# as different linetypes
compprior1fig2 <- ggplot(compdf1, aes(x = t, y = rel)) + geom_line(aes(linetype = comp)) + 
  ylab("Reliability") + scale_x_continuous(breaks=seq(0, 10, by=2), minor_breaks=0:10) +
  guides(linetype = guide_legend(title=NULL)) + coord_cartesian(ylim = c(0, 1))
pdf("compprior1fig2.pdf", width = 6, height = 4)
print(compprior1fig2)
dev.off()

Crel1[which.min(abs(Crel1$rel - 0.5)),]
Hrel1[which.min(abs(Hrel1$rel - 0.5)),]
Mrel1[which.min(abs(Mrel1$rel - 0.5)),]
Prel1[which.min(abs(Prel1$rel - 0.5)),]
