# --------------------------------------------------------------------------- #
# CBM systems paper - code for example figs - 3 exemplary failure histories
# --------------------------------------------------------------------------- #

source("cbm-sys-other.R")

# exemplary failure history
brcompfts <- list(C1 = NA, C2 = 6, C3 = 7, C4 = NA, H = 8,  M = NA, P1 = NA, P2 = 3, P3 = 4, P4 = NA)
# failure history plot
brfstdf <- compfts2df(brcompfts)
plotfts(brcompfts)
#plotfts(brcompfts,10) 
#plotfts(brcompfts,6) 

# gnow & sysrelnow histories
br1ghist1 <- gnowhist(br, brctypes, brcompfts, br1n0y0, br1beta, c(0,2,4,6,8), hor = 10, seqlen = 1001)
br1ghist1$tnow <- as.factor(br1ghist1$tnow)
br1ghist1a <- br1ghist1
names(br1ghist1a)[c(3,4)] <- c("Reliability", "Unit cost rate")
br1ghist1a <- melt(br1ghist1a, c("tnow", "t"))
br1ghist1a <- subset(br1ghist1a, !is.na(value))
ghist1fig2r <- ggplot(subset(br1ghist1a, variable == "Reliability"), aes(x = t, y = value)) +
  geom_line(aes(linetype = tnow)) + coord_cartesian(xlim = c(0, 18), ylim = c(0, 1)) +
  xlab(expression(t)) + ylab(expression(paste(R[sys]^(t[now]), (t)))) + guides(linetype = 'none') + 
  scale_x_continuous(breaks = seq(0, 18, by = 2), minor_breaks = 0:18) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), minor_breaks = seq(0, 1, by = 0.1)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ghist1fig2g <- ggplot(subset(br1ghist1a, variable == "Unit cost rate"), aes(x = t, y = value)) +
  geom_line(aes(linetype = tnow)) + coord_cartesian(xlim = c(0, 18), ylim = c(0, 2.2)) +
  xlab(expression(t)) + ylab(expression(paste(g[sys]^(t[now]), (t)))) +
  guides(linetype = guide_legend(title=expression(t[now]))) +
  theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  scale_x_continuous(breaks = seq(0, 18, by = 2), minor_breaks = 0:18) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5), minor_breaks = seq(0, 3, by = 0.25))
pdf("ghist1fig2.pdf", width = 6, height = 4)
grid.arrange(ghist1fig2r, ghist1fig2g, nrow = 2, heights = c(1, 1.5))
dev.off()


# -------------------------------------------------------------------------------------------------

# tauhist figures

br1taus1fine <- taustarhist(br, brctypes, brcompfts, br1n0y0, br1beta, seq(0,8,by=0.1), hor=4, seqlen=401)
names(br1taus1fine)[2:6] <- c(taustar = expression(tau['*']^(t[now])),
                              tstar = expression(t['*']^(t[now])),
                              cstar = expression(g['*']^(t[now])),
                              ctotal = expression(g[total]^(t[now])),
                              relstar = expression(paste(R[sys]^(t[now]),(t['*']^(t[now])))))
tauhist1fig2 <- ggplot(melt(br1taus1fine[-6], "tnow"), aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(group = variable))  + geom_point(aes(group = variable), size = 0.15) +
  facet_wrap(~ variable, nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(axis.title.y = element_blank())
pdf("tauhist1fig2.pdf", width = 6, height = 4)
tauhist1fig2
dev.off()
#br1taus1fine[40:41,]

tauhist1fig2R <- ggplot(melt(br1taus1fine[-5], "tnow"), aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(group = variable))  + geom_point(aes(group = variable), size = 0.15) +
  facet_wrap(~ variable, nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(axis.title.y = element_blank())
pdf("tauhist1fig2R.pdf", width = 6, height = 4)
tauhist1fig2R
dev.off()

br1taus1finepr <- taustarhist(br, brctypes, brcompfts, br1n0y0, br1beta, seq(0,8,by=0.1), hor=4, seqlen=401, prior=T)
names(br1taus1finepr)[2:6] <- c(taustar = expression(tau['*']^(t[now])),
                                tstar = expression(t['*']^(t[now])),
                                cstar = expression(g['*']^(t[now])),
                                ctotal = expression(g[total]^(t[now])),
                                relstar = expression(paste(R[sys]^(t[now]),(t['*']^(t[now])))))
br1taus1prpo <- rbind(data.frame(melt(br1taus1fine[-6], "tnow"), variable2 = "Parameter update"),
                      data.frame(melt(br1taus1finepr[-6], "tnow"), variable2 = "No parameter update"))
tauhist1fig3 <- ggplot(br1taus1prpo, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable2), size = 0.15) + #theme_bw() +
  facet_wrap(~ variable, nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(axis.title.y = element_blank()) + bottomlegend + ijarcol
pdf("tauhist1fig3.pdf", width = 6, height = 4.5)
tauhist1fig3
dev.off()

br1taus1prpoR <- rbind(data.frame(melt(br1taus1fine[-5], "tnow"), variable2 = "Parameter update"),
                       data.frame(melt(br1taus1finepr[-5], "tnow"), variable2 = "No parameter update"))
tauhist1fig3R <- ggplot(br1taus1prpoR, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable2), size = 0.15) + #theme_bw() +
  facet_wrap(~ variable, nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(axis.title.y = element_blank()) + bottomlegend + ijarcol
pdf("tauhist1fig3R.pdf", width = 6, height = 4.5)
tauhist1fig3R
dev.off()

# -------------------------------------------------------------------------

# earlier failures = fitting failures
brcompfts2 <- list(C1 = NA, C2 = 1, C3 = 2, C4 = NA, H = 8,  M = NA, P1 = NA, P2 = 0.5, P3 = 1.5, P4 = NA)
plotfts(brcompfts2, 2)
br1taus2fine <- taustarhist(br, brctypes, brcompfts2, br1n0y0, br1beta, seq(0, 2, by=0.05), hor=4, seqlen=401)
names(br1taus2fine)[2:6] <- c(taustar = expression(tau['*']^(t[now])),
                              tstar = expression(t['*']^(t[now])),
                              cstar = expression(g['*']^(t[now])),
                              ctotal = expression(g[total]^(t[now])),
                              relstar = expression(paste(R[sys]^(t[now]),(t['*']^(t[now])))))
br1taus2finepr <- taustarhist(br, brctypes, brcompfts2, br1n0y0, br1beta, seq(0, 2, by=0.05), hor=4, seqlen=401, prior=T)
names(br1taus2finepr)[2:6] <- c(taustar = expression(tau['*']^(t[now])),
                                tstar = expression(t['*']^(t[now])),
                                cstar = expression(g['*']^(t[now])),
                                ctotal = expression(g[total]^(t[now])),
                                relstar = expression(paste(R[sys]^(t[now]),(t['*']^(t[now])))))
br1taus2prpo <- rbind(data.frame(melt(br1taus2fine[-6], "tnow"), variable2 = "Parameter update"),
                      data.frame(melt(br1taus2finepr[-6], "tnow"), variable2 = "No parameter update"))
tauhist2fig3 <- ggplot(br1taus2prpo, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable2), size = 0.15) +
  facet_wrap(~ variable, nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(axis.title.y = element_blank()) + bottomlegend+ ijarcol
pdf("tauhist2fig3.pdf", width = 6, height = 4.5)
tauhist2fig3
dev.off()

br1taus2prpoR <- rbind(data.frame(melt(br1taus2fine[-5], "tnow"), variable2 = "Parameter update"),
                       data.frame(melt(br1taus2finepr[-5], "tnow"), variable2 = "No parameter update"))
tauhist2fig3R <- ggplot(br1taus2prpoR, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable2), size = 0.15) +
  facet_wrap(~ variable, nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(axis.title.y = element_blank()) + bottomlegend+ ijarcol
pdf("tauhist2fig3R.pdf", width = 6, height = 4.5)
tauhist2fig3R
dev.off()

# much later failures (not in paper)
brcompfts3 <- list(C1 = NA, C2 = 10, C3 = 12, C4 = NA, H = 14,  M = NA, P1 = NA, P2 = 8, P3 = 9, P4 = NA)
plotfts(brcompfts3, 14)
br1taus3fine <- taustarhist(br, brctypes, brcompfts3, br1n0y0, br1beta, seq(0, 14, by=0.2), hor=4, seqlen=401)
names(br1taus3fine)[2:5] <- c(taustar = expression(tau['*']^(t[now])),
                              tstar = expression(t['*']^(t[now])),
                              cstar = expression(g['*']^(t[now])),
                              ctotal = expression(g[total]^(t[now])))
br1taus3finepr <- taustarhist(br, brctypes, brcompfts3, br1n0y0, br1beta, seq(0, 14, by=0.2), hor=4, seqlen=401, prior=T)
names(br1taus3finepr)[2:5] <- c(taustar = expression(tau['*']^(t[now])),
                                tstar = expression(t['*']^(t[now])),
                                cstar = expression(g['*']^(t[now])),
                                ctotal = expression(g[total]^(t[now])))
br1taus3prpo <- rbind(data.frame(melt(br1taus3fine, "tnow"), variable2 = "Parameter update"),
                      data.frame(melt(br1taus3finepr, "tnow"), variable2 = "No parameter update"))
tauhist3fig3 <- ggplot(br1taus3prpo, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + geom_point(aes(colour = variable2), size = 0.15) +
  scale_x_continuous(breaks = seq(0, 14, by = 4), minor_breaks = seq(0, 14, by = 2)) +
  facet_wrap(~ variable, nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(axis.title.y = element_blank()) + bottomlegend + ijarcol
pdf("tauhist3fig3.pdf", width = 6, height = 4.5)
tauhist3fig3
dev.off()

br1taus3finepr[(br1taus3finepr$tnow > 8.2) & (br1taus3finepr$tnow < 9.2),]

# extremely early failures
brcompfts4 <- list(C1 = NA, C2 = 0.2, C3 = 0.3, C4 = NA, H = 0.5,  M = NA, P1 = NA, P2 = 0.05, P3 = 0.15, P4 = NA)
plotfts(brcompfts4, 0.5)
br1taus4fine <- taustarhist(br, brctypes, brcompfts4, br1n0y0, br1beta, seq(0, 0.5, by=0.01), hor=4, seqlen=401)
names(br1taus4fine)[2:6] <- c(taustar = expression(tau['*']^(t[now])),
                              tstar = expression(t['*']^(t[now])),
                              cstar = expression(g['*']^(t[now])),
                              ctotal = expression(g[total]^(t[now])),
                              relstar = expression(paste(R[sys]^(t[now]),(t['*']^(t[now])))))
br1taus4finepr <- taustarhist(br, brctypes, brcompfts4, br1n0y0, br1beta, seq(0, 0.5, by=0.01), hor=4, seqlen=401, prior=T)
names(br1taus4finepr)[2:6] <- c(taustar = expression(tau['*']^(t[now])),
                                tstar = expression(t['*']^(t[now])),
                                cstar = expression(g['*']^(t[now])),
                                ctotal = expression(g[total]^(t[now])),
                                relstar = expression(paste(R[sys]^(t[now]),(t['*']^(t[now])))))
br1taus4prpo <- rbind(data.frame(melt(br1taus4fine[-6], "tnow"), variable2 = "Parameter update"),
                      data.frame(melt(br1taus4finepr[-6], "tnow"), variable2 = "No parameter update"))
tauhist4fig3 <- ggplot(br1taus4prpo, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + geom_point(aes(colour = variable2), size = 0.15) +
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.1), minor_breaks = seq(0, 0.5, by = 0.05)) +
  facet_wrap(~ variable, nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(axis.title.y = element_blank()) + bottomlegend + ijarcol
pdf("tauhist4fig3.pdf", width = 6, height = 4.5)
tauhist4fig3
dev.off()

br1taus4prpoR <- rbind(data.frame(melt(br1taus4fine[-5], "tnow"), variable2 = "Parameter update"),
                       data.frame(melt(br1taus4finepr[-5], "tnow"), variable2 = "No parameter update"))
tauhist4fig3R <- ggplot(br1taus4prpoR, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + geom_point(aes(colour = variable2), size = 0.15) +
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.1), minor_breaks = seq(0, 0.5, by = 0.05)) +
  facet_wrap(~ variable, nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(axis.title.y = element_blank()) + bottomlegend + ijarcol
pdf("tauhist4fig3R.pdf", width = 6, height = 4.5)
tauhist4fig3R
dev.off()



#