# --------------------------------------------------------------------------- #
# CBM systems paper - code for example figs - 3 exemplary failure histories
# using total unit cost rate now!
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
br1ghist2 <- gnowhist(br, brctypes, brcompfts, br1n0y0, br1beta, c(0,2,4,6,8), hor = 10,
                      seqlen = 1001, onecycle = FALSE)
br1ghist2$tnow <- as.factor(br1ghist2$tnow)
br1ghist2a <- br1ghist2
names(br1ghist2a)[c(3,4)] <- c("Reliability", "Unit cost rate")
br1ghist2a <- melt(br1ghist2a, c("tnow", "t"))
br1ghist2a <- subset(br1ghist2a, !is.na(value))
ghist1fig2rT <- ggplot(subset(br1ghist2a, variable == "Reliability"), aes(x = t, y = value)) +
  geom_line(aes(linetype = tnow)) + coord_cartesian(xlim = c(0, 18), ylim = c(0, 1)) +
  xlab(expression(t)) + ylab(expression(paste(R[sys]^(t[now]), (t)))) + guides(linetype = 'none') + 
  scale_x_continuous(breaks = seq(0, 18, by = 2), minor_breaks = 0:18) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), minor_breaks = seq(0, 1, by = 0.1)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ghist1fig2gT <- ggplot(subset(br1ghist2a, variable == "Unit cost rate"), aes(x = t, y = value)) +
  geom_line(aes(linetype = tnow)) + coord_cartesian(xlim = c(0, 18), ylim = c(0, 0.5)) +
  xlab(expression(t)) + ylab(expression(paste(g[sys]^(t[now]), (t)))) +
  guides(linetype = guide_legend(title=expression(t[now]))) +
  theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  scale_x_continuous(breaks = seq(0, 18, by = 2), minor_breaks = 0:18) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5), minor_breaks = seq(0, 3, by = 0.25))
pdf("ghist1fig2T.pdf", width = 6, height = 4)
grid.arrange(ghist1fig2rT, ghist1fig2gT, nrow = 2, heights = c(1, 1.5))
dev.off()


# -------------------------------------------------------------------------------------------------

br1taus1fineT <- taustarhist(br, brctypes, brcompfts, br1n0y0, br1beta, seq(0,8,by=0.1),
                             hor=4, seqlen=401, onecycle = FALSE)
names(br1taus1fineT)[2:5] <- c(taustar = expression(tau['*']^(t[now])),
                               tstar = expression(t['*']^(t[now])),
                               cstar = expression(g['*']^(t[now])),
                               ctotal = expression(g[total]^(t[now])))
tauhist1fig2T <- ggplot(melt(br1taus1fineT, "tnow"), aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(group = variable))  + geom_point(aes(group = variable), size = 0.15) +
  facet_wrap(~ variable, nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(axis.title.y = element_blank())
pdf("tauhist1fig2T.pdf", width = 6, height = 4)
tauhist1fig2T
dev.off()
#br1taus1fine[40:41,]
br1taus1fineprT <- taustarhist(br, brctypes, brcompfts, br1n0y0, br1beta, seq(0,8,by=0.1),
                               hor=4, seqlen=401, prior = TRUE, onecycle = FALSE)
names(br1taus1fineprT)[2:5] <- c(taustar = expression(tau['*']^(t[now])),
                                 tstar = expression(t['*']^(t[now])),
                                 cstar = expression(g['*']^(t[now])),
                                 ctotal = expression(g[total]^(t[now])))
br1taus1prpoT <- rbind(data.frame(melt(br1taus1fineT, "tnow"), variable2 = "Parameter update"),
                       data.frame(melt(br1taus1fineprT, "tnow"), variable2 = "No parameter update"))
tauhist1fig3T <- ggplot(br1taus1prpoT, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable2), size = 0.15) + #theme_bw() +
  facet_wrap(~ variable, nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(axis.title.y = element_blank()) + bottomlegend + ijarcol
pdf("tauhist1fig3T.pdf", width = 6, height = 4.5)
tauhist1fig3T
dev.off()


# -------------------------------------------------------------------------------------------------

# earlier failures - are more or less like expected???
brcompfts2 <- list(C1 = NA, C2 = 1, C3 = 2, C4 = NA, H = 8,  M = NA, P1 = NA, P2 = 0.5, P3 = 1.5, P4 = NA)
plotfts(brcompfts2, 2)
br1taus2fineT <- taustarhist(br, brctypes, brcompfts2, br1n0y0, br1beta, seq(0, 2, by=0.05),
                             hor=4, seqlen=401, onecycle = FALSE)
names(br1taus2fineT)[2:5] <- c(taustar = expression(tau['*']^(t[now])),
                               tstar = expression(t['*']^(t[now])),
                               cstar = expression(g['*']^(t[now])),
                               ctotal = expression(g[total]^(t[now])))
br1taus2fineprT <- taustarhist(br, brctypes, brcompfts2, br1n0y0, br1beta, seq(0, 2, by=0.05),
                               hor=4, seqlen=401, prior=T, onecycle = FALSE)
names(br1taus2fineprT)[2:5] <- c(taustar = expression(tau['*']^(t[now])),
                                 tstar = expression(t['*']^(t[now])),
                                 cstar = expression(g['*']^(t[now])),
                                 ctotal = expression(g[total]^(t[now])))
br1taus2prpoT <- rbind(data.frame(melt(br1taus2fineT, "tnow"), variable2 = "Parameter update"),
                       data.frame(melt(br1taus2fineprT, "tnow"), variable2 = "No parameter update"))
tauhist2fig3T <- ggplot(br1taus2prpoT, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable2), size = 0.15) +
  facet_wrap(~ variable, nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(axis.title.y = element_blank()) + bottomlegend+ ijarcol
pdf("tauhist2fig3T.pdf", width = 6, height = 4.5)
tauhist2fig3T
dev.off()

# -------------------------------------------------------------------------------------------------

# later failures *** do extremely early failures instead???
brcompfts3 <- list(C1 = NA, C2 = 10, C3 = 12, C4 = NA, H = 14,  M = NA, P1 = NA, P2 = 8, P3 = 9, P4 = NA)
plotfts(brcompfts3, 14)
br1taus3fineT <- taustarhist(br, brctypes, brcompfts3, br1n0y0, br1beta, seq(0, 14, by=0.2),
                             hor=4, seqlen=401, onecycle = FALSE)
names(br1taus3fineT)[2:5] <- c(taustar = expression(tau['*']^(t[now])),
                               tstar = expression(t['*']^(t[now])),
                               cstar = expression(g['*']^(t[now])),
                               ctotal = expression(g[total]^(t[now])))
br1taus3fineprT <- taustarhist(br, brctypes, brcompfts3, br1n0y0, br1beta, seq(0, 14, by=0.2),
                               hor=4, seqlen=401, prior=T, onecycle = FALSE)
names(br1taus3fineprT)[2:5] <- c(taustar = expression(tau['*']^(t[now])),
                                 tstar = expression(t['*']^(t[now])),
                                 cstar = expression(g['*']^(t[now])),
                                 ctotal = expression(g[total]^(t[now])))
br1taus3prpoT <- rbind(data.frame(melt(br1taus3fineT, "tnow"), variable2 = "Parameter update"),
                       data.frame(melt(br1taus3fineprT, "tnow"), variable2 = "No parameter update"))
tauhist3fig3T <- ggplot(br1taus3prpoT, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + geom_point(aes(colour = variable2), size = 0.15) +
  scale_x_continuous(breaks = seq(0, 14, by = 4), minor_breaks = seq(0, 14, by = 2)) +
  facet_wrap(~ variable, nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(axis.title.y = element_blank()) + bottomlegend + ijarcol
pdf("tauhist3fig3T.pdf", width = 6, height = 4.5)
tauhist3fig3T
dev.off()

#br1taus3fineprT[(br1taus3fineprT$tnow > 8.2) & (br1taus3fineprT$tnow < 9.2),]

#
