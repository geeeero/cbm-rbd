# --------------------------------------------------------------------------- #
# cbm systems paper - plot definitions
# --------------------------------------------------------------------------- #

bottomlegend <- theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank())
rightlegend <- theme(legend.title = element_blank())
nolegend <- guides(linetype="none", fill="none", group="none", colour="none")

minirel <- list(theme(axis.title = element_blank()),
                coord_cartesian(ylim = c(0, 1)),
                scale_x_continuous(breaks=c(0, 5, 10, 15, 20)),
                scale_y_continuous(breaks=c(0, 0.5, 1)))
minirelsize <- 1.5

#ijarcols <- c("#b2df8a", "#1f78b4")
ijarcols <- c("#1f78b4", "#b2df8a")
ijarcol <- scale_colour_manual(values = ijarcols)
ijarfill <- scale_fill_manual(values = ijarcols)

tuered <- rgb(0.839,0.000,0.290)
tueblue <- rgb(0.000,0.400,0.800)
tueyellow <- rgb(1.000,0.867,0.000)

tuegreen <- rgb(0.000,0.675,0.510)
tuewarmred <- rgb(0.969,0.192,0.192)
tueorange <- rgb(1.000,0.604,0.000)
tuedarkblue <- rgb(0.063,0.063,0.451)

#