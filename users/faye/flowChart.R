library(grid)
library(Gmisc)

grid.newpage()

### build the boxes
fishCom <- boxGrob("\n\n  Fish community change  \n through time\n\n", x=0.7, y=0.75, box_gp = gpar(fill = "#abdda4"))
socEc <- boxGrob("\nSocio-economic \nfactors\n", x=0.85, y=0.25,box_gp = gpar(fill = "#d53e4f"))
climCh <- boxGrob("\nEnvironmental factors\n including temperature\n", x=0.2, y=0.85, box_gp = gpar(fill = "#fdae61"))
fish <- boxGrob("\n    Fishing practices    \n", x=0.15, y=0.5, box_gp = gpar(fill = "#3288bd"))
fishMan <- boxGrob("\nFisheries\n management policies\n", x=0.3, y=0.15, box_gp = gpar(fill = "#5e4fa2"))

### get arrows
connectGrob(fishCom, socEc, "vertical", lty_gp = gpar(lwd=1, col="black", fill="black"))
connectGrob(fishCom, fishMan, "vertical", lty_gp = gpar(lwd=2, col="black", fill="black"))
connectGrob(fishCom, fish, "horizontal", lty_gp = gpar(lwd=2, col="black", fill="black"))
connectGrob(fishMan, fishCom, "vertical")
connectGrob(fish, fishCom, "horizontal")
connectGrob(climCh, fish, "vertical")
connectGrob(climCh, fishCom, "horizontal")
connectGrob(socEc, fish, "horizontal", lty_gp = gpar(lwd=2, col="black", fill="black"))
connectGrob(fish, socEc, "horizontal")
connectGrob(fishMan, socEc, "horizontal")
connectGrob(fishMan, fish, "-")
connectGrob(fish, fishMan, "vertical", lty_gp = gpar(lwd=2, col="black", fill="black"))



### show the plot
fishCom
socEc
climCh
fishMan
fish