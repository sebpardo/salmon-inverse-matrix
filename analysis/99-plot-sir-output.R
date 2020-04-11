# Function to plot SIR output

plot_sir_output <- function(siroutput, siryears, posteriors, simuldat, S1, S2, Prg) {

old.par <- par(no.readonly = TRUE)
  
par(mfcol = c(length(siroutput),3), mar = c(0,1,0,1.3), oma = c(6, 4, 4.5, 1))
#lapply(posteriors, function (x) 
for (i in 1:length(posteriors))  {
  ds <- density(posteriors[[i]]$S1)
  plot(ds, main = "", xlim = c(-0.02, 0.4), xaxs = "i", axes = FALSE, bty = "l", col = "red")
  axis(1, labels = FALSE, tck = 0.05)
  #  axis(2, labels = FALSE, tck = 0)
  lines(density(S1), col = "#55555570")
  abline(v = simuldat[simuldat$year == siryears[i], "S1_true"], col = "blue")
  qs <- quantile(posteriors[[i]]$S1, c(0.05, 0.25, 0.5, 0.75, 0.95))
  yp <- max(ds$y)/4
  segments(qs[1],yp,qs[5],yp, col = "red")
  segments(qs[1],yp,qs[5],yp, col = "darkred")
  segments(qs[2],yp,qs[4],yp, lwd =2, col = "darkred")
  points(qs[3], yp, col = "darkred")
  mtext(siryears[i], cex = 1.2, line = -2.5, adj = -0.3) # adj =-0.3 - -0.4)
  
}
axis(1, tck = 0.05, cex.axis = 1.3)

for (i in 1:length(posteriors))  {
  ds <- density(posteriors[[i]]$S2)
  plot(ds, main = "", xlim = c(-0.02, 0.9), xaxs = "i",axes = FALSE, bty = "l", col = "red")
  axis(1, labels = FALSE, tck = 0.05)
  #  axis(2, labels = FALSE, tck = 0)
  lines(density(S2), col = "#55555570")
  abline(v = simuldat[simuldat$year == siryears[i], "S2_true"], col = "blue")
  qs <- quantile(posteriors[[i]]$S2, c(0.05, 0.25, 0.5, 0.75, 0.95))
  yp <- max(ds$y)/4
  segments(qs[1],yp,qs[5],yp, col = "red")
  segments(qs[1],yp,qs[5],yp, col = "darkred")
  segments(qs[2],yp,qs[4],yp, lwd =2, col = "darkred")
  points(qs[3], yp, col = "darkred")
  #mtext(siryears[i], line = -1.5, adj =0.9, cex = 0.8)
}
axis(1, tck = 0.05, cex.axis = 1.3)

for (i in 1:length(posteriors))  {
  ds <- density(posteriors[[i]]$Prg)
  plot(ds, main = "", xlim = c(-0.02, 1.05), xaxs = "i",axes = FALSE, bty = "l", col = "red")
  axis(1, labels = FALSE, tck = 0.05)
  #  axis(2, labels = FALSE, tck = 0)
  lines(density(Prg), col = "#55555570")
  abline(v = simuldat[simuldat$year == siryears[i], "Pr_true"], col = "blue")
  qs <- quantile(posteriors[[i]]$Prg, c(0.05, 0.25, 0.5, 0.75, 0.95))
  yp <- max(ds$y)/4
  segments(qs[1],yp,qs[5],yp, col = "red")
  segments(qs[1],yp,qs[5],yp, col = "darkred")
  segments(qs[2],yp,qs[4],yp, lwd =2, col = "darkred")
  points(qs[3], yp, col = "darkred")
  #mtext(siryears[i], line = -1.5, adj =0.9, cex = 0.8)
}
axis(1, tck = 0.05, cex.axis = 1.3)

# Three columns
mtext(text = "Smolt to first year\n at sea survival (S1)",
      outer = TRUE, side = 1, adj = 0.07, line = 4)
mtext(text = "Survival probability of\n 2nd year at sea (S2)",
      outer = TRUE, side = 1, adj = 0.5, line = 4)
mtext(text = "Proportion returning\n as grilse (Prg)",
      outer = TRUE, side = 1, adj = 0.93, line = 4)
# 
mtext(text = "Marine surivival estimates using a\nsample-importance-resampling algorithm", 
      outer = TRUE, side = 3, cex = 1.3)
mtext(text = paste("Scenario", parse_number(deparse(substitute(simuldat))), "\n"),
      outer = TRUE, side = 3, cex = 1.3, adj = 0)
on.exit(par(old.par))
}
