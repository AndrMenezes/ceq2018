rm(list = ls())

# Especifica??es ----------------------------------------------------------

library(fitdistrplus)
library(xtable)
library(qcc)

wd    <- "C:/Users/User/Dropbox/5° Série/Controle Estatístico de Qualidade/trabalho2"
# wd    <- "/home/andrefbm/Dropbox/5° Série/Controle Estatístico de Qualidade/trabalho2"
FF  <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}
setwd(wd)
source(paste0(wd, "/bcc-brcc.R"))


# Function ----------------------------------------------------------------
myppplot <- function(y, nome)
{
  n       <- length(y)
  sy      <- sort(y)
  emp     <- ppoints(n)  
  
  mle1    <- c(mean(y), sd(y))
  mle2    <- coef(fitdist(y, "beta"))
  teo1    <- pnorm(q = sy, mean = mle1[1], sd = mle1[2])  
  teo2    <- pbeta(q = sy, shape1 = mle2[1], shape2 = mle2[2])
  ks1     <- ks.test(x = y, "pnorm", mean = mle1[1], sd = mle1[2])$p.value
  ks2     <- ks.test(x = y, "pbeta", shape1 = mle2[1], shape2 = mle2[2])$p.value
  teo     <- list(teo1, teo2)
  tit     <- c("Normal", "Beta")
  ks      <- c(ks1, ks2)
  for(i in 1:2)
  {
    pdf(paste0("ppplot-", i, "-", nome, ".pdf"))
    par(mar = c(3.2, 3.2, 1.0, 1.0), cex = 1.8)
    plot(x = teo[[i]], y = emp,  type = 'p', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', pch = 3,
         xlim = c(0, 1), ylim = c(0,1), col = 1, bty = 'l')
    box()
    abline(0, 1, lwd = 1)
    abline(h = seq(0, 1, by = 0.25), v = seq(0, 1, by = 0.25), col = "lightgray", lty = "dotted")
    axis(1, seq(0, 1, by = 0.25), FF(seq(0, 1, by = 0.25), 2))
    axis(2, seq(0, 1, by = 0.25), FF(seq(0, 1, by = 0.25), 2))
    mtext("Probabilidades teóricas", side = 1, line = 2.0, cex = 1.8)
    mtext("Probabilidades empíricas", side = 2.1, line =2.2, cex = 1.8)
    mtext(tit[i], side = 3, cex = 1.8)
    legend("topleft", legend = paste0("valor-p = ", FF(ks[i], 4)), bty = "n")
    graphics.off()  
  }
}

# Manufacturing process of frozen orange juice ----------------------------
data("orangejuice")
orangejuice$y <- orangejuice$D / orangejuice$size
orangejuice2  <- subset(orangejuice, subset = orangejuice$trial == TRUE)
y             <- orangejuice2$y
myppplot(y = y, nome = "orange")

# Contaminated peanut by toxic substances in 34 batches -------------------
peanuts <- read.table(file = paste0(wd, "/dados/peanut.txt"), header = T, skip = 4, sep = "\t")
y       <- peanuts$y
myppplot(y = y, nome = "peanuts")

# Ammonia oxidation -------------------------------------------------------
amonia <- read.table(paste0(wd, "/dados/amonia.txt"), sep = "\t", header = T, skip = 7)
y      <- amonia$y
myppplot(y = y, nome = "amonia")

# Transitors data ---------------------------------------------------------
transistors <- read.table(paste0(wd, "/dados/transistors.txt"), sep = "", header = T, skip = 2)
y           <- transistors$x / 1000
myppplot(y = y, nome = "transistors")




