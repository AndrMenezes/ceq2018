rm(list = ls())

# Especifica??es ----------------------------------------------------------

library(fitdistrplus)
library(mle.tools)
library(qcc)
library(xtable)

# wd    <- "C:/Users/User/Dropbox/5? S?rie/Controle Estat?stico de Qualidade/trabalho2"
wd    <- "/home/andrefbm/Dropbox/5° Série/Controle Estatístico de Qualidade/trabalho2"
cores <- c("#E41A1C", "#377EB8", "#4DAF4A",  "#FF7F00", "#A65628")
setwd(wd)
source(paste0(wd, "/bcc-brcc.R"))

# Manufacturing process of frozen orange juice ----------------------------
data("orangejuice")
orangejuice$y <- orangejuice$D / orangejuice$size
orangejuice2  <- subset(orangejuice, subset = orangejuice$trial == TRUE)

## Calculating the control limits

CL.bcc        <- bcc(y = orangejuice2$y, alpha = 0.05)
CL.bcc.bias   <- bcc(y = orangejuice2$y, bias = TRUE, alpha = 0.05)
CL.shewhart   <- pcontrol(prop = orangejuice2$y, n = 50, correction = "none")
CL.ryan       <- pcontrol(prop = orangejuice2$y, n = 50, correction = "ryan")
CL.chen       <- pcontrol(prop = orangejuice2$y, n = 50, correction = "chen")
CL.joekes     <- pcontrol(prop = orangejuice2$y, n = 50, correction = "joekes")
CLs           <- matrix(c(CL.bcc, CL.bcc.bias, CL.shewhart, CL.ryan, CL.chen, CL.joekes), ncol = 3, byrow = T)
rownames(CLs) <- c("BCC", "BCC com corre??o", "Shewhart", "Ryan", "Chen", "Joekes")

## Plotting the control limits and the data
y  <- orangejuice2$y
t  <- 1:length(y)

for(i in 1:nrow(CLs))
{
  CL        <- CLs[i, ]
  R         <- range(y, CL)
  twarnings <- which(y < CL[1] | y > CL[3])
  ywarnings <- y[twarnings]
  
  pdf(paste0("orangeCL-", i, ".pdf"))
  par(mar =  c(3.2, 3.2, 1.0, 0.5), cex = 1.6)
  plot(t, y, type = "o", ylim = R, cex = 0.8, pch = 15, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  abline(h = CL, col = "#0080ff", lty = 2, lwd = 2)
  points(x = twarnings, y = ywarnings, col = "red", pch = 15)
  axis(1, seq(min(t) +1, max(t), by = 2), FF(seq(min(t) + 1, max(t), by = 2), 0))
  axis(2, seq(R[1], R[2], length.out = 5),FF(seq(R[1], R[2], length.out = 5), 2))
  mtext("", side = 1, line = 2.0, cex = 1.8)
  mtext("Propor??o de embalagens n?o conformes", side = 2.1, line = 2.2, cex = 1.8)
  mtext(rownames(CLs)[i], side = 3, cex = 1.8)
  graphics.off()
}

tab           <- FF(CLs)
colnames(tab) <- c("LCL", "CL" ,"UCL")
print(xtable(tab), include.rownames = T, include.colnames = T)


# Contaminated peanut by toxic substances in 34 batches -------------------
peanuts <- read.table(file = paste0(wd, "/dados/peanut.txt"), header = T, skip = 4, sep = "\t")
head(peanuts)

## Calculating the control limits
CL.bcc        <- bcc(y = peanuts$y, alpha = 0.05)
CL.shewhart   <- pcontrol(prop = peanuts$y, n = 34, correction = "none")
CL.ryan       <- pcontrol(prop = peanuts$y, n = 34, correction = "ryan")
CL.chen       <- pcontrol(prop = peanuts$y, n = 34, correction = "chen")
CL.joekes     <- pcontrol(prop = peanuts$y, n = 34, correction = "joekes")
# res.brcc      <- brcc(y = peanuts$y, X = peanuts$x, Z = peanuts$x, alpha = 0.05, diagnostic = F)

CLs           <- matrix(c(CL.bcc, CL.shewhart, CL.ryan, CL.chen, CL.joekes), ncol = 3, byrow = T)
rownames(CLs) <- c("BCC", "Shewhart", "Ryan", "Chen", "Joekes")

## Plotting the control limits and the data
y  <- peanuts$y
t  <- 1:length(y)

for(i in 1:nrow(CLs))
{
  CL        <- CLs[i, ]
  R         <- range(y, CL)
  twarnings <- which(y < CL[1] | y > CL[3])
  ywarnings <- y[twarnings]
  
  pdf(paste0("peanutCL-", i, ".pdf"))
  par(mar =  c(3.2, 3.2, 1.0, 0.5), cex = 1.6)
  plot(t, y, type = "o", ylim = R, cex = 0.8, pch = 15, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  abline(h = CL, col = "#0080ff", lty = 2, lwd = 2)
  points(x = twarnings, y = ywarnings, col = "red", pch = 15)
  axis(1, seq(min(t) +1, max(t), by = 2), FF(seq(min(t) + 1, max(t), by = 2), 0))
  axis(2, seq(R[1], R[2], length.out = 5),FF(seq(R[1], R[2], length.out = 5), 4))
  mtext("", side = 1, line = 2.0, cex = 1.8)
  mtext("Propor??o de amendoins n?o contamindados", side = 2.1, line = 2.2, cex = 1.8)
  mtext(rownames(CLs)[i], side = 3, cex = 1.8)
  graphics.off()
}

tab           <- FF(CLs)
colnames(tab) <- c("LCL", "CL" ,"UCL")
print(xtable(tab), include.rownames = T, include.colnames = T)

# Ammonia oxidation -------------------------------------------------------
amonia <- read.table(paste0(wd, "/dados/amonia.txt"), sep = "\t", header = T, skip = 7)
head(amonia)

## Calculating the control limits
CL.bcc        <- bcc(y = amonia$y, alpha = 0.05)
CL.shewhart   <- pcontrol(prop = amonia$y, n = 21, correction = "none")
CL.ryan       <- pcontrol(prop = amonia$y, n = 21, correction = "ryan")
CL.chen       <- pcontrol(prop = amonia$y, n = 21, correction = "chen")
CL.joekes     <- pcontrol(prop = amonia$y, n = 21, correction = "joekes")
CL.seno       <- pcontrol(prop = amonia$y, n = 21, correction = "seno")

CLs           <- matrix(c(CL.bcc, CL.shewhart, CL.ryan, CL.chen, CL.joekes), ncol = 3, byrow = T)
rownames(CLs) <- c("BCC", "Shewhart", "Ryan", "Chen", "Joekes")

## Plotting the control limits and the data
y  <- amonia$y
t  <- 1:length(y)

for(i in 1:nrow(CLs))
{
  CL        <- CLs[i, ]
  R         <- range(y, CL)
  twarnings <- which(y < CL[1] | y > CL[3])
  ywarnings <- y[twarnings]
  
  pdf(paste0("amoniaCL-", i, ".pdf"))
  par(mar =  c(3.2, 3.2, 1.0, 0.5), cex = 1.6)
  plot(t, y, type = "o", ylim = R, cex = 0.8, pch = 15, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  abline(h = CL, col = "#0080ff", lty = 2, lwd = 2)
  points(x = twarnings, y = ywarnings, col = "red", pch = 15)
  axis(1, seq(min(t) +1, max(t), by = 2), FF(seq(min(t) + 1, max(t), by = 2), 0))
  axis(2, seq(R[1], R[2], length.out = 5),FF(seq(R[1], R[2], length.out = 5), 4))
  mtext("", side = 1, line = 2.0, cex = 1.8)
  mtext("Propor??o de am?nia n?o convertida", side = 2.1, line = 2.2, cex = 1.8)
  mtext(rownames(CLs)[i], side = 3, cex = 1.8)
  graphics.off()
}

tab           <- FF(CLs)
colnames(tab) <- c("LCL", "CL" ,"UCL")
print(xtable(tab), include.rownames = T, include.colnames = T)


# Transitors data ---------------------------------------------------------
transistors   <- read.table(paste0(wd, "/dados/transistors.txt"), sep = "", header = T, skip = 2)
transistors$y <- transitors$x / 1000

## Calculating the control limits
CL.bcc        <- bcc(y = transistors$y, alpha = 0.05)
CL.shewhart   <- pcontrol(prop = transistors$y, n = 1000, correction = "none")
CL.ryan       <- pcontrol(prop = transistors$y, n = 1000, correction = "ryan")
CL.chen       <- pcontrol(prop = transistors$y, n = 1000, correction = "chen")
CL.joekes     <- pcontrol(prop = transistors$y, n = 1000, correction = "joekes")
CL.seno       <- pcontrol(prop = transistors$y, n = 1000, correction = "seno")

CLs           <- matrix(c(CL.bcc, CL.shewhart, CL.ryan, CL.chen, CL.joekes), ncol = 3, byrow = T)
rownames(CLs) <- c("BCC", "Shewhart", "Ryan", "Chen", "Joekes")

## Plotting the control limits and the data
y  <- transistors$y
t  <- 1:length(y)

for(i in 1:nrow(CLs))
{
  CL        <- CLs[i, ]
  R         <- range(y, CL)
  twarnings <- which(y < CL[1] | y > CL[3])
  ywarnings <- y[twarnings]
  
  pdf(paste0("transistorsCL-", i, ".pdf"))
  par(mar =  c(3.2, 3.2, 1.0, 0.5), cex = 1.6)
  plot(t, y, type = "o", ylim = R, cex = 0.8, pch = 15, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  abline(h = CL, col = "#0080ff", lty = 2, lwd = 2)
  points(x = twarnings, y = ywarnings, col = "red", pch = 15)
  axis(1, seq(min(t) +1, max(t), by = 2), FF(seq(min(t) + 1, max(t), by = 2), 0))
  axis(2, seq(R[1], R[2], length.out = 5),FF(seq(R[1], R[2], length.out = 5), 4))
  mtext("", side = 1, line = 2.0, cex = 1.8)
  mtext("Propor??o de transistores n?o conformes", side = 2.1, line = 2.2, cex = 1.8)
  mtext(rownames(CLs)[i], side = 3, cex = 1.8)
  graphics.off()
}

tab           <- FF(CLs)
colnames(tab) <- c("LCL", "CL" ,"UCL")
print(xtable(tab), include.rownames = T, include.colnames = T)

# Circuit boards data -----------------------------------------------------
data("circuit")
circuit$y <- circuit$x / circuit$size
circuit2  <- subset(circuit, subset = circuit$trial == TRUE)

## Calculating the control limits
CL.bcc        <- bcc(y = circuit2$y, alpha = 0.05)
CL.bcc.bias   <- bcc(y = circuit2$y, bias = TRUE, alpha = 0.05)
CL.shewhart   <- pcontrol(prop = circuit2$y, n = 100, correction = "none")
CL.ryan       <- pcontrol(prop = circuit2$y, n = 100, correction = "ryan")
CL.chen       <- pcontrol(prop = circuit2$y, n = 100, correction = "chen")
CL.joekes     <- pcontrol(prop = circuit2$y, n = 100, correction = "joekes")



# Tire manufacturing process ----------------------------------------------
tire    <- read.table(paste0(wd, "/dados/tire.csv"), sep = ",", header = T)
tire$I1 <- tire$x1 * tire$x2
tire$I2 <- tire$x1 * tire$x4
tire$I3 <- tire$x2 * tire$x5
head(tire)

## Calculating the control limits
res.brcc    <- brcc(y = tire$y, X = tire[, c("x1", "x2", "I1", "I2", "I3")], Z = tire[, c("x1", "I1")])
CL.brcc     <- res.brcc$crtlim
CL.bcc      <- bcc(y = tire$y)
CL.bcc.bias <- bcc(y = tire$y, bias = TRUE)
CL.rcc      <- rcc(model = y ~ x1 + x2 + I1 + I2 + I3, data = tire)



## Plotting the control limits and the data
y         <- tire$y
t         <- 1:length(y)
R         <- range(y, CL.bcc, CL.brcc, CL.rcc)
twarnings <- which(y < CL.brcc[, 1] | y > CL.brcc[, 2])
ywarnings <- y[twarnings]

pdf("tireCL.pdf", width = 9)
par(mar =  c(3.2, 3.2, 2.0, 0.5), cex = 1.6)
plot(t, y, type = "p", ylim = R, cex = 0.8, pch = 15, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
lines(t, CL.brcc[, 1], lwd = 2, col = "gray12")
lines(t, CL.brcc[, 2], lwd = 2, col = "gray12")
lines(t, CL.rcc[, 1],  lwd = 2, lty = 4, col = "#4DAF4A")
lines(t, CL.rcc[, 2],  lwd = 2, lty = 4, col = "#4DAF4A")
abline(h = CL.bcc, col = "#0080ff", lty = 2, lwd = 2)
points(x = twarnings, y = ywarnings, col = "red", pch = 15)
axis(1, seq(min(t)+1, max(t), by = 2), FF(seq(min(t) + 1, max(t), by = 2), 0))
axis(2, seq(R[1], R[2], length.out = 5),FF(seq(R[1], R[2], length.out = 5), 4))
mtext("", side = 1, line = 2.0, cex = 1.8)
mtext("Proporção de massa não convertida", side = 2.1, line = 2.2, cex = 1.8)
mtext("Amostras", side = 1.1, line = 2.2, cex = 1.8)
legend("bottom", c("BRCC", "BCC", "RCC"), lty = c(1, 2, 4), lwd = 2, col = c("gray12", "#0080ff", "#4DAF4A"),
       inset = c(0,1), xpd = TRUE, horiz = TRUE, bty = "n")
graphics.off()

tab <- res.brcc$model
tab <- FF(tab)
tab <- cbind(tab, paste0("(", tab[, 3], ", ", tab[, 4], ")"))
tab <- tab[, c(1, 2, 6, 5)]
rownames(tab) <- c(paste0("$\\beta_", 0:5, "$"), paste0("$\\gamma_", 0:2, "$"))
colnames(tab) <- c("Estimativa", "Erro padr?o", "IC 95\\%", "valor-\\emph{p}")
print(xtable(tab), include.rownames = T, include.colnames = T, sanitize.text.function = force)



# Air relative humidity ---------------------------------------------------
humidity         <- read.table(file = paste0(wd, "/dados/umidade83767.csv"), header = T, sep = ";", stringsAsFactors = F)
humidity$inverno <- with(humidity, ifelse(estacaoano == "Inverno", 1, 0))
humidity$verao   <- with(humidity, ifelse(estacaoano == "Verao", 1, 0))
humidity$outono  <- with(humidity, ifelse(estacaoano == "Outono", 1, 0))
humidity$y       <- humidity$umidade.relativa.media / 100
head(humidity)

## Calculating the control limits
res.brcc    <- brcc(y = humidity$y, X = humidity[, c("inverno", "verao", "outono")], Z = humidity[, c("inverno", "verao", "outono")])
CL.brcc     <- res.brcc$crtlim
CL.bcc      <- bcc(y = humidity$y)
CL.bcc.bias <- bcc(y = humidity$y, bias = TRUE)
CL.rcc      <- rcc(model = y ~ inverno + verao + outono, data = humidity)

## Plotting the control limits and the data
y         <- humidity$y
t         <- 1:length(y)
R         <- range(y, CL.bcc, CL.brcc, CL.rcc)
twarnings <- which(y < CL.brcc[, 1] | y > CL.brcc[, 2])
ywarnings <- y[twarnings]

pdf("humidityCL.pdf", width = 9)
par(mar =  c(3.2, 3.2, 2.0, 0.5), cex = 1.6)
plot(t, y, type = "p", ylim = R, cex = 0.7, pch = 15, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
lines(t, CL.brcc[, 1], lwd = 2, col = "gray12")
lines(t, CL.brcc[, 2], lwd = 2, col = "gray12")
lines(t, CL.rcc[, 1],  lwd = 2, lty = 1, col = "#4DAF4A")
lines(t, CL.rcc[, 2],  lwd = 2, lty = 1, col = "#4DAF4A")
abline(h = CL.bcc, col = "#0080ff", lty = 2, lwd = 2)
points(x = twarnings, y = ywarnings, col = "red", pch = 15, cex = 0.7)
axis(1, seq(min(t), max(t), l = 6), FF(seq(min(t), max(t), l = 6), 0))
axis(2, seq(R[1], R[2], length.out = 5),FF(seq(R[1], R[2], length.out = 5), 2))
mtext("", side = 1, line = 2.0, cex = 1.8)
mtext("Umidade relativa", side = 2.1, line = 2.2, cex = 1.8)
mtext("Amostras", side = 1.1, line = 2.2, cex = 1.8)
legend("bottom", c("BRCC", "BCC", "RCC"), lty = c(1, 2, 4), lwd = 2, col = c("gray12", "#0080ff", "#4DAF4A"),
       inset = c(0,1), xpd = TRUE, horiz = TRUE, bty = "n")
graphics.off()

tab <- res.brcc$model
tab <- FF(tab)
tab <- cbind(tab, paste0("(", tab[, 3], ", ", tab[, 4], ")"))
tab <- tab[, c(1, 2, 6, 5)]
rownames(tab) <- c(paste0("$\\beta_", 0:3, "$"), paste0("$\\gamma_", 0:3, "$"))
colnames(tab) <- c("Estimativa", "Erro padr?o", "IC 95\\%", "valor-\\emph{p}")
print(xtable(tab), include.rownames = T, include.colnames = T, sanitize.text.function = force)


## Identificando observa??es fora de controle
out.of.control   <- humidity[twarnings, c(2, 3, 6)]
nrow(out.of.control)
mean.by.season   <- with(humidity, tapply(y, estacaoano, mean))
median.by.season <- with(humidity, tapply(y, estacaoano, median))
max.by.season    <- with(humidity, tapply(y, estacaoano, max))
min.by.season    <- with(humidity, tapply(y, estacaoano, min))
sd.by.season     <- with(humidity, tapply(y, estacaoano, sd))
tab              <- cbind(mean.by.season, median.by.season, max.by.season, sd.by.season) * 100
cat(paste0("\\multirow{4}{*}{", FF(tab[1, ], 2), "}"), sep = " & ")
cat(paste0("\\multirow{2}{*}{", FF(tab[2, ], 2), "}"), sep = " & ")
cat(FF(tab[3, ], 2), sep = " & ")
cat(FF(tab[4, ], 2), sep = " & ")

xtable(out.of.control[order(out.of.control$estacaoano), ])


