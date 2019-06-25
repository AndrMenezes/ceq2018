rm(list = ls())

source("my_plot_qc.R")

pkgs <- c("ggplot2", "qcc", "xtable", "goftest", "tidyr")
sapply(pkgs, require, character.only = T)

FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}

# Dados da fase 1 ---------------------------------------------------------
x1          <- c(2020, 2000, 2025, 2010, 2005, 2010, 2010, 2010, 2005, 2000, 2000, 2000, 2000, 1980, 2000, 2000, 2000, 2000, 1998, 2000, 2000, 1998, 2002, 2000, 2002, 2000, 2002, 2002, 2005, 2005, 2000, 2003, 2000, 2002, 1996, 2000, 2004, 2002, 2000, 2000, 2005, 2000, 2000, 2005, 2000, 1998, 2000, 1998, 2000, 2010, 2000, 2002, 2000, 2025, 2025, 2000, 2010, 2005, 2010, 2020, 2005, 2000, 2000, 2000, 2005, 2010, 2010, 2000, 2005, 2005, 2000, 2000, 1998, 2000, 2000, 1995, 2000, 2000, 1993, 1998, 2000, 2005, 2005, 2000, 2000, 2000, 2005, 2010, 2010, 2005, 2005, 2000, 2000, 2005, 2000, 2000, 2000, 1998, 2000, 1995)
dados.fase1 <- matrix(x1, ncol = 5, byrow = T)
dados.fase1

## Estimando os par?metros do processo e contruindo os gr?ficos
aux1 <- my.plot.qc(dados = dados.fase1)
ggsave(filename = "fase1-media.pdf", plot = aux1$plots$media,     device = "pdf", height = 5, width = 9)
ggsave(filename = "fase1-R.pdf",     plot = aux1$plots$amplitude, device = "pdf", height = 5, width = 9)

auxboot <- my.plot.qc(dados = dados.fase1, cimedia = "bootstrap", ciR = "bootstrap")
ggsave(filename = "fase1-media-boot.pdf", plot = auxboot$plots$media,     device = "pdf", height = 5, width = 9)
ggsave(filename = "fase1-R-boot.pdf",     plot = auxboot$plots$amplitude, device = "pdf", height = 5, width = 9)

## Retirando as obs 1, 3, 11, 12 e 16, pois foram influenciadas por causas especiais
dados.fase1 <- dados.fase1[-c(1, 3, 11, 12, 16),]
aux2        <- my.plot.qc(dados = dados.fase1)
aux2$plots$amplitude <- aux2$plots$amplitude + scale_color_manual(values = c("black", "red"))
ggsave(filename = "fase1-media2.pdf", plot = aux2$plots$media,     device = "pdf", height = 5, width = 9)
ggsave(filename = "fase1-R2.pdf",     plot = aux2$plots$amplitude, device = "pdf", height = 5, width = 9)


# Dados da fase 2 ---------------------------------------------------------
x2          <- c(2000, 2000, 1980, 1980, 1995, 1990, 2005, 2000, 1990, 1995, 2010, 2010, 2000, 2020, 2000, 2000, 2005, 2000, 2010, 1990, 2000, 2005, 2005, 2005, 2005, 2020, 2010, 2010, 2010, 2020, 2005, 2020, 2020, 2020, 2010, 2000, 2005, 2005, 2000, 2010, 2000, 2003, 2010, 2010, 2010, 1990, 2000, 2010, 2010, 2000, 2005, 2005, 2000, 2020, 2010, 2010, 2000, 2005, 2010, 2025, 2005, 2000, 2010, 2005, 2010, 2010, 2010, 2020, 2010, 2000, 2005, 2000, 2010, 2005, 2010)
dados.fase2 <- matrix(x2, ncol = 5, byrow = T)
todos       <- rbind(dados.fase1, dados.fase2)

aux3 <- my.plot.qc(dados = todos, parms = aux2$parms)
p1 <- aux3$plots$media + geom_vline(xintercept = 15.5, size = 0.8, col = "blue")
p2 <- aux3$plots$amplitude + geom_vline(xintercept = 15.5, size = 0.8, col = "blue")
ggsave(filename = "fase2-media.pdf", plot = p1, device = "pdf", height = 6, width = 12)
ggsave(filename = "fase2-R.pdf",     plot = p2, device = "pdf", height = 6, width = 12)



# Normalidade das médias das m amostras -----------------------------------
testes <- function(x, p, ...)
{
  xsw       <- shapiro.test(x)
  xks       <- ks.test(x  = x, p, ...)
  xcvm      <- cvm.test(x = x, null = p, ...)
  xad       <- ad.test(x  = x, null = p, ...)
  statistic <- c(xsw$statistic, xks$statistic, xcvm$statistic, xad$statistic)
  pvalue    <- c(xsw$p.value, xks$p.value, xcvm$p.value, xad$p.value)
  paste(FF(statistic, 4), " (", FF(pvalue, 4), ")", sep = "")
}
myppplot <- function(x, mle, nome, j)
{
  n       <- length(x)
  sx      <- sort(x)
  emp     <- ppoints(n)  
  teo     <- pnorm(q = sx, mle[1], mle[2])
  pdf(paste0("ppplot-", nome, ".pdf"))
  par(mar = c(3.2, 3.2, 0.5, 0.5), cex = 1.8)
  plot(x = teo, y = emp,  type = 'p', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', xlim = c(0, 1), ylim = c(0,1), col = 1, bty = 'l');box();abline(0, 1, lwd = 1);
  abline(h=seq(0, 1, by = 0.25), v=seq(0, 1, by = 0.25), col = "lightgray", lty = "dotted")
  axis(1, seq(0, 1, by = 0.25),FF(seq(0, 1, by = 0.25), 2))
  axis(2, seq(0, 1, by = 0.25),FF(seq(0, 1, by = 0.25), 2))
  mtext("Probabilidades teóricas", side = 1, line = 2.0, cex = 1.8)
  mtext("Probabilidades empíricas", side = 2.1, line =2.2, cex = 1.8)
  text(x = 0.5, y = 0.98, labels = paste0("Fase ", j), cex = 1.2)
  graphics.off()  
}

medias.fase1 <- rowMeans(dados.fase1)
medias.fase2 <- rowMeans(dados.fase2)

t1 <- testes(medias.fase1, "pnorm", mean = mean(medias.fase1), sd = sd(medias.fase1))
t2 <- testes(medias.fase2, "pnorm", mean = mean(medias.fase2), sd = sd(medias.fase2))
tab <- cbind(t1, t2)
rownames(tab) <- c("SW", "KS", "CvM", "AD")
colnames(tab) <- c("Fase I", "Fase II")
xtable(tab)

myppplot(x = medias.fase1, mle = c(mean(medias.fase1), sd(medias.fase1)), nome = "fase1", j = 1)
myppplot(x = medias.fase2, mle = c(mean(medias.fase2), sd(medias.fase2)), nome = "fase2", j = 2)



# Curva característica de operação ----------------------------------------
myqcc  <- qcc(data = dados.fase1, type = "xbar")
myqcc2 <- qcc(dados.fase1, type="xbar", newdata=dados.fase2)
mybeta <- oc.curves(myqcc)

enes  <- c(1, 5, 10, 15, 20)
k     <- seq(0, 5, by = 0.05)  
L     <- 3
betas <- matrix(ncol = length(enes), nrow = length(k))

for(i in 1:length(enes))
{
  n          <- enes[i]
  betas[, i] <- pnorm(L - k * sqrt(n)) - pnorm(-L - k * sqrt(n))
}

rownames(betas) <- k
colnames(betas) <- paste0("n=", enes)
betas["1.5", ]

df    <- cbind(data.frame(betas), k)
df    <- df %>% gather(key = "n", "beta", -k)
df$n1 <- as.numeric(gsub("n.", "", df$n, fixed = T))
df$n  <- factor(df$n)
levels(df$n) <- paste0("n = ", enes[c(1, 3, 4, 5, 2)])
df$n <- with(df, reorder(n, n1, mean))

df %>% ggplot(aes(x = k, y = beta, col = n)) +
  geom_line(size = 1.1) +
  labs(x = "k", y = expression(beta), col = "") +
  theme_bw() +
  theme(text = element_text(size = 16))
ggsave("oc-curve.pdf", device = "pdf", width = 11, height = 6)

