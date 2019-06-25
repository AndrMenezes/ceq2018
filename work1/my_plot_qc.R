my.plot.qc <- function(dados, parms = NULL, cimedia = "3sigma", ciR = "normal")
{
  # Constantes
  A2s <- c(1.88, 1.023, 0.729, 0.577, 0.483, 0.419, 0.373, 0.337, 0.308, 0.285, 0.266, 0.249, 0.235, 0.223, 0.212, 0.203, 0.194, 0.187, 0.18, 0.173, 0.167, 0.162, 0.157, 0.153)
  d2s <- c(1.128, 1.693, 2.059, 2.326, 2.534, 2.704, 2.847, 2.97, 3.078, 3.173, 3.258, 3.336, 3.407, 3.472, 3.532, 3.588, 3.64, 3.689, 3.735, 3.778, 3.819, 3.858, 3.895, 3.931)
  d3s <- c(0.853, 0.888, 0.88, 0.864, 0.848, 0.833, 0.82, 0.808, 0.797, 0.787, 0.778, 0.77, 0.763, 0.756, 0.75, 0.744, 0.739, 0.734, 0.729, 0.724, 0.72, 0.716, 0.712, 0.708)
  D3s <- c(0, 0, 0, 0, 0, 0.076, 0.136, 0.184, 0.223, 0.256, 0.283, 0.307, 0.328, 0.347, 0.363, 0.378, 0.391, 0.403, 0.415, 0.425, 0.434, 0.443, 0.451, 0.459)
  D4s <- c(3.267, 2.574, 2.282, 2.114, 2.004, 1.924, 1.864, 1.816, 1.777, 1.744, 1.717, 1.693, 1.672, 1.653, 1.637, 1.622, 1.608, 1.597, 1.585, 1.575, 1.566, 1.557, 1.548, 1.541)
  
  # Valores amostrais
  n          <- ncol(dados)
  medias     <- rowMeans(dados)
  Rs         <- apply(dados, 1, function(u) diff(range(u)))
  
  # Parâmetros
  if(is.null(parms))
  {
    media.geral <- mean(medias)
    R.media     <- mean(Rs)
    
    # Constantes
    d2          <- d2s[n - 1]
    d3          <- d3s[n - 1]
    A2          <- A2s[n - 1]
    D3          <- D3s[n - 1]
    D4          <- D4s[n - 1]
    
    # Erro padrão
    Sd0         <- R.media / d2
    Sm1         <- A2 * R.media
    Sm2         <- Sd0 / sqrt(n)
    
    # Limites de confiança média 
    if(cimedia == "3sigma")    CI.media <- media.geral + c(-1, 1) * 3 * Sm2
    if(cimedia == "amplitude") CI.media <- media.geral + c(-1, 1) * Sm1
    if(cimedia == "bootstrap")
    {
      set.seed(1212)
      media.boot  <- colMeans(replicate(10000, sample(medias, replace = T)))
      CI.media    <- quantile(media.boot, probs = c(0.001, 0.999))
    }
    
    # Limites de confiança amplitude (R) 
    if(ciR == "normal") CI.R <- c(D3 * R.media, D4 * R.media)
    if(ciR == "bootstrap")
    {
      set.seed(1212)
      R.boot <- apply(replicate(10000, sample(Rs, replace = T)), 2, function(u) diff(range(u)))
      CI.R   <- quantile(R.boot, probs = c(0.001, 0.999))
    }
  }
  else 
  {
    media.geral <- parms[1, 1]
    CI.media    <- parms[-1, 1]
    R.media     <- parms[1, 2]
    CI.R        <- parms[-1, 2]
  }
  
  # Gráficos
  df  <- data.frame(Amostra = 1:length(medias), Medias = medias, Amplitudes = Rs,
                    idmedia = (medias >= CI.media[1]) & (medias <= CI.media[2]),
                    idR = (Rs >= CI.R[1]) & (Rs <= CI.R[2]))
  
  ## Média
  pmedia <- ggplot(data = df, aes(x = Amostra, y = Medias)) +
    geom_line(linetype = 2, size = 0.4) +
    geom_point(aes(col = idmedia), size = 2) + 
    geom_hline(yintercept = CI.media, linetype = 2, size = 0.4) +
    geom_hline(yintercept = media.geral, size = 0.4) +
    guides(col = FALSE) +
    scale_color_manual(values = c("red", "black")) +
    scale_x_continuous(breaks = df$Amostra) +
    scale_y_continuous(breaks = seq(min(df$Medias), max(df$Medias), l = 5), 
                       labels = round(seq(min(df$Medias), max(df$Medias), l = 5)),
                       sec.axis = dup_axis(name = "", breaks = c(media.geral, CI.media),
                                           labels = c("CL", "LCL", "UCL"))) +  # paste(c("CL", "LCL", "UCL"), " = ", FF(c(media.geral, CI.media2), 3))
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(), text = element_text(size = 16),
          axis.title.y.right = element_blank()) +
    labs(y = "Média amostral", x = "Amostra") 
  
  ## Amplitude
  pamplitude <- ggplot(data = df, aes(x = Amostra, y = Amplitudes)) +
    geom_line(linetype = 2, size = 0.4) +
    geom_point(aes(col = idR), size = 2) + 
    geom_hline(yintercept = CI.R, linetype = 2, size = 0.4) +
    geom_hline(yintercept = R.media, size = 0.4) +
    guides(col = FALSE) +
    scale_color_manual(values = c("red", "black")) +
    scale_x_continuous(breaks = df$Amostra) +
    scale_y_continuous(breaks = seq(min(df$Amplitudes), max(df$Amplitudes), l = 5), 
                       labels = round(seq(min(df$Amplitudes), max(df$Amplitudes), l = 5)),
                       sec.axis = dup_axis(name = "", breaks = c(R.media, CI.R),
                                           labels = c("CL", "LCL", "UCL"))) +  # paste(c("CL", "LCL", "UCL"), " = ", FF(c(media.geral, CI.media2), 3))
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(), text = element_text(size = 16),
          axis.title.y.right = element_blank()) +
    labs(y = "Amplitude amostral", x = "Amostra") 
  
  parms           <- cbind(c(media.geral, CI.media), c(R.media, CI.R))
  colnames(parms) <- c("Media", "R")
  rownames(parms) <- c("CL", "LCL", "UCL")
  plots           <- list(media = pmedia, amplitude = pamplitude)
  return(list(parms = parms, plots = plots))
}