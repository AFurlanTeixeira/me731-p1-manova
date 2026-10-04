# =============================================================================
# p1_manova.R — ME731 (Unicamp, 2s2026) — Projeto 1
# MANOVA a um fator: medidas morfométricas de pinguins (Palmer Archipelago)
#
# Como executar (a partir da RAIZ do repositório):
#     Rscript R/p1_manova.R
# ou, no RStudio, abra o projeto na raiz e use "Source".
#
# Saídas:
#   results/tables/*.csv     tabelas citadas no relatório
#   results/figures/*.png    figuras citadas no relatório
#   results/log_execucao.txt registro completo da execução + sessionInfo()
#
# Dependências: apenas R base (>= 4.0). Nenhum pacote externo.
# Referência-guia: Johnson & Wichern (4a ed.), Seções 4.6, 6.4, 6.5 e 6.9.
# =============================================================================

# ---- 0. Configuração e reprodutibilidade -----------------------------------
.wd_original <- getwd()   # restaurado ao final (evita efeito colateral de source())
local({
  # Garante que o diretório de trabalho é a raiz do repositório, quer o script
  # seja chamado via Rscript, quer via source() no RStudio.
  if (!file.exists("R/funcoes.R")) {
    arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    caminho <- if (length(arg)) sub("^--file=", "", arg[1]) else
      tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(caminho)) setwd(file.path(dirname(normalizePath(caminho)), ".."))
  }
  if (!file.exists("R/funcoes.R"))
    stop("Execute o script a partir da raiz do repositório (pasta que contém R/ e data/).")
})

# Os arquivos estão em UTF-8 (acentos). Em sessões com locale não-UTF-8
# (ex.: LANG=C em Linux, ou R < 4.2 no Windows) o source() falha ao converter
# os acentos. Tentamos ativar um locale UTF-8 antes de prosseguir.
if (!isTRUE(l10n_info()[["UTF-8"]])) {
  for (loc in c("C.UTF-8", "en_US.UTF-8", "pt_BR.UTF-8", "Portuguese_Brazil.utf8"))
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  if (!isTRUE(l10n_info()[["UTF-8"]]))
    warning("Locale não-UTF-8: acentos podem aparecer incorretos nas saídas.")
}
source("R/funcoes.R", encoding = "UTF-8")
stopifnot(exists("salvar_tabela"), exists("sqpc_manova"))  # source() completo

SEMENTE <- 731          # semente única para todo o processo aleatório
ALFA    <- 0.05
R_PERM  <- 9999         # réplicas do teste de permutação
R_SIM   <- 20000        # réplicas do Monte Carlo (EP de MC ~ 0,0015 em alfa = 0,05)

dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
options(digits = 6, scipen = 3, OutDec = ".")

while (sink.number() > 0) sink()            # limpa sinks de execuções anteriores
log_con <- file("results/log_execucao.txt", open = "wt", encoding = "UTF-8")
sink(log_con, split = TRUE)                 # tudo o que é impresso vai p/ o log

secao <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n",
                           strrep("=", 78), "\n", sep = "")

# Paleta categórica (validada para daltonismo; ver README) + símbolos
# distintos por grupo, de modo que a identidade nunca dependa só da cor.
CORES <- c(Adelie = "#2a78d6", Chinstrap = "#eb6834", Gentoo = "#1baf7a")
PCH   <- c(Adelie = 16, Chinstrap = 17, Gentoo = 15)

abrir_png <- function(nome, w = 2000, h = 1500)
  png(file.path("results/figures", nome), width = w, height = h, res = 250)

# ---- 1. Leitura e verificação de integridade dos dados ---------------------
secao("1. LEITURA DOS DADOS")
arq <- "data/raw/penguins.csv"
MD5_ESPERADO <- "a06a0210251465a86fb970018292304d"
md5 <- unname(tools::md5sum(arq))
cat("MD5 do arquivo:", md5, "\n")
if (!identical(md5, MD5_ESPERADO))
  warning("MD5 diferente do esperado: o arquivo de dados foi alterado ",
          "(ou convertido para CRLF pelo Git). Resultados podem divergir.")

dados <- read.csv(arq, stringsAsFactors = FALSE, na.strings = c("NA", ""))
cat("Dimensão bruta:", dim(dados), "\n")
str(dados)

# ---- 2. Pré-processamento ---------------------------------------------------
secao("2. PRÉ-PROCESSAMENTO")
VARS <- c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g")
ROTULOS <- c(bill_length_mm = "Comprimento do bico (mm)",
             bill_depth_mm = "Profundidade do bico (mm)",
             flipper_length_mm = "Comprimento da nadadeira (mm)",
             body_mass_g = "Massa corporal (g)")

faltantes <- colSums(is.na(dados[, c("species", VARS, "sex")]))
print(faltantes)
salvar_tabela(data.frame(n_faltantes = faltantes), "t01_dados_faltantes.csv")

# Casos completos nas 4 respostas (o sexo NÃO entra na análise principal,
# então não descartamos as 9 linhas com sexo faltante).
cc <- complete.cases(dados[, VARS])
df <- dados[cc, c("species", "sex", "island", "year", VARS)]
df$species <- factor(df$species, levels = c("Adelie", "Chinstrap", "Gentoo"))
rownames(df) <- paste0("obs", which(cc))      # rastreia a linha original
X <- as.matrix(df[, VARS])
grupo <- df$species

n <- nrow(X); p <- ncol(X); g <- nlevels(grupo); n_l <- table(grupo)
cat(sprintf("n = %d, p = %d, g = %d  |  n > p: %s ; n_l - 1 >= p em todos: %s\n",
            n, p, g, n > p, all(n_l - 1 >= p)))
print(n_l)
stopifnot(n >= 30, p >= 4, n > p, n - g >= p)  # requisitos do PDD e de W p.d.

# Tabela descritiva por grupo
desc <- do.call(rbind, lapply(levels(grupo), function(l) {
  Xl <- X[grupo == l, , drop = FALSE]
  data.frame(especie = l, variavel = VARS, n = nrow(Xl),
             media = colMeans(Xl), dp = apply(Xl, 2, sd),
             min = apply(Xl, 2, min), max = apply(Xl, 2, max))
}))
rownames(desc) <- NULL
print(desc, digits = 4)
salvar_tabela(desc, "t02_descritiva_por_especie.csv")

# Matrizes de correlação dentro de cada grupo (estrutura de dependência)
for (l in levels(grupo)) {
  cat("\nCorrelação dentro de", l, ":\n")
  R_l <- cor(X[grupo == l, ])
  print(round(R_l, 3))
  salvar_tabela(R_l, sprintf("t03_correlacao_%s.csv", l))
}

# Figura 1: matriz de dispersão colorida por espécie (com símbolos)
abrir_png("f01_matriz_dispersao.png", 2200, 2200)
pairs(X, labels = sub(" \\(", "\n(", ROTULOS[VARS]),
      col = CORES[grupo], pch = PCH[grupo], cex = 0.6, gap = 0.4,
      oma = c(3, 3, 9, 3), main = "Medidas morfométricas por espécie")
# legenda sobre uma camada que cobre a figura inteira
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("top", inset = 0.005, horiz = TRUE, bty = "n",
       legend = levels(grupo), col = CORES, pch = PCH, cex = 0.9)
invisible(dev.off())

# Figura 2: boxplots por variável
abrir_png("f02_boxplots.png", 2400, 1200)
op <- par(mfrow = c(1, 4), mar = c(4, 4.5, 3, 1))
for (v in VARS) {
  boxplot(X[, v] ~ grupo, col = adjustcolor(CORES, 0.35), border = CORES,
          main = ROTULOS[v], xlab = "", ylab = "", las = 1, cex.main = 0.9)
}
par(op); invisible(dev.off())

# ---- 3. Diagnósticos ANTES do ajuste ----------------------------------------
secao("3. DIAGNÓSTICOS PRÉ-AJUSTE")

# 3a. Normalidade univariada por grupo: QQ-plots + Shapiro-Wilk (descritivo)
abrir_png("f03_qqplots_univariados.png", 2400, 1800)
op <- par(mfrow = c(g, p), mar = c(3.5, 3.5, 2.5, 0.8), mgp = c(2.2, 0.7, 0))
sw <- list()
for (l in levels(grupo)) for (v in VARS) {
  x <- X[grupo == l, v]
  qqnorm(x, main = paste(l, "-", sub(" \\(.*", "", ROTULOS[v])),
         col = CORES[l], pch = PCH[l], cex = 0.6, cex.main = 0.85, las = 1)
  qqline(x, col = "grey40", lty = 2)
  sw[[length(sw) + 1]] <- data.frame(especie = l, variavel = v,
                                     W = shapiro.test(x)$statistic,
                                     p_valor = shapiro.test(x)$p.value)
}
par(op); invisible(dev.off())
sw <- do.call(rbind, sw); rownames(sw) <- NULL
print(sw, digits = 4)
salvar_tabela(sw, "t04_shapiro_univariado.csv")

# 3b. Normalidade multivariada por grupo: Mardia + gráfico qui-quadrado
mard <- t(sapply(levels(grupo), function(l) mardia(X[grupo == l, ])))
print(round(mard, 4))
salvar_tabela(mard, "t05_mardia_por_especie.csv")

abrir_png("f04_quiquadrado_por_especie.png", 2400, 900)
op <- par(mfrow = c(1, g), mar = c(4.5, 4.5, 3, 1))
rQ <- numeric(0)
for (l in levels(grupo)) {
  Xl <- X[grupo == l, ]
  d2 <- dist_mahalanobis2(Xl)
  rQ[l] <- grafico_quiquadrado(d2, p, main = paste("Gráfico qui-quadrado:", l),
                               col = CORES[l], pch = PCH[l])
}
par(op); invisible(dev.off())
cat("Correlação no gráfico qui-quadrado por espécie:\n"); print(round(rQ, 4))

# 3c. Observações atípicas: d^2 > chi^2_p(0,001) dentro do próprio grupo
lim_out <- qchisq(1 - 0.001, p)
atip <- do.call(rbind, lapply(levels(grupo), function(l) {
  Xl <- X[grupo == l, ]; d2 <- dist_mahalanobis2(Xl)
  idx <- which(d2 > lim_out)
  if (!length(idx)) return(NULL)
  data.frame(obs = names(d2)[idx], especie = l, d2 = d2[idx], Xl[idx, , drop = FALSE])
}))
cat(sprintf("Limite chi2_%d(0,001) = %.3f; atípicos encontrados: %d\n",
            p, lim_out, if (is.null(atip)) 0 else nrow(atip)))
if (!is.null(atip)) { print(atip); salvar_tabela(atip, "t06_atipicos.csv") }

# 3d. Homogeneidade das matrizes de covariâncias: M de Box + |S_l|
bm <- box_m(X, grupo)
cat(sprintf("Box M = %.3f; u = %.5f; C = (1-u)M = %.3f; gl = %d; p-valor = %.3g\n",
            bm$M, bm$u, bm$C, bm$gl, bm$p_valor))
var_gen <- data.frame(log_det_S = c(bm$log_det_S_l, pooled = bm$log_det_S_pooled))
print(var_gen)
salvar_tabela(data.frame(M = bm$M, u = bm$u, C = bm$C, gl = bm$gl,
                         p_valor = bm$p_valor), "t07_box_m.csv")
salvar_tabela(var_gen, "t08_variancia_generalizada.csv")
# Razão entre desvios-padrão por variável (heurística: > 2 preocupa mais)
razao_dp <- apply(sapply(levels(grupo), function(l) apply(X[grupo == l, ], 2, sd)),
                  1, function(s) max(s) / min(s))
cat("Razão max/min dos desvios-padrão entre grupos:\n"); print(round(razao_dp, 3))

# ---- 4. Ajuste da MANOVA ----------------------------------------------------
secao("4. MANOVA A UM FATOR")
s <- sqpc_manova(X, grupo)
cat("Matriz B (entre):\n");  print(round(s$B, 2))
cat("Matriz W (dentro):\n"); print(round(s$W, 2))
stopifnot(isTRUE(all.equal(s$B + s$W, s$T, check.attributes = FALSE)))  # (6-36)
salvar_tabela(s$B, "t09_matriz_B.csv"); salvar_tabela(s$W, "t10_matriz_W.csv")

est <- estatisticas_manova(s)
stopifnot(isTRUE(all.equal(est$wilks, est$wilks_det)))  # (6-38) = prod 1/(1+l)
tw <- teste_wilks(est$wilks, n, p, g)
cat(sprintf("Autovalores de W^-1 B: %s\n", paste(round(est$autovalores, 4), collapse = ", ")))
cat(sprintf("Lambda de Wilks = %.6f\n", est$wilks))
print(tw)

# Conferência com a implementação de referência do R (stats::manova)
fit <- manova(X ~ grupo)
smf <- lapply(c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"),
              function(t) summary(fit, test = t)$stats[1, ])
tab_testes <- do.call(rbind, smf)
rownames(tab_testes) <- c("Wilks", "Pillai", "Hotelling-Lawley", "Roy")
print(tab_testes)
cat("Obs.: para Roy, o F reportado é um LIMITE SUPERIOR, logo o p-valor é um limite inferior.\n")
cat("Obs.: p-valores < 1e-16 devem ser reportados como 'p < 0,001'; abaixo da\n",
    "     precisão da máquina (ex.: 4.94e-324 = underflow) não têm significado.\n")
# (tolerância RELATIVA; obs.: o R reporta para Roy o maior autovalor
#  lambda_1 de W^{-1}B, e não theta = lambda_1/(1+lambda_1))
confere <- function(a, b) isTRUE(all.equal(unname(a), unname(b), tolerance = 1e-8))
stopifnot(confere(tab_testes["Wilks", 2], est$wilks),
          confere(tab_testes["Pillai", 2], est$pillai),
          confere(tab_testes["Hotelling-Lawley", 2], est$hotelling),
          confere(tab_testes["Roy", 2], est$autovalores[1]))
# Para g = 3 a aproximação F de Rao (usada pelo R) é EXATA e deve coincidir
# com a Tabela 6.3 de J&W:
stopifnot(confere(tab_testes["Wilks", "approx F"], tw$exato["F"]))
salvar_tabela(tab_testes, "t11_testes_manova.csv")

eta2_mult <- 1 - est$wilks^(1 / min(p, g - 1))
cat(sprintf("eta^2 multivariado = 1 - Lambda^(1/s) = %.4f\n", eta2_mult))
salvar_tabela(data.frame(
  estatistica = c("Wilks", "F_exato (Tab. 6.3)", "gl1", "gl2", "p_exato",
                  "Bartlett qui2 (6-39)", "gl", "p_Bartlett", "eta2_mult"),
  valor = c(est$wilks, tw$exato, tw$bartlett, eta2_mult)),
  "t12_wilks_exato_bartlett.csv")

# Propriedade de invariância: Lambda não muda sob X -> XA + 1c' (A não
# singular). Ex.: trocar gramas por quilogramas e milímetros por centímetros.
A <- diag(c(0.1, 0.1, 0.1, 0.001))
s_inv <- sqpc_manova(X %*% A + 5, grupo)
L_inv <- estatisticas_manova(s_inv)$wilks
cat(sprintf("Invariância: Lambda original = %.10f; transformado = %.10f\n",
            est$wilks, L_inv))
stopifnot(confere(L_inv, est$wilks))

# ---- 5. Análises de acompanhamento (pós-teste global) ----------------------
secao("5. COMPARAÇÕES SIMULTÂNEAS (BONFERRONI, RESULT 6.5)")
icb <- ic_bonferroni(s, alpha = ALFA)
cat(sprintf("m = %d afirmações; t crítico = t_{%d}(%.5f) = %.4f\n",
            attr(icb, "m"), n - g, ALFA / (2 * attr(icb, "m")), attr(icb, "t_critico")))
print(icb, digits = 4)
salvar_tabela(icb, "t13_ic_bonferroni.csv")

# Figura 5: intervalos de Bonferroni padronizados pelo d.p. combinado,
# para que variáveis com unidades diferentes caibam num mesmo eixo.
sp <- sqrt(diag(s$W) / (n - g))
icb$dif_pad <- icb$diferenca / sp[icb$variavel]
icb$li_pad  <- icb$li / sp[icb$variavel]
icb$ls_pad  <- icb$ls / sp[icb$variavel]
abrir_png("f05_ic_bonferroni.png", 2200, 1500)
op <- par(mar = c(4.5, 15, 3, 1))
y <- rev(seq_len(nrow(icb)))
plot(icb$dif_pad, y, xlim = range(c(icb$li_pad, icb$ls_pad, 0)), yaxt = "n",
     pch = 19, xlab = "Diferença de médias / d.p. combinado (unidades de DP)",
     ylab = "", main = sprintf("ICs simultâneos de Bonferroni (%d%%, m = %d)",
                               100 * (1 - ALFA), attr(icb, "m")), las = 1)
segments(icb$li_pad, y, icb$ls_pad, y, lwd = 2)
abline(v = 0, lty = 2, col = "grey40")
axis(2, at = y, labels = paste(icb$contraste, "|",
                               sub(" \\(.*", "", ROTULOS[icb$variavel])),
     las = 1, cex.axis = 0.7)
par(op); invisible(dev.off())

# Figura 6: variáveis canônicas (descritivo; J&W Seção 11.7)
vc <- variaveis_canonicas(X, s)
cat("Coeficientes canônicos (normalizados por S_pooled):\n"); print(round(vc$coef, 4))
cat("Proporção da separação por eixo:", round(vc$prop, 4), "\n")
# Coeficientes padronizados (multiplicados pelo d.p. combinado de cada
# variável): comparáveis entre variáveis com unidades distintas.
coef_pad <- vc$coef * sqrt(diag(s$W) / (n - g))
cat("Coeficientes canônicos padronizados:\n"); print(round(coef_pad, 4))
salvar_tabela(cbind(vc$coef, setNames(as.data.frame(coef_pad),
                                      paste0(colnames(coef_pad), "_pad"))),
              "t14_coef_canonicos.csv")
abrir_png("f06_variaveis_canonicas.png", 1800, 1500)
plot(vc$escores, col = CORES[grupo], pch = PCH[grupo], cex = 0.7, las = 1,
     xlab = sprintf("CAN1 (%.1f%% da separação)", 100 * vc$prop[1]),
     ylab = sprintf("CAN2 (%.1f%% da separação)", 100 * vc$prop[2]),
     main = "Escores nas variáveis canônicas discriminantes")
cen <- apply(vc$escores, 2, function(z) tapply(z, grupo, mean))
points(cen, pch = 4, cex = 2, lwd = 3)
text(cen, labels = rownames(cen), pos = 3, offset = 1.1, font = 2, cex = 0.9)
legend("topright", legend = levels(grupo), col = CORES, pch = PCH, bty = "n")
invisible(dev.off())

# ---- 6. Diagnósticos DEPOIS do ajuste (resíduos) ---------------------------
secao("6. DIAGNÓSTICOS PÓS-AJUSTE (RESÍDUOS)")
E <- X - s$medias[as.character(grupo), ]        # e_lj = x_lj - xbar_l
# Os resíduos são padronizados por S_pooled, ou seja, SOB o modelo ajustado
# (Sigma comum). Como as Sigma_l diferem (Box M), parte da curtose dos
# resíduos reflete a heterogeneidade, e não apenas a não-normalidade.
Sp <- s$W / (n - g)
d2_res <- dist_mahalanobis2(E, centro = rep(0, p), S = Sp)
abrir_png("f07_residuos_quiquadrado.png", 2400, 1100)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
rQ_res <- grafico_quiquadrado(d2_res, p, main = "Resíduos: gráfico qui-quadrado",
                              col = CORES[grupo], pch = PCH[grupo])
plot(seq_along(d2_res), d2_res, col = CORES[grupo], pch = PCH[grupo], cex = 0.7,
     xlab = "Índice da observação (ordem do arquivo)", ylab = expression(d^2),
     main = "Distâncias dos resíduos por observação", las = 1)
abline(h = lim_out, lty = 2, col = "grey40")
par(op); invisible(dev.off())
cat(sprintf("Correlação no gráfico qui-quadrado dos resíduos: %.4f\n", rQ_res))
mard_res <- mardia(E)
print(round(mard_res, 4))
salvar_tabela(t(mard_res), "t15_mardia_residuos.csv")

# ---- 7. Robustez e sensibilidade -------------------------------------------
secao("7. ROBUSTEZ E SENSIBILIDADE")

# 7a. Teste de permutação para Wilks (livre de normalidade)
perm <- permutacao_wilks(X, grupo, R = R_PERM, semente = SEMENTE)
cat(sprintf("Permutação (R = %d, semente = %d): Lambda_obs = %.5f; min(Lambda_perm) = %.5f; p = %.5f\n",
            perm$R, perm$semente, perm$L_obs, min(perm$L_perm), perm$p_valor))
abrir_png("f08_permutacao.png", 1800, 1200)
hist(perm$L_perm, breaks = 50, col = "grey80", border = "white", las = 1,
     xlim = range(c(perm$L_perm, perm$L_obs)),
     main = sprintf("Distribuição de permutação de Lambda (R = %d)", R_PERM),
     xlab = expression(Lambda^"*"), ylab = "Frequência")
abline(v = perm$L_obs, lwd = 2, col = "#2a78d6")
text(perm$L_obs, par("usr")[4] * 0.9, "observado", pos = 4, col = "grey20")
invisible(dev.off())

# 7b. Sensibilidade: MANOVA com transformação log (estabiliza variâncias)
s_log <- sqpc_manova(log(X), grupo)
L_log <- estatisticas_manova(s_log)$wilks
p_log <- teste_wilks(L_log, n, p, g)$exato[["p_valor"]]
bm_log <- box_m(log(X), grupo)
cat(sprintf("Log-escala: Lambda = %.5f; Box C = %.2f (p = %.3g)\n",
            L_log, bm_log$C, bm_log$p_valor))

# 7c. Sensibilidade: excluindo atípicos detectados em 3c
if (!is.null(atip)) {
  keep <- !(rownames(X) %in% atip$obs)
  s_sem <- sqpc_manova(X[keep, ], grupo[keep])
  L_sem <- estatisticas_manova(s_sem)$wilks
  p_sem <- teste_wilks(L_sem, sum(keep), p, g)$exato[["p_valor"]]
  cat(sprintf("Sem %d atípicos: Lambda = %.5f\n", sum(!keep), L_sem))
} else { L_sem <- NA; p_sem <- NA }

# 7d. Fonte de heterogeneidade: sexo. MANOVA a dois fatores (J&W Seção 6.6)
#     com espécie, sexo e interação, nos n = 333 casos com sexo conhecido.
df2 <- df[!is.na(df$sex), ]
X2 <- as.matrix(df2[, VARS])
fit2 <- manova(X2 ~ species * sex, data = df2)
cat("\nMANOVA dois fatores (espécie x sexo), Wilks — somas de quadrados sequenciais (tipo I):\n")
print(summary(fit2, test = "Wilks"))
tab2 <- summary(fit2, test = "Wilks")$stats
salvar_tabela(tab2, "t16_manova_especie_sexo.csv")
bm_sex <- sapply(split(seq_len(nrow(df2)), df2$sex), function(i)
  box_m(X2[i, ], df2$species[i])$p_valor)
cat("Box M por sexo (p-valores):\n"); print(signif(bm_sex, 3))
mard_cel <- t(sapply(split(seq_len(nrow(df2)), interaction(df2$species, df2$sex)),
                     function(i) mardia(X2[i, ])[c("n", "p_assimetria", "p_curtose")]))
cat("Mardia por célula espécie x sexo:\n"); print(round(mard_cel, 4))
salvar_tabela(mard_cel, "t17_mardia_por_celula.csv")

salvar_tabela(data.frame(
  analise = c("Principal", "Permutação", "Log", "Sem atípicos"),
  wilks = c(est$wilks, perm$L_obs, L_log, L_sem),
  p_valor = c(tw$exato[["p_valor"]], perm$p_valor, p_log, p_sem)),
  "t18_sensibilidade.csv")

# ---- 8. Monte Carlo: tamanho dos testes sob H0 ------------------------------
# Tamanho empírico do teste de Wilks com os MESMOS n_l, p e g dos dados, em
# quatro cenários, todos sob H0 (médias iguais):
#   (a) normal, Sigma comum (= S_pooled), F exato da Tabela 6.3
#   (b) idem, aproximação de Bartlett (6-39)
#   (c) normal, Sigma_l = S_l observadas (heterogeneidade real)
#   (d) bootstrap dos resíduos de cada grupo: forma (não-normalidade) E
#       heterogeneidade reais dos dados, médias igualadas por construção
secao("8. ESTUDO DE MONTE CARLO (TAMANHO EMPÍRICO)")
set.seed(SEMENTE)
L_chol <- chol(Sp)                      # chol() devolve U com U'U = Sp
gr_sim <- rep(factor(levels(grupo), levels = levels(grupo)), times = n_l)
Sl_chol <- lapply(bm$S_l, chol)
res_l <- split.data.frame(X - s$medias[as.character(grupo), ], grupo)
rejeita <- function(Xs) teste_wilks(estatisticas_manova(sqpc_manova(Xs, gr_sim))$wilks,
                                    n, p, g)
sim_rej <- t(replicate(R_SIM, {
  X0 <- matrix(rnorm(n * p), n, p) %*% L_chol
  X1 <- do.call(rbind, lapply(levels(grupo), function(l)
    matrix(rnorm(n_l[l] * p), n_l[l], p) %*% Sl_chol[[l]]))
  X2b <- do.call(rbind, lapply(levels(grupo), function(l)
    res_l[[l]][sample.int(n_l[l], n_l[l], replace = TRUE), , drop = FALSE]))
  r0 <- rejeita(X0)
  c(a_exato_sigma_comum    = r0$exato[["p_valor"]] < ALFA,
    b_bartlett_sigma_comum = r0$bartlett[["p_valor"]] < ALFA,
    c_exato_sigma_heterog  = rejeita(X1)$exato[["p_valor"]] < ALFA,
    d_exato_bootstrap_res  = rejeita(X2b)$exato[["p_valor"]] < ALFA)
}))
tam <- colMeans(sim_rej)
ep_mc <- sqrt(tam * (1 - tam) / R_SIM)
tab_mc <- data.frame(tamanho_empirico = tam, ic95_li = tam - 1.96 * ep_mc,
                     ic95_ls = tam + 1.96 * ep_mc, R = R_SIM)
print(round(tab_mc, 4))
cat(sprintf("Nível nominal = %.2f\n", ALFA))  # sem tempo de execução: log determinístico
salvar_tabela(tab_mc, "t19_monte_carlo_tamanho.csv")

# ---- 9. Registro do ambiente ------------------------------------------------
secao("9. AMBIENTE")
print(sessionInfo())
cat("\nExecução concluída sem erros.\n")
sink()
close(log_con)
setwd(.wd_original)
