# =============================================================================
# p1_manova_completo.R — ME731 (Unicamp, 2s2026) — Projeto 1
# MANOVA a um fator: medidas morfométricas de pinguins (Palmer Archipelago)
# Discente: Arthur Agostinho Furlan Teixeira (RA 164363)
#
# Arquivo único com todo o código do projeto: as funções auxiliares (Parte A)
# e a análise completa (Parte B). Reproduz todas as tabelas e figuras do
# relatório. Repositório: https://github.com/AFurlanTeixeira/me731-p1-manova
#
# Como executar:
#     Rscript p1_manova_completo.R
# ou, no RStudio, abra este arquivo e use "Source".
#
# Dados: o script procura penguins.csv na pasta deste arquivo (ou em
# data/raw/ do repositório). Se não encontrar, baixa a cópia original do
# pacote palmerpenguins e confere o MD5 antes de usar.
#
# Saídas, gravadas na pasta deste arquivo:
#   results/tables/*.csv     tabelas citadas no relatório (t01 a t22)
#   results/figures/*.png    figuras citadas no relatório (f01 a f08)
#   results/log_execucao.txt registro completo da execução + sessionInfo()
#
# Dependências: apenas R base (>= 4.0). Nenhum pacote externo.
# Tempo de execução: cerca de dois minutos e meio.
# Referência-guia: Johnson & Wichern (4a ed.), Seções 4.6, 6.4, 6.5 e 6.9.
# =============================================================================

# ---- 0. Configuração e reprodutibilidade -----------------------------------
.wd_original <- getwd()   # restaurado ao final (evita efeito colateral de source())
local({
  # Usa como pasta de trabalho a pasta deste arquivo, quer o script seja
  # chamado via Rscript, quer via source() no RStudio.
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  caminho <- if (length(arg)) sub("^--file=", "", arg[1]) else
    tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(caminho)) setwd(dirname(normalizePath(caminho)))
})

# O arquivo está em UTF-8 (acentos). Em sessões com locale não-UTF-8
# (ex.: LANG=C em Linux) tentamos ativar um locale UTF-8 antes de prosseguir.
if (!isTRUE(l10n_info()[["UTF-8"]])) {
  for (loc in c("C.UTF-8", "en_US.UTF-8", "pt_BR.UTF-8", "Portuguese_Brazil.utf8"))
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  if (!isTRUE(l10n_info()[["UTF-8"]]))
    warning("Locale não-UTF-8: acentos podem aparecer incorretos nas saídas.")
}

# #############################################################################
# PARTE A — FUNÇÕES AUXILIARES
# #############################################################################
# -----------------------------------------------------------------------------
# Matrizes de somas de quadrados e produtos cruzados (SQPC)
# [JW] Seção 6.4, eqs. (6-36) e (6-37):  T = B + W
#   B = sum_l n_l (xbar_l - xbar)(xbar_l - xbar)'      (entre grupos, g-1 g.l.)
#   W = sum_l sum_j (x_lj - xbar_l)(x_lj - xbar_l)'    (dentro,      n-g g.l.)
# -----------------------------------------------------------------------------
sqpc_manova <- function(X, grupo) {
  X <- as.matrix(X)
  grupo <- droplevels(as.factor(grupo))
  xbar <- colMeans(X)
  niveis <- levels(grupo)
  p <- ncol(X)
  B <- matrix(0, p, p, dimnames = list(colnames(X), colnames(X)))
  W <- B
  medias <- matrix(NA_real_, length(niveis), p,
                   dimnames = list(niveis, colnames(X)))
  n_l <- integer(length(niveis)); names(n_l) <- niveis
  for (k in seq_along(niveis)) {
    Xl <- X[grupo == niveis[k], , drop = FALSE]
    n_l[k] <- nrow(Xl)
    m_l <- colMeans(Xl)
    medias[k, ] <- m_l
    d <- m_l - xbar
    B <- B + n_l[k] * tcrossprod(d)
    Xc <- sweep(Xl, 2, m_l)
    W <- W + crossprod(Xc)
  }
  Tt <- crossprod(sweep(X, 2, xbar))
  list(B = B, W = W, T = Tt, medias = medias, xbar = xbar, n_l = n_l,
       n = nrow(X), p = p, g = length(niveis))
}

# -----------------------------------------------------------------------------
# Estatísticas de teste da MANOVA a partir dos autovalores de W^{-1}B
# [JW] nota de rodapé 2 da p. 322 e Seção 6.9 (p. 357)
#   Wilks:  Lambda* = |W| / |B + W| = prod 1/(1+lambda_i)        (6-38)
#   Pillai: V = tr[B (B+W)^{-1}]      = sum lambda_i/(1+lambda_i)
#   Lawley-Hotelling: U = tr[B W^{-1}] = sum lambda_i
#   Roy:    theta = max autovalor de B (B+W)^{-1} = lambda_1/(1+lambda_1)
#   (Obs.: a 4a ed., p. 357, imprime "W (B+W)^{-1}" para Roy; o correto é
#    B (B+W)^{-1}. Ver docs/02_derivacoes_referencias.md.)
# -----------------------------------------------------------------------------
estatisticas_manova <- function(s) {
  # autovalores de W^{-1}B via problema simétrico equivalente
  # W^{-1/2} B W^{-1/2}, numericamente estável (W é p.d. quando n - g >= p)
  eW <- eigen(s$W, symmetric = TRUE)
  W_mhalf <- eW$vectors %*% diag(1 / sqrt(eW$values)) %*% t(eW$vectors)
  M <- W_mhalf %*% s$B %*% W_mhalf
  lambda <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
  r <- min(s$p, s$g - 1)              # posto de B
  lambda <- pmax(lambda[seq_len(r)], 0)
  list(
    autovalores = lambda,
    wilks = prod(1 / (1 + lambda)),
    wilks_det = det(s$W) / det(s$B + s$W),   # conferência direta de (6-38)
    pillai = sum(lambda / (1 + lambda)),
    hotelling = sum(lambda),
    roy = lambda[1] / (1 + lambda[1])
  )
}

# -----------------------------------------------------------------------------
# Distribuição EXATA de Wilks para g = 3  ([JW] Tabela 6.3, p. 323):
#   ((n - p - 2)/p) * (1 - sqrt(L))/sqrt(L)  ~  F(2p, 2(n - p - 2))
# Aproximação de Bartlett ([JW] eqs. (6-39)-(6-40)):
#   -(n - 1 - (p + g)/2) ln L  ~  chi^2_{p(g-1)}   (aprox., n grande)
# -----------------------------------------------------------------------------
teste_wilks <- function(L, n, p, g) {
  out <- list()
  if (g == 3) {
    F_ex <- ((n - p - 2) / p) * (1 - sqrt(L)) / sqrt(L)
    df1 <- 2 * p; df2 <- 2 * (n - p - 2)
    out$exato <- c(F = F_ex, df1 = df1, df2 = df2,
                   p_valor = pf(F_ex, df1, df2, lower.tail = FALSE))
  }
  chi <- -(n - 1 - (p + g) / 2) * log(L)
  gl <- p * (g - 1)
  out$bartlett <- c(qui2 = chi, gl = gl,
                    p_valor = pchisq(chi, gl, lower.tail = FALSE))
  out
}

# -----------------------------------------------------------------------------
# Teste M de Box para igualdade de matrizes de covariâncias
# [JW6] Seção 6.6 ("Testing for Equality of Covariance Matrices");
# [MKB] Cap. 5 (TRV para igualdade de covariâncias; correção de Box).
#   M = (n - g) ln|S_pooled| - sum (n_l - 1) ln|S_l|
#   u = [sum 1/(n_l - 1) - 1/(n - g)] (2p^2 + 3p - 1) / [6 (p + 1)(g - 1)]
#   C = (1 - u) M  ~  chi^2 com  p(p+1)(g-1)/2  g.l.  (aprox.)
# ATENÇÃO: é muito sensível à não-normalidade (curtose) — ver discussão.
# -----------------------------------------------------------------------------
box_m <- function(X, grupo) {
  X <- as.matrix(X)
  grupo <- droplevels(as.factor(grupo))
  p <- ncol(X); g <- nlevels(grupo); n <- nrow(X)
  S_l <- lapply(split.data.frame(X, grupo), cov)
  n_l <- as.vector(table(grupo))
  S_p <- Reduce(`+`, Map(function(S, m) (m - 1) * S, S_l, n_l)) / (n - g)
  logdet <- function(A) as.numeric(determinant(A, logarithm = TRUE)$modulus)
  M <- (n - g) * logdet(S_p) - sum((n_l - 1) * vapply(S_l, logdet, 0))
  u <- (sum(1 / (n_l - 1)) - 1 / (n - g)) *
       (2 * p^2 + 3 * p - 1) / (6 * (p + 1) * (g - 1))
  C <- (1 - u) * M
  gl <- p * (p + 1) * (g - 1) / 2
  list(M = M, u = u, C = C, gl = gl,
       p_valor = pchisq(C, gl, lower.tail = FALSE),
       S_l = S_l, S_pooled = S_p,
       log_det_S_l = vapply(S_l, logdet, 0), log_det_S_pooled = logdet(S_p))
}

# -----------------------------------------------------------------------------
# M de Box com valor crítico por bootstrap dos resíduos agrupados
# (Zhang & Boos, 1992, JASA 87, 425-429). Os resíduos de cada grupo
# (x_lj - xbar_l) são reunidos num só conjunto e reamostrados com reposição
# para todos os grupos. Isso impõe H0 (mesma distribuição de resíduos, logo
# mesma covariância) mas preserva a curtose observada, que é justamente o que
# distorce a referência qui-quadrado do M de Box.
# p-valor = (1 + #{C* >= C_obs}) / (R + 1)
# -----------------------------------------------------------------------------
box_m_bootstrap <- function(X, grupo, R = 9999, semente = 731) {
  X <- as.matrix(X)
  grupo <- droplevels(as.factor(grupo))
  n_l <- table(grupo)
  C_obs <- box_m(X, grupo)$C
  medias <- do.call(rbind, lapply(split.data.frame(X, grupo), colMeans))
  E <- X - medias[as.character(grupo), ]
  gr <- rep(factor(levels(grupo), levels = levels(grupo)), times = n_l)
  set.seed(semente)
  C_boot <- vapply(seq_len(R), function(r)
    box_m(E[sample.int(nrow(E), nrow(E), replace = TRUE), , drop = FALSE], gr)$C, 0)
  list(C_obs = C_obs, C_boot = C_boot, R = R, semente = semente,
       quantil_95 = unname(quantile(C_boot, 0.95)),
       p_valor = (1 + sum(C_boot >= C_obs)) / (R + 1))
}

# -----------------------------------------------------------------------------
# Assimetria e curtose multivariadas de Mardia (1970)
# [MKB] Seção 1.8 (medidas b_{1,p} e b_{2,p}) e Cap. 5 (testes de
# multinormalidade). Usa S_n (divisor n), como na definição original.
#   g_ij = (x_i - xbar)' S_n^{-1} (x_j - xbar)
#   b1 = n^{-2} sum_ij g_ij^3 ;  n b1 / 6 ~ chi^2_{p(p+1)(p+2)/6}
#   b2 = n^{-1} sum_i g_ii^2 ;  (b2 - E b2) / sqrt(8p(p+2)/n) ~ N(0,1),
#        com E b2 = p(p+2)(n-1)/(n+1)
# -----------------------------------------------------------------------------
mardia <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  Xc <- sweep(X, 2, colMeans(X))
  Sn <- crossprod(Xc) / n
  G <- Xc %*% solve(Sn, t(Xc))
  b1 <- sum(G^3) / n^2
  b2 <- sum(diag(G)^2) / n
  est_ass <- n * b1 / 6
  gl_ass <- p * (p + 1) * (p + 2) / 6
  # média exata de b2 sob normalidade: p(p+2)(n-1)/(n+1) (Mardia, 1970);
  # a média assintótica p(p+2) vicia o z para n moderado (ex.: n = 68)
  e_b2 <- p * (p + 2) * (n - 1) / (n + 1)
  z_curt <- (b2 - e_b2) / sqrt(8 * p * (p + 2) / n)
  c(n = n, b1p = b1, qui2_assimetria = est_ass, gl = gl_ass,
    p_assimetria = pchisq(est_ass, gl_ass, lower.tail = FALSE),
    b2p = b2, esperado_b2p = e_b2, z_curtose = z_curt,
    p_curtose = 2 * pnorm(-abs(z_curt)))
}

# -----------------------------------------------------------------------------
# Distâncias de Mahalanobis quadradas e gráfico qui-quadrado
# [JW] Seção 4.6, eq. (4-32) e Exemplo 4.13; Result 4.7:
#   (X - mu)' Sigma^{-1} (X - mu) ~ chi^2_p  sob normalidade
# Quantis usados: chi^2_p((j - 1/2)/n), como em [JW] p. 197.
# -----------------------------------------------------------------------------
dist_mahalanobis2 <- function(X, centro = colMeans(X), S = cov(X)) {
  stats::mahalanobis(as.matrix(X), centro, S)
}

grafico_quiquadrado <- function(d2, p, main = "", col = "black", pch = 19,
                                rotular_n = 3) {
  n <- length(d2)
  ord <- order(d2)
  # cores/símbolos por ponto precisam seguir a MESMA ordenação de d2
  if (length(col) == n) col <- col[ord]
  if (length(pch) == n) pch <- pch[ord]
  q <- qchisq((seq_len(n) - 0.5) / n, df = p)
  plot(q, d2[ord], xlab = bquote(chi[.(p)]^2 * " quantil  ((j - 1/2)/n)"),
       ylab = expression(d[j]^2 ~ "ordenado"), main = main,
       pch = pch, col = col, cex = 0.8, las = 1)
  abline(0, 1, lty = 2, col = "grey40")
  if (rotular_n > 0) {
    top <- tail(seq_len(n), rotular_n)
    text(q[top], d2[ord][top], labels = names(d2)[ord][top] %||% ord[top],
         pos = 2, cex = 0.7, col = "grey20")
  }
  # correlação dos pontos do gráfico (análoga ao r_Q de [JW] (4-31))
  invisible(cor(q, d2[ord]))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# -----------------------------------------------------------------------------
# Intervalos simultâneos de Bonferroni para tau_ki - tau_li
# [JW] Result 6.5 (p. 329) e eq. (6-42): m = p g (g - 1)/2 afirmações
#   (xbar_ki - xbar_li) +- t_{n-g}(alpha / (p g (g-1))) *
#                          sqrt( w_ii/(n - g) * (1/n_k + 1/n_l) )
# -----------------------------------------------------------------------------
ic_bonferroni <- function(s, alpha = 0.05) {
  niveis <- rownames(s$medias)
  m <- s$p * s$g * (s$g - 1) / 2
  tcrit <- qt(1 - alpha / (2 * m), df = s$n - s$g)
  pares <- utils::combn(niveis, 2)
  res <- list()
  for (j in seq_len(ncol(pares))) {
    k <- pares[1, j]; l <- pares[2, j]
    for (i in seq_len(s$p)) {
      dif <- s$medias[k, i] - s$medias[l, i]
      ep <- sqrt(s$W[i, i] / (s$n - s$g) * (1 / s$n_l[k] + 1 / s$n_l[l]))
      res[[length(res) + 1]] <- data.frame(
        variavel = colnames(s$medias)[i], contraste = paste(k, "-", l),
        diferenca = dif, ep = ep, li = dif - tcrit * ep, ls = dif + tcrit * ep,
        exclui_zero = (dif - tcrit * ep > 0) | (dif + tcrit * ep < 0))
    }
  }
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  attr(out, "m") <- m
  attr(out, "t_critico") <- tcrit
  out
}

# -----------------------------------------------------------------------------
# Teste de permutação para Lambda de Wilks (não-paramétrico; não exige
# normalidade, mas exige permutabilidade sob H0 -> mesma distribuição nos
# grupos, o que inclui covariâncias iguais). T = B + W é invariante às
# permutações, logo basta recalcular |W|.
# p-valor = (1 + #{Lambda_perm <= Lambda_obs}) / (R + 1)
# -----------------------------------------------------------------------------
permutacao_wilks <- function(X, grupo, R = 9999, semente = 731) {
  X <- as.matrix(X)
  grupo <- droplevels(as.factor(grupo))
  logdet <- function(A) as.numeric(determinant(A, logarithm = TRUE)$modulus)
  Xc <- sweep(X, 2, colMeans(X))
  ldT <- logdet(crossprod(Xc))
  ldW <- function(gr) {
    W <- matrix(0, ncol(X), ncol(X))
    for (lv in levels(gr)) {
      Xl <- X[gr == lv, , drop = FALSE]
      W <- W + crossprod(sweep(Xl, 2, colMeans(Xl)))
    }
    logdet(W)
  }
  L_obs <- exp(ldW(grupo) - ldT)
  set.seed(semente)
  L_perm <- vapply(seq_len(R), function(r) exp(ldW(sample(grupo)) - ldT), 0)
  list(L_obs = L_obs, L_perm = L_perm, R = R, semente = semente,
       p_valor = (1 + sum(L_perm <= L_obs)) / (R + 1))
}

# -----------------------------------------------------------------------------
# Variáveis canônicas discriminantes (direções de W^{-1}B)
# [JW] Seção 11.7 (Fisher para g populações) — usadas aqui apenas como
# ferramenta DESCRITIVA para visualizar a separação detectada pela MANOVA.
# Normalização: a' (W/(n-g)) a = 1  ([JW] Seção 11.7).
# -----------------------------------------------------------------------------
variaveis_canonicas <- function(X, s) {
  Sp <- s$W / (s$n - s$g)
  eS <- eigen(Sp, symmetric = TRUE)
  S_mhalf <- eS$vectors %*% diag(1 / sqrt(eS$values)) %*% t(eS$vectors)
  Bs <- s$B / (s$n - s$g)
  M <- S_mhalf %*% Bs %*% S_mhalf
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  r <- min(s$p, s$g - 1)
  A <- S_mhalf %*% e$vectors[, seq_len(r), drop = FALSE]
  dimnames(A) <- list(colnames(X), paste0("CAN", seq_len(r)))
  escores <- sweep(as.matrix(X), 2, s$xbar) %*% A
  list(coef = A, escores = escores,
       prop = e$values[seq_len(r)] / sum(e$values[seq_len(r)]))
}

# -----------------------------------------------------------------------------
# Utilitário: salvar tabela CSV com arredondamento controlado
# -----------------------------------------------------------------------------
salvar_tabela <- function(df, nome, dir = "results/tables", digitos = 4) {
  df <- as.data.frame(df)
  num <- vapply(df, is.numeric, TRUE)
  # formatação explícita evita artefatos de ponto flutuante no CSV
  # (ex.: 8.61970999999999e-244). Obs.: p-valores < 1e-16 não têm
  # significado numérico além de "praticamente zero" (precisão da máquina).
  txt <- vapply(df, is.character, TRUE) | vapply(df, is.factor, TRUE)
  df[num] <- lapply(df[num], function(v)
    ifelse(is.na(v), NA, trimws(formatC(v, digits = digitos + 2, format = "g"))))
  utils::write.csv(df, file.path(dir, nome), row.names = TRUE,
                   quote = which(txt))  # aspas só em colunas de texto (e nomes)
  invisible(df)
}

# #############################################################################
# PARTE B — ANÁLISE
# #############################################################################

R_BOOT_BOX <- 9999   # réplicas do bootstrap do M de Box (Seção 3e)
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
URL_DADOS <- paste0("https://raw.githubusercontent.com/allisonhorst/palmerpenguins/",
                    "main/inst/extdata/penguins.csv")
MD5_ESPERADO <- "a06a0210251465a86fb970018292304d"
candidatos <- c("penguins.csv", "data/raw/penguins.csv", "../data/raw/penguins.csv")
arq <- candidatos[file.exists(candidatos)][1]
if (is.na(arq)) {
  arq <- "penguins.csv"
  cat("penguins.csv não encontrado; baixando de", URL_DADOS, "\n")
  utils::download.file(URL_DADOS, arq, mode = "wb", quiet = TRUE)
}
cat("Arquivo de dados:", arq, "\n")
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

# 3e. M de Box com referência por bootstrap (Zhang & Boos, 1992): separa
#     heterogeneidade de curtose, porque a distribuição de referência é gerada
#     com a curtose real dos resíduos e sob covariância comum.
bmb <- box_m_bootstrap(X, grupo, R = R_BOOT_BOX, semente = SEMENTE)
cat(sprintf("Box M bootstrap (R = %d): C_obs = %.2f; quantil 95%% bootstrap = %.2f (qui2: %.2f); p = %.5f\n",
            bmb$R, bmb$C_obs, bmb$quantil_95, qchisq(0.95, bm$gl), bmb$p_valor))
salvar_tabela(data.frame(C_obs = bmb$C_obs, quantil_95_bootstrap = bmb$quantil_95,
                         quantil_95_qui2 = qchisq(0.95, bm$gl),
                         max_C_bootstrap = max(bmb$C_boot), p_valor = bmb$p_valor,
                         R = bmb$R), "t21_box_m_bootstrap.csv")

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

# eta^2 univariado de cada resposta, b_ii / t_ii (ANOVA de cada variável),
# como termo de comparação concreto para o eta^2 multivariado.
eta2_uni <- diag(s$B) / diag(s$T)
cat("eta^2 univariado (b_ii / t_ii):\n"); print(round(eta2_uni, 4))
salvar_tabela(data.frame(eta2 = c(eta2_uni, multivariado = eta2_mult)),
              "t20_eta2.csv")

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

# Cenários que separam VOLUME e ORIENTAÇÃO da heterogeneidade de (c).
# Rodam num laço próprio, com semente própria, para não alterar (a)-(d).
#   (e) Sigma_l = c_l S_pooled, com |Sigma_l| = |S_l|: volumes observados,
#       mesma forma e orientação (o caso coberto pela heurística univariada)
#   (f) Sigma_l = S_l / k_l,   com |Sigma_l| = |S_pooled|: formas e
#       orientações observadas, volumes iguais
c_l <- exp((bm$log_det_S_l - bm$log_det_S_pooled) / p)
cat("Fatores de volume c_l = (|S_l|/|S_pooled|)^(1/p):\n"); print(round(c_l, 4))
Sl_vol_chol <- lapply(levels(grupo), function(l) chol(c_l[[l]] * Sp))
Sl_ori_chol <- lapply(levels(grupo), function(l) chol(bm$S_l[[l]] / c_l[[l]]))
names(Sl_vol_chol) <- names(Sl_ori_chol) <- levels(grupo)
gera <- function(chols) do.call(rbind, lapply(levels(grupo), function(l)
  matrix(rnorm(n_l[l] * p), n_l[l], p) %*% chols[[l]]))
set.seed(SEMENTE + 1)
sim_rej2 <- t(replicate(R_SIM, c(
  e_exato_volume_sem_orientacao = rejeita(gera(Sl_vol_chol))$exato[["p_valor"]] < ALFA,
  f_exato_orientacao_sem_volume = rejeita(gera(Sl_ori_chol))$exato[["p_valor"]] < ALFA)))
tam2 <- colMeans(sim_rej2)
ep2 <- sqrt(tam2 * (1 - tam2) / R_SIM)
tab_mc2 <- data.frame(tamanho_empirico = tam2, ic95_li = tam2 - 1.96 * ep2,
                      ic95_ls = tam2 + 1.96 * ep2, R = R_SIM)
print(round(tab_mc2, 4))
salvar_tabela(tab_mc2, "t22_monte_carlo_volume_orientacao.csv")

# ---- 9. Registro do ambiente ------------------------------------------------
secao("9. AMBIENTE")
print(sessionInfo())
cat("\nExecução concluída sem erros.\n")
sink()
close(log_con)
setwd(.wd_original)
