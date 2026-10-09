# =============================================================================
# funcoes.R — Funções auxiliares do P1 (ME731): MANOVA a um fator
#
# Tudo implementado em R base (pacotes 'stats', 'graphics', 'grDevices',
# 'utils', 'tools'), sem dependências externas, para maximizar a
# reprodutibilidade. Cada função cita a fonte teórica correspondente:
#   [JW]  Johnson & Wichern, Applied Multivariate Statistical Analysis, 4a ed.
#   [JW6] Johnson & Wichern, 6a ed. (bibliografia [1.A] do PDD)
#   [MKB] Mardia, Kent & Bibby, Multivariate Analysis (1979)
# =============================================================================

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
