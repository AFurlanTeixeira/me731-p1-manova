# Derivações e resultados teóricos a apresentar

Legenda de confiabilidade das referências:

- ✅ **J&W 4ª ed.**: seção, página e equação conferidas no PDF do livro-texto (`Livro_ME731.pdf`, 4ª ed., 1998).
- 🔶 **J&W 6ª ed.** e **MKB** (Mardia, Kent & Bibby, 1979)

A coluna "Nível" sugere o que fazer com cada item: **D** = demonstrar no texto, **E** = enunciar e citar, **V** = verificar numericamente no script.

---

## A. Fundamentos (pré-requisitos)

| # | Resultado | Onde | Nível |
|---|---|---|---|
| A1 | Se X ~ Nₚ(μ, Σ), então (X−μ)'Σ⁻¹(X−μ) ~ χ²ₚ | ✅ J&W Result 4.7 (Seç. 4.2) | E |
| A2 | Combinações lineares de normais multivariadas são normais | ✅ J&W Result 4.3 (Seç. 4.2) | E |
| A3 | EMV de μ e Σ; lema de maximização (1/\|Σ\|ᵇ) exp(−tr(Σ⁻¹B)/2) | ✅ J&W Results 4.10 e 4.11 (Seç. 4.3) | E (usado em C1) |
| A4 | X̄ ~ Nₚ(μ, Σ/n); (n−1)S ~ Wₙ₋₁(Σ); X̄ e S independentes | ✅ J&W Seç. 4.4, (4-23) · ✅ MKB Seç. 3.4 | E |
| A5 | Wishart: definição W_m(Σ) = dist. de Σ ZⱼZⱼ'; graus de liberdade se somam; CAC' ~ W_m(CΣC') | ✅ J&W (4-22), (4-24) · ✅ MKB Seç. 3.4 | E |
| A6 | Densidade de Wishart (só existe se n > p) | ✅ J&W (4-25) | E (justifica n − g ≥ p) |
| A7 | LGN e TCL multivariados (base da robustez assintótica) | ✅ J&W Seç. 4.5, (4-26), (4-27); Result 4.13 (TCL) | E |
| A8 | TRV geral: −2 ln Λ ≈ χ²_{ν−ν₀} | ✅ J&W (5-16) e Result 5.2 (Seç. 5.3) | E |

## B. O modelo MANOVA a um fator

| # | Resultado | Onde | Nível |
|---|---|---|---|
| B1 | Modelo X_ℓj = μ + τ_ℓ + e_ℓj, com Σ n_ℓ τ_ℓ = 0; e_ℓj iid Nₚ(0, Σ) | ✅ J&W (6-34), Seç. 6.4, p. 314–320 · ✅ MKB Cap. 12 (Seç. 12.3, one-way) | E |
| B2 | Suposições 1–3 (independência, Σ comum, normalidade) e relaxamento da 3 via TCL | ✅ J&W p. 314 | E |
| B3 | Decomposição x_ℓj = x̄ + (x̄_ℓ − x̄) + (x_ℓj − x̄_ℓ) | ✅ J&W (6-35) | E |
| B4 | **T = B + W** | ✅ J&W (6-36), (6-37) | **D** + V |
| B5 | W = Σ (n_ℓ − 1) S_ℓ, generalização do S_pooled de duas amostras | ✅ J&W (6-37) | D (uma linha) |
| B6 | Tabela MANOVA com g.l. g−1, n−g, n−1 | ✅ J&W p. 322 | E |
| B7 | MANOVA como modelo linear multivariado X = Zβ + E; TRV e sua distribuição | ✅ J&W Seç. 7.7 (Result 7.11) e Suplemento 7A | E |

**Esboço de B4.** Escreva x_ℓj − x̄ = (x_ℓj − x̄_ℓ) + (x̄_ℓ − x̄) e expanda o produto externo. A soma em j dos termos cruzados é zero, porque Σ_j (x_ℓj − x̄_ℓ) = 0 (J&W p. 321).

## C. Distribuição sob H₀ e o teste

| # | Resultado | Onde | Nível |
|---|---|---|---|
| C1 | **Λ* = \|W\|/\|B+W\| é função do TRV**: Λ_TRV = (Λ*)^{n/2} | ✅ J&W (6-38) e p. 322 ("related to the likelihood ratio criterion") · ✅ MKB Seç. 12.3 | **D** |
| C2 | W ~ Wₚ(n−g, Σ) sempre; sob H₀, B ~ Wₚ(g−1, Σ) e é independente de W | ✅ consequência de A4–A5 (W soma Wisharts independentes) · ✅ MKB Seç. 3.4 (teorema de Cochran) | D (W) / E (B) |
| C3 | Λ* ~ Λ(p, n−g, g−1) de Wilks; produto de betas independentes | ✅ MKB Seç. 3.7 | E |
| C4 | Λ* = ∏ 1/(1+λᵢ), λᵢ autovalores de W⁻¹B, s = min(p, g−1) | ✅ J&W nota 2, p. 322 | **D** + V |
| C5 | Pillai tr[B(B+W)⁻¹], Lawley–Hotelling tr[BW⁻¹], Roy (maior raiz) | ✅ J&W Seç. 6.9, p. 357 · ✅ MKB Seç. 12.3 (união-interseção → Roy) | E + V |
| C6 | **Exato para g = 3**: ((n−p−2)/p)·(1−√Λ*)/√Λ* ~ F_{2p, 2(n−p−2)} | ✅ J&W Tabela 6.3, p. 323 | E + V |
| C7 | **Bartlett**: −(n−1−(p+g)/2) ln Λ* ≈ χ²_{p(g−1)} | ✅ J&W (6-39), (6-40) | E + V (Monte Carlo) |
| C8 | **Invariância afim**: Λ* invariante sob X → AX + c, A não singular | consequência de (6-38) e de \|A'MA\| = \|A\|²\|M\| | **D** + V |

**Esboço de C1.** Sob H₁ (médias livres e Σ comum), o EMV é Σ̂ = W/n. Sob H₀ (média comum), Σ̂₀ = (B+W)/n = T/n. Pelo Result 4.10, o máximo da verossimilhança é proporcional a |Σ̂|^{−n/2}. Logo Λ_TRV = (|W|/|B+W|)^{n/2}. É o mesmo raciocínio de (5-13)–(5-14), que dá a relação T² ↔ Λ para uma amostra.

**Esboço de C4.** |W|/|B+W| = 1/|I + W⁻¹B| = ∏ 1/(1+λᵢ). Os λᵢ não nulos são no máximo s = posto(B) = min(p, g−1).

**Esboço de C8.** Com Y = XA', B_Y = A B_X A' e W_Y = A W_X A'. Então |W_Y|/|B_Y+W_Y| = |A|²|W_X| / (|A|²|B_X+W_X|) = Λ*_X. Consequência prática: trocar g por kg ou mm por cm não altera nada (o script verifica). Já variáveis canônicas e ICs dependem da escala.

> ⚠️ **Erro de impressão na 4ª ed. (p. 357):** Roy aparece como "maximum eigenvalue of W(B+W)⁻¹". O correto é o maior autovalor de **B**(B+W)⁻¹, que vale θ₁ = λ₁/(1+λ₁). Com p > s, como aqui, o maior deles vale 1 em qualquer amostra.  Atenção também: o R (`summary.manova`) reporta para Roy o próprio λ₁ = 15,02, e não θ₁ = 0,938.

## D. Diagnóstico das suposições

| # | Resultado | Onde | Nível |
|---|---|---|---|
| D1 | Q-Q plot e r_Q (correlação no Q-Q) | ✅ J&W Seç. 4.6, (4-29)–(4-31), Tabela 4.2 | E |
| D2 | **Gráfico qui-quadrado** com d²ⱼ = (xⱼ−x̄)'S⁻¹(xⱼ−x̄) contra χ²ₚ((j−½)/n) | ✅ J&W (4-32), Exemplo 4.13 | E (justificar com A1) |
| D3 | Detecção de atípicos multivariados | ✅ J&W Seç. 4.7 | E |
| D4 | Box–Cox (log = λ 0) | ✅ J&W Seç. 4.8, (4-34)–(4-35) | E (sensibilidade) |
| D5 | **M de Box**: C = (1−u)M ≈ χ²_{p(p+1)(g−1)/2} | ✅ J&W **6ª ed.**, Seç. 6.6 ("Testing for Equality of Covariance Matrices"). **Não está na 4ª ed.** · Box (1949) · ✅ MKB Cap. 5 (TRV para Σ₁ = … = Σ_g) | E |
| D6 | Assimetria b₁,ₚ e curtose b₂,ₚ de Mardia e seus testes (o script usa E b₂,ₚ = p(p+2)(n−1)/(n+1), a média exata) | ✅ MKB Seç. 1.8 (definições) e Cap. 5 (testes de multinormalidade) · Mardia (1970) | E |
| D7 | Robustez: T² e MANOVA pouco afetados por não-normalidade leve com n grande | ✅ J&W p. 313 (fim da Seç. 6.3) e p. 314 (TCL) | E |
| D8 | Pillai é um pouco mais robusto à **não-normalidade** (J&W) e à heterogeneidade de Σ (Olson) | ✅ J&W p. 357 (não-normalidade) · Olson (1974) (Σ heterogêneas) | E |

## E. Análise de acompanhamento

| # | Resultado | Onde | Nível |
|---|---|---|---|
| E1 | Var(τ̂_kᵢ − τ̂_ℓᵢ) = (1/n_k + 1/n_ℓ)σᵢᵢ, estimada por wᵢᵢ/(n−g) | ✅ J&W Seç. 6.5, p. 329 | D (curto) |
| E2 | **ICs de Bonferroni simultâneos**, m = pg(g−1)/2 | ✅ J&W Result 6.5, (6-42); desigualdade (5-28) | E + V |
| E3 | Variáveis canônicas de Fisher: autovetores de W⁻¹B (descritivo) | ✅ J&W Seç. 11.7 · ✅ MKB Seç. 12.5 (canonical variables, test of dimensionality) | E |
| E4 | MANOVA a dois fatores (extensão: espécie × sexo) | ✅ J&W Seç. 6.6, (6-43)–(6-46) | E |
| E5 | Estratégia: teste global → Bonferroni; cuidado com efeitos diluídos | ✅ J&W Seç. 6.9, p. 357–358 | E |

---

## Prioridade (se o espaço apertar)

1. **Indispensáveis:** B1, B2, B4, C1, C4, C6/C7, D2, D5, E2.
2. **Diferenciais de rigor:** C2, C8, a correção de Roy, a explicação da Monte Carlo para C7.
3. **Opcionais:** C3, A6, E3, E4.

## Cuidados conceituais (erros comuns nessa técnica)

1. **Normalidade é por grupo** (ou dos resíduos), não da amostra agregada. A amostra agregada de três espécies é uma mistura, e testá-la como se fosse uma única normal é um erro.
2. **Rejeitar o M de Box não proíbe a MANOVA.** O M de Box é muito sensível à curtose e, com n grande, detecta diferenças pequenas. O que importa é o efeito da violação sobre o tamanho do teste. A regra univariada (o maior n_ℓ com a maior variância torna o teste conservador) é só uma tendência no caso multivariado. Aqui, apesar de Adelie ter n e |S| maiores, a Monte Carlo deu tamanho 0,056, levemente liberal. Os cenários (e) e (f) da Monte Carlo (t22) mostram por quê: só volumes diferentes dão 0,048, e só formas diferentes dão 0,059.
3. **"p-valor minúsculo" não é tamanho de efeito.** Reporte η² multivariado e intervalos.
4. **ANOVAs univariadas depois da MANOVA não "explicam" o efeito multivariado.** Elas ignoram as correlações. Use ICs simultâneos e variáveis canônicas.
5. **A MANOVA não é "outra técnica" em relação à regressão multivariada.** É o mesmo modelo linear com regressoras indicadoras (J&W Seç. 7.7).
6. **Independência não se testa com esses dados.** Argumente pelo delineamento.
