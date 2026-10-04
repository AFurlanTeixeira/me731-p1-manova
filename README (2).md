# ME731 — Projeto 1: MANOVA a um fator

Aspectos metodológicos da **Análise de Variância Multivariada (MANOVA) a um fator**, ilustrados com medidas morfométricas de pinguins do Arquipélago Palmer (Antártida).

ME731 — Análise Multivariada, Unicamp, 2º semestre de 2026. Prof. Aluísio de Souza Pinheiro.

| | |
|---|---|
| **Técnica** | MANOVA a um fator, teste Λ de Wilks (J&W, Seção 6.4) |
| **Dados** | `palmerpenguins` (Horst, Hill & Gorman, 2020), licença CC0 |
| **Tamanho** | n = 342 casos completos, p = 4 respostas contínuas, g = 3 espécies (151 / 68 / 123) |
| **Requisitos do PDD** | n ≥ 30 ✔, p ≥ 4 ✔, n > p ✔ (e n − g ≥ p, para que **W** seja positiva definida) |
| **Linguagem** | R (≥ 4.0), **apenas R base**, sem nenhum pacote externo |

## Pergunta

Os vetores de médias de (comprimento do bico, profundidade do bico, comprimento da nadadeira, massa corporal) são iguais nas três espécies?

H₀: τ₁ = τ₂ = τ₃ = 0 no modelo **X**ₗⱼ = **μ** + **τ**ₗ + **e**ₗⱼ, com **e**ₗⱼ ~ Nₚ(**0**, **Σ**).

## Estrutura do repositório

```
.
├── R/
│   ├── p1_manova.R        # script principal: roda a análise de ponta a ponta
│   └── funcoes.R          # funções implementadas à mão (B, W, Wilks, Box M, Mardia...)
├── data/
│   ├── raw/penguins.csv   # dados originais, sem nenhuma alteração (MD5 conferido no script)
│   └── README.md          # origem, licença e dicionário de variáveis
├── results/
│   ├── tables/            # t01–t19 (.csv), geradas pelo script
│   ├── figures/           # f01–f08 (.png), geradas pelo script
│   └── log_execucao.txt   # toda a saída do console + sessionInfo()
├── docs/
│   ├── 01_roteiro_secoes.md          # o que cada seção do relatório precisa demonstrar
│   ├── 02_derivacoes_referencias.md  # resultados teóricos, com seção/equação do J&W e do MKB
│   ├── 03_discussao_critica.md       # vantagens, limitações, extensões (com os números obtidos)
│   └── 04_revisao_codigo.md          # problemas encontrados na revisão e como foram resolvidos
├── .gitattributes         # impede conversão CRLF dos dados (o MD5 continuaria válido)
├── .gitignore
└── LICENSE
```

## Como reproduzir

```bash
git clone <url-do-repositorio>
cd me731-p1-manova
Rscript R/p1_manova.R          # ~40 s; regenera results/ inteiro
```

No RStudio: abra a pasta do repositório e dê *Source* em `R/p1_manova.R` (o script localiza a raiz sozinho).

Garantias de reprodutibilidade:

- **Sem dependências**: só pacotes que já vêm com o R (`stats`, `graphics`, `grDevices`, `utils`, `tools`).
- **Semente única** (`SEMENTE <- 731`) para o teste de permutação e o estudo de Monte Carlo; duas execuções geram tabelas, figuras e log byte a byte idênticos (testado também num clone novo).
- **Integridade dos dados**: o MD5 de `penguins.csv` é verificado na leitura.
- **Autoverificação**: `stopifnot()` confere B + W = T, Λ via determinantes = Λ via autovalores, os valores implementados à mão contra `stats::manova()` e a invariância afim de Λ.
- **Locale**: se a sessão não for UTF-8 (ex.: `LANG=C`), o script tenta ativar um locale UTF-8 antes de ler os arquivos acentuados.
- Testado com `options(warn = 2)` (warnings viram erros): zero warnings.

## Resultados principais

| Quantidade | Valor |
|---|---|
| Λ de Wilks | 0,01879 |
| F exato (J&W Tabela 6.3, g = 3) | F(8, 672) = 528,9; p < 0,001 |
| Bartlett (6-39) | χ²₈ = 1341,5; p < 0,001 |
| Permutação (9.999 réplicas) | p = 0,0001 (mínimo possível com R = 9.999) |
| η² multivariado = 1 − Λ^(1/s) | 0,863 |
| M de Box | C = 76,8; gl = 20; p ≈ 1,4·10⁻⁸ → **Σ comum é rejeitada** |
| Tamanho empírico do teste exato sob H₀ (Monte Carlo, R = 20.000) | 0,048 (normal, Σ comum) · 0,056 (Σₗ reais) · 0,052 (bootstrap dos resíduos); nominal 0,05 |

A rejeição de H₀ é trivial: as espécies são visivelmente diferentes. O interesse do projeto está na metodologia, ou seja, nos diagnósticos, na violação de Σ comum, no confundimento com o sexo e no comportamento dos testes sob H₀. Veja `docs/03_discussao_critica.md`.

## Referências

- Johnson, R. A. & Wichern, D. W. (1998). *Applied Multivariate Statistical Analysis*, 4ª ed. Prentice-Hall. Roteiro do curso.
- Johnson, R. A. & Wichern, D. W. (2007). *Applied Multivariate Statistical Analysis*, 6ª ed. Pearson.
- Mardia, K. V., Kent, J. T. & Bibby, J. M. (1979). *Multivariate Analysis*. Academic Press.
- Box, G. E. P. (1949). A general distribution theory for a class of likelihood criteria. *Biometrika*, 36, 317–346.
- Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis with applications. *Biometrika*, 57, 519–530.
- Olson, C. L. (1974). Comparative robustness of six tests in multivariate analysis of variance. *JASA*, 69, 894–908.
- Gorman, K. B., Williams, T. D. & Fraser, W. R. (2014). Ecological sexual dimorphism and environmental variability within a community of Antarctic penguins (genus *Pygoscelis*). *PLoS ONE*, 9(3), e90081.
- Horst, A. M., Hill, A. P. & Gorman, K. B. (2020). *palmerpenguins: Palmer Archipelago (Antarctica) penguin data*. R package. doi:10.5281/zenodo.3960218.

## Licença

Código: MIT (ver `LICENSE`). Dados: CC0 (ver `data/README.md`).
