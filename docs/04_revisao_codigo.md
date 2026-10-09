# Revisão do código: problemas encontrados e correções

Registro da revisão em duas rodadas: autorrevisão e uma auditoria independente que reimplementou as fórmulas. Serve de evidência de rigor e de guia se você alterar o script.

## Fórmulas conferidas por implementação independente

B e W (B + W = T); Wilks, Pillai, Lawley–Hotelling e Roy; F exato da Tabela 6.3; Bartlett (6-39); M de Box (log-determinante via autovalores, concordância de 1e-9); Mardia (laço duplo explícito); Bonferroni (t₃₃₉(0,05/24) = 2,8848); p-valor de permutação; variáveis canônicas (A'S_pA = I); quantis do gráfico qui-quadrado. Além disso, o script confere com `stopifnot()` os valores feitos à mão contra `stats::manova()`.

## Erros corrigidos

| # | Problema | Tipo | Correção |
|---|---|---|---|
| 1 | `source()` falhava em silêncio em locale não-UTF-8 (`LANG=C`): nenhuma função carregava e o erro só aparecia adiante | reprodutibilidade | o script ativa um locale UTF-8 e confere com `stopifnot(exists(...))` |
| 2 | No gráfico qui-quadrado, cores e símbolos não seguiam a ordenação de d², o que atribuía espécies erradas aos pontos (obs15, que é Adelie, aparecia como Chinstrap) | **erro de figura** | `col` e `pch` reordenados com `ord` em `grafico_quiquadrado()` |
| 3 | Comparação com `manova()` usava tolerância absoluta; falhava para Lawley–Hotelling ≈ 17 | falso alarme | tolerância relativa (`all.equal`) |
| 4 | O R reporta para Roy λ₁ = 15,02, e não θ₁ = λ₁/(1+λ₁) | interpretação | conferência ajustada e nota no código e nos docs |
| 5 | `line = 1` passado a `pairs()` vazava para os painéis (44 warnings) | warnings | removido; o script roda com `options(warn = 2)` |
| 6 | Legenda da matriz de dispersão sobreposta ao título | figura | legenda numa camada própria |
| 7 | CSVs com artefatos de ponto flutuante (8.61970999999999e-244) e números entre aspas | apresentação | `formatC` e aspas só em colunas de texto |
| 8 | **Monte Carlo com R = 2.000 sustentava a conclusão "teste conservador" (0,042)**. Com R = 20.000 o tamanho é 0,056 [0,053; 0,059] | **conclusão errada** | R_SIM = 20.000, IC de Monte Carlo, texto reescrito |
| 9 | Monte Carlo citada como evidência de robustez à não-normalidade, mas gerava dados normais | argumento inválido | novo cenário: bootstrap dos resíduos reais (0,052) |
| 10 | Curtose de Mardia com média assintótica p(p+2), que vicia o z para n = 68 | precisão | média exata p(p+2)(n−1)/(n+1) |
| 11 | Pillai atribuído a J&W como "robusto à heterogeneidade de Σ"; o livro diz "não-normalidade" | citação | corrigido (J&W para não-normalidade, Olson para Σ) |
| 12 | Sensibilidade (log, sem atípico) sem p-valor | completude | p-valores exatos adicionados (t18) |
| 13 | `source()` mudava o diretório de trabalho do usuário | efeito colateral | o diretório original é restaurado ao final |
| 14 | Resíduos padronizados por S_pooled apesar de as Σ_ℓ diferirem | interpretação | comentário no código e no roteiro |

## Testes de reprodutibilidade feitos

- Execução limpa com `Rscript` a partir da raiz, de outro diretório e via `source()`.
- `LANG=C LC_ALL=C` e `options(warn = 2)`: zero warnings.
- Duas execuções produzem tabelas, figuras e log byte a byte idênticos; um clone novo regenera exatamente os mesmos arquivos.
- Cópia independente em outro diretório produz tabelas idênticas.

## Limitações conhecidas (não corrigidas)

- Se o script der erro no meio, o `sink()` fica aberto na sessão interativa. A próxima execução o fecha automaticamente; ou rode `sink()` à mão.
- Os PNGs podem diferir em bytes entre sistemas operacionais (fontes, antialiasing). As tabelas não diferem.
- Os resultados numéricos podem variar na última casa com outra BLAS/LAPACK. A Monte Carlo usa `rnorm` e `sample`, estáveis entre versões do R ≥ 3.6 (`sample.kind = "Rejection"`).
