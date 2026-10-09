# Roteiro de seções do relatório P1

Cada seção abaixo diz **o que precisa demonstrar** e **qual evidência do repositório usar**. O texto é seu: este roteiro é um esqueleto de argumentação, não um rascunho para colar (ver item (n) do PDD).

Critérios do PDD que o roteiro cobre: apresentação formal do modelo, ilustração em dados (n ≥ 30, p ≥ 4, n > p), discussão crítica, rigor estatístico com referências, coerência textual ("o quê, por quê e como") e programa que reproduz tudo.

---

## 1. Introdução (≈ ½ página)

**Precisa demonstrar:** qual problema a MANOVA resolve e por que um teste multivariado é preferível a p ANOVAs.

- Pergunta: comparar g vetores de médias p-dimensionais.
- Dois argumentos: (i) controle do erro tipo I global; (ii) uso da correlação entre respostas. Cite J&W Seção 6.9 e o Exemplo 6.14, em que duas ANOVAs não rejeitam e o T² rejeita.
- Antecipe o dado (pinguins) e a estrutura do texto.

## 2. Modelo e formulação matricial (≈ 1½ página)

**Precisa demonstrar:** domínio da formulação, das hipóteses e da decomposição em somas de quadrados e produtos cruzados.

- Modelo X_ℓj = μ + τ_ℓ + e_ℓj, com a restrição Σ n_ℓ τ_ℓ = 0 (J&W (6-34)). Estimadores pela decomposição (6-35).
- As três suposições de J&W p. 314: independência, Σ comum, normalidade multivariada.
- Decomposição T = B + W (6-36), com graus de liberdade g−1, n−g e n−1. Tabela MANOVA (p. 322).
- Forma de modelo linear multivariado: **X** = **Z** **β** + **E**, em que **Z** é a matriz de delineamento com indicadoras. Mostra que a MANOVA é um caso particular da regressão multivariada (J&W Seção 7.7).
- **Derivação curta obrigatória:** B + W = T (expansão do produto cruzado; o termo cruzado se anula porque Σ_j (x_ℓj − x̄_ℓ) = 0).

## 3. Distribuições e teste (≈ 2 páginas). É o núcleo teórico.

**Precisa demonstrar:** de onde vem Λ, qual é a sua distribuição e por que os testes são equivalentes ou não.

- Wishart: definição (4-22), propriedades (4-24), (n−1)S ~ W (4-23).
- Sob normalidade: W ~ W_p(n−g, Σ). Sob H₀: B ~ W_p(g−1, Σ), independente de W.
- Λ* = |W| / |B + W| (6-38) como TRV: Λ_TRV = (Λ*)^(n/2). Derivar pela maximização da verossimilhança (Result 4.10 / 4.11), como em (5-13)–(5-14) para T².
- Λ* em termos dos autovalores de W⁻¹B (nota de rodapé 2, p. 322) e as outras três estatísticas (p. 357).
- Distribuição exata para g = 3 (Tabela 6.3) e aproximação de Bartlett (6-39)/(6-40). Ligação com o Result 5.2 (−2 ln Λ ~ χ² assintoticamente).
- **Propriedade de invariância:** Λ* não muda sob X → AX + c, com A não singular. Demonstração de duas linhas; o script a verifica numericamente.

## 4. Dados e análise exploratória (≈ 1 página)

**Precisa demonstrar:** que os dados atendem aos requisitos do PDD e que você conhece a estrutura deles antes de modelar.

- Origem, licença, n_ℓ, faltantes (t01), por que island e year ficam de fora.
- Tabela t02 e figuras f01–f02.
- Correlações dentro dos grupos (t03): justificam a abordagem multivariada.

## 5. Verificação das suposições, antes do ajuste (≈ 1½ página)

**Precisa demonstrar:** que cada suposição foi checada com o diagnóstico adequado e que você sabe interpretar o resultado.

| Suposição | Diagnóstico | Evidência |
|---|---|---|
| Normalidade multivariada **por grupo** | Q-Q univariados, gráfico qui-quadrado (J&W 4.6, (4-32)), Mardia | f03, f04, t04, t05 |
| Atípicos | d² > χ²₄(0,001) dentro do grupo (J&W 4.7) | t06 (obs294, Chinstrap) |
| Σ comum | M de Box, \|S_ℓ\|, razão de desvios-padrão | t07, t08 |
| Independência | argumento de delineamento, não um teste | texto |

Ponto central: o M de Box rejeita. Discuta se a rejeição vem de heterogeneidade real ou da sensibilidade do teste à curtose. Os log|S_ℓ| diferem (Adelie 17,33 contra Gentoo 16,04), então a heterogeneidade existe. A curtose de Mardia em Chinstrap (p < 0,001) também contamina o M.

## 6. Resultados da MANOVA (≈ 1 página)

**Precisa demonstrar:** cálculo correto e leitura correta, não só "p < 0,05".

- B e W (t09, t10), autovalores 15,02 e 2,32, as quatro estatísticas (t11), F exato e Bartlett (t12).
- Reporte p < 0,001, não "4·10⁻²⁸⁴": p-valores abaixo de ~10⁻¹⁶ não têm significado numérico.
- **Tamanho de efeito** η² = 0,863.
- Conferência: valores feitos à mão iguais aos de `stats::manova()`.

## 7. Análise de acompanhamento (≈ 1 página)

**Precisa demonstrar:** como localizar as diferenças sem inflar o erro tipo I.

- ICs simultâneos de Bonferroni, Result 6.5, m = 12 (t13, f05). Adelie e Chinstrap não diferem em profundidade do bico nem em massa.
- Variáveis canônicas (f06, t14): o 1º eixo (86,6%) separa Gentoo, puxado pela profundidade do bico; o 2º separa Chinstrap pelo comprimento do bico. Use os coeficientes **padronizados**.

## 8. Diagnóstico pós-ajuste (≈ ½ página)

**Precisa demonstrar:** que os resíduos e_ℓj = x_ℓj − x̄_ℓ foram examinados (J&W p. 320: os resíduos verificam as suposições).

- f07 e t15: Mardia nos resíduos agrupados (curtose p = 0,016); obs294 com d² ≈ 29. Os resíduos são padronizados por S_pooled. Como as Σ_ℓ diferem, parte da curtose reflete a heterogeneidade e não só a não-normalidade.

## 9. Robustez e sensibilidade (≈ 1 página)

**Precisa demonstrar:** que a conclusão não depende de uma suposição violada.

- Permutação (f08), escala log, remoção do atípico (t18).
- Monte Carlo (t19, R = 20.000, IC 95% de Monte Carlo): o teste exato tem tamanho 0,048 [0,045; 0,051] sob as suposições, 0,056 [0,053; 0,059] com as Σ_ℓ reais (levemente **liberal**) e 0,051 [0,048; 0,055] no bootstrap dos resíduos (forma e heterogeneidade reais). Conclusão: tamanho próximo do nominal. A heurística univariada ("maior n com maior variância → conservador") **não** se confirmou aqui. Vale discutir isso.
- MANOVA espécie × sexo (t16): a interação é significativa. Mardia rejeita em 2 das 6 células (Chinstrap-fêmea e Adelie-macho, t17), e o M de Box ainda rejeita dentro de cada sexo. O sexo explica **parte** da não-normalidade e da heterogeneidade. Lembre também que as células, com n entre 34 e 73, têm menos poder.

## 10. Discussão crítica (≈ 1½ página). É exigida explicitamente pelo PDD.

Ver `03_discussao_critica.md`. Estrutura: vantagens, limitações gerais do método, limitações **deste** estudo e extensões ou alternativas.

## 11. Conclusão (≈ ⅓ página)

Responda o quê, por quê e como (item 5 do PDD).

## Referências e Apêndice

- Referências: as do README.
- Apêndice: como rodar o código (`Rscript R/p1_manova.R`), correspondência entre tabelas/figuras e seções, e `sessionInfo`.

---

### Checklist final antes de enviar (prazo: 09/10/2026, 15h59; o PDD recomenda enviar 2 h antes)

- [ ] Todo número no texto confere com `results/` (rode o script de novo antes de fechar o texto).
- [ ] Toda figura e tabela citada no texto existe e está numerada.
- [ ] O script enviado é o `.R`, não R Markdown: o item (n) do PDD pede para evitar R Markdown.
- [ ] Autoavaliação (A1) preenchida.
