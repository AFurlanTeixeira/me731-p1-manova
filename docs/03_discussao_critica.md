# Discussão crítica: insumos

Pontos para a seção de discussão, com os números que os sustentam. Organize e redija com suas palavras.

## 1. Vantagens do método

- **Controle do erro tipo I global:** um único teste no lugar de p = 4 ANOVAs (J&W Seç. 6.9).
- **Usa a estrutura de correlação:** detecta diferenças em direções que nenhuma variável isolada mostra (J&W Exemplo 6.14). Aqui as correlações dentro do grupo vão de 0,31 a 0,72 (t03).
- **Base teórica completa:** TRV, distribuição exata para g = 3, aproximação de Bartlett e quatro estatísticas que coincidem quando s = 1.
- **Invariância afim:** Λ* não depende das unidades (verificado: g→kg e mm→cm dão Λ* idêntico).
- **Robustez assintótica** à não-normalidade, via TCL (J&W p. 314). No bootstrap dos resíduos reais (não-normais e heterogêneos), o tamanho empírico foi 0,052 [0,048; 0,055].
- Gera subprodutos interpretáveis: variáveis canônicas e ICs simultâneos.

## 2. Limitações e desvantagens gerais

- **Σ comum é forte e o teste para ela é ruim.** O M de Box confunde heterogeneidade com curtose, e o efeito da violação depende de como n_ℓ e |S_ℓ| se combinam.
- **Responde só "existe diferença?"** O teste global não diz onde está a diferença, e o acompanhamento (Bonferroni) é conservador quando m cresce (aqui m = 12 e t crítico = 2,88).
- **Quatro estatísticas sem ordenação uniforme de poder** quando s > 1 (J&W p. 357): Roy domina com um único autovalor grande; Pillai é um pouco mais robusto (J&W p. 357; Olson, 1974). A escolha deve ser feita **antes** de ver os dados.
- **Exige n − g ≥ p** para W ser invertível. Não serve para p grande diante de n (dados de alta dimensão).
- **Sensível a atípicos:** médias e SQPC não são robustas.
- **Diluição:** muitas respostas inertes reduzem o poder para um efeito concentrado em uma variável (J&W p. 358).
- **Só lida com diferenças de locação:** não detecta grupos que diferem só em dispersão.

## 3. Limitações deste estudo, com evidências

| Problema | Evidência | Consequência |
|---|---|---|
| Σ_ℓ heterogêneas | Box C = 76,8, p ≈ 1,4·10⁻⁸; log\|S\| de 16,0 a 17,3 (t07, t08) | Monte Carlo (R = 20.000): tamanho 0,056 [0,053; 0,059] com as Σ_ℓ reais, levemente liberal, apesar de a heurística univariada prever um teste conservador. Com p-valor de ordem 10⁻²⁸⁴ a conclusão não muda, mas em um efeito marginal faria diferença. |
| Não-normalidade em Chinstrap e Gentoo | Mardia: assimetria p = 0,019 e 0,009; curtose Chinstrap p < 0,001 (t05) | Atenuada pelo TCL (n_ℓ ≥ 68). Permutação dá p = 0,0001. |
| **Confundimento com sexo** | Interação espécie×sexo: Λ = 0,889, p < 0,001 (t16). Por célula espécie×sexo, Mardia rejeita em 2 das 6 (Chinstrap-fêmea e Adelie-macho, t17), e o M de Box ainda rejeita dentro de cada sexo (p < 0,001). | Cada espécie é uma **mistura** de dois sexos, o que explica **parte** da não-normalidade e da heterogeneidade (as células têm menos poder, com n de 34 a 73). Modelo mais adequado: dois fatores. O efeito de espécie depende do sexo. |
| Atípico | obs294 (Chinstrap, bico de 58 mm), d² = 25,6 contra limite 18,5 | Removê-lo muda Λ de 0,01879 para 0,01867: irrelevante. |
| Testes de normalidade múltiplos | 12 Shapiro-Wilk sem correção (t04) | Uso **descritivo**. Não decida a partir de um único p < 0,05. |
| Medidas discretizadas | nadadeira em mm inteiros; massa em múltiplos de 25 g | Empates, que afetam Shapiro-Wilk e Q-Q. É um efeito pequeno. |
| Independência | 3 anos de coleta (2007–2009), ninhos amostrados por ilha | Possível dependência (mesmo indivíduo ou ninho em anos diferentes, efeito de ilha). Confira `Individual ID` em `palmerpenguins::penguins_raw` antes de afirmar qualquer coisa. Island é parcialmente confundida com espécie (Gentoo só em Biscoe, Chinstrap só em Dream). |
| H₀ obviamente falsa | Λ = 0,019, η² = 0,86 | A inferência é pouco informativa em si. Por isso a Monte Carlo sob H₀ é parte importante do projeto metodológico. |
| Unidades e escalas distintas | massa em g ≫ mm | Não afeta Λ (invariância). Afeta coeficientes canônicos brutos: use os padronizados. |

## 4. Extensões, alternativas e melhorias

1. **MANOVA a dois fatores (espécie × sexo)**, J&W Seç. 6.6. Já feita como sensibilidade. Atenção: com delineamento desbalanceado as SQPC sequenciais (tipo I) dependem da ordem dos fatores.
2. **MANCOVA** com massa corporal como covariável: compara a *forma* do bico e da nadadeira descontado o tamanho.
3. **Testes que não exigem Σ comum:** versões multivariadas de Welch/James (James, 1954; Johansen, 1980) ou o teste de Nel & Van der Merwe (1986) para g = 2.
4. **Pillai como estatística principal** quando Σ_ℓ diferem (Olson, 1974).
5. **Testes de permutação e bootstrap:** já incluídos (permutação). Atenção: a permutação também supõe permutabilidade, então não corrige a heterogeneidade de Σ.
6. **Estimação robusta** (MCD, S-estimadores) para B e W, reduzindo a influência de atípicos.
7. **Testes não-paramétricos por postos** (multivariados de Kruskal–Wallis; Puri & Sen, 1971). Estão na lista de técnicas do PDD.
8. **PERMANOVA** (Anderson, 2001): baseada em distâncias, para dados não-euclidianos.
9. **Alta dimensão** (p ≥ n): MANOVA regularizada ou testes do tipo Bai–Saranadasa.
10. **Análise discriminante** (J&W Cap. 11): passa de "diferem?" para "quão bem se classificam?". É o passo natural depois de f06.

## 5. Possíveis perguntas do professor (prepare as respostas)

- Por que Λ de Wilks e não Pillai, dado que Σ comum foi rejeitada?
- Por que o tamanho foi levemente liberal se o maior grupo tem a maior variância generalizada?
- O que muda se o maior grupo tiver a menor variância generalizada?
- Por que a distribuição é exata para g = 3? (Λ* é função de um único quadrado de beta; ver MKB Seç. 3.7.)
- Por que o teste de permutação não resolve a heterogeneidade de Σ?
- Com quatro estatísticas, como evitar escolher a que dá menor p-valor?
