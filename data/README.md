# Dados

## Origem

`raw/penguins.csv` é cópia **inalterada** de
`https://raw.githubusercontent.com/allisonhorst/palmerpenguins/main/inst/extdata/penguins.csv`
(baixado em 2026-10-04).

- MD5: `a06a0210251465a86fb970018292304d`
- SHA-256: `f204db2c753b0937caac3cb35258562c14f073e4bbc76be24b4c51ce22767a93`
- Linhas: 344 observações + cabeçalho; 8 colunas.

Dados coletados por Dra. Kristen Gorman e pelo Palmer Station LTER (2007–2009).
Licença: **CC0 1.0** (domínio público). Citação: Horst, Hill & Gorman (2020); Gorman, Williams & Fraser (2014).

## Dicionário

| Coluna | Tipo | Uso na análise |
|---|---|---|
| `species` | categórica (Adelie, Chinstrap, Gentoo) | **fator** (g = 3) |
| `island` | categórica (Biscoe, Dream, Torgersen) | não usada (confundida com espécie: Gentoo só ocorre em Biscoe) |
| `bill_length_mm` | contínua, mm | resposta X₁ |
| `bill_depth_mm` | contínua, mm | resposta X₂ |
| `flipper_length_mm` | inteira, mm | resposta X₃ (medida arredondada ao mm) |
| `body_mass_g` | inteira, g | resposta X₄ (todas as medidas são múltiplas de 25 g) |
| `sex` | categórica (female, male), 11 NA | só na análise de sensibilidade (MANOVA espécie × sexo) |
| `year` | 2007, 2008, 2009 | não usada (ver nota sobre independência em `docs/03`) |

## Faltantes

2 linhas sem nenhuma das 4 medidas → excluídas (n = 342).
Mais 9 linhas sem `sex` são mantidas na análise principal e excluídas apenas na análise espécie × sexo (n = 333).
