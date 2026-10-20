# cpt
Implementação em C++ da Teoria de Potência CPT

## C++

A biblioteca é o header `cpt.h` (C++20), sem dependências e sem alocação dinâmica.

É necessário copiar o mingw na pasta do projeto para a compilação do arquivo de teste final.psimsch.

## Python

O módulo `cpt` é um wrapper para o `cpt.h`: os cálculos são feitos pelo mesmo código C++ usado no firmware.

Instalação (requer um compilador C++20 e CMake ≥ 3.18):

```sh
pip install .            # ou: pip install .[test] para incluir o pytest
```

Uso:

```python
import numpy as np
import cpt

c = cpt.CPT(channels=3, sampling_rate=12000, freq=60)

# u e i com shape (N, 3); para monofásico, (N,) ou (N, 1)
r = c.feed_block(u, i)
r['ia']                  # corrente ativa por fase, shape (N, 3)
r['totals']['P']         # potência ativa total, shape (N,)

# Só algumas grandezas, uma vez por ciclo da rede
r = c.feed_block(u, i, fields=['t.P', 't.Q', 't.V'], step=12000 // 60)
r['n']                   # número da amostra de cada linha

# Uma amostra por vez (API anterior, bem mais lenta)
r = c.feed([va, vb, vc], [ia, ib, ic])
```

O estado é mantido entre chamadas, então um sinal longo pode ser processado em blocos; `reset()` volta ao estado inicial.

Grandezas por fase: `u i p P w W U ia ir iv iba ibr iua iur`. Totais: `P Q N V A rf uf df pf` (em `fields`, com prefixo `t.`). No monofásico, `iba ibr iua iur N uf` não são definidas.

O template `CPT` é resolvido em tempo de compilação, então o módulo traz um conjunto fixo de configurações: 1 e 3 fases, `double` e `float`, nas taxas abaixo. `cpt.configurations()` lista todas. Para incluir outra, acrescente uma linha `bind_rate<taxa, freq>` no fim de `python/_cpt.cpp`.

| Rede | Taxas de amostragem (Hz) |
| --- | --- |
| 60 Hz | 7680, 12000, 15360, 30720, 61440 |
| 50 Hz | 6400, 12000, 12800, 25600, 51200 |

Testes do wrapper:

```sh
pytest
```
