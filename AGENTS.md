# Apunte EDO–EDP — notas para agentes

Apunte en español de métodos numéricos para EDO y EDP. Capítulos `NN-*.md`
numerados en orden pedagógico, figuras PNG en `figs/`, HTML en `html/`.

## Comandos

- **Regenerar figuras**: `uv run --with matplotlib --with numpy figs/generar_figuras.py`
  (cada función `fig_*` devuelve un check que se verifica con assert; agregar
  la función nueva a la lista `checks` del `__main__`).
- **Regenerar HTML**: `./build_html.sh`

## Gotchas

- `pandoc --mathjax` **sin URL** en Debian apunta a
  `/usr/share/javascript/mathjax` (MathJax 2 local, no instalado) y las
  fórmulas quedan como LaTeX crudo. `build_html.sh` ya usa la URL de CDN
  de MathJax 3 — no regenerar a mano sin esa URL.
- Los `src="figs/..."` del HTML se resuelven relativo a `html/`. Existe un
  symlink `html/figs -> ../figs` que hace que las imágenes carguen al abrir
  los HTML directo en el navegador — no borrarlo ni regenerar `html/` de
  forma que lo pise.
- Python del sistema no tiene matplotlib ni `python` (solo `python3`): usar
  siempre `uv run --with ...`.
- Pandoc no cierra math inline si el `$` de cierre va seguido de un dígito:
  `$\sim$50` renderiza como `$$50` literal y se come el `\sim`. Escribir
  `$\sim 50$` (el número dentro del math).
- El ancho de columna lo fija `style.css` (`--measure: 44rem` ≈ 790px) e
  `img{max-width:100%}`: una figura ancha y baja (p. ej. 2 paneles lado a
  lado, 1425x510) se reduce a ~55% — preferir apilar paneles en vertical
  (~960px de ancho) o subir fuentes.
- `build_html.sh` también genera `index.html` desde `README.md`
  (reescribe links `.md` → `.html` con sed) y la nav anterior/siguiente
  de cada capítulo. El CSS vive en `style.css` (raíz) y el script lo
  copia a `html/` — editar el de la raíz, no el de `html/`.

## Convenciones

- Notación: PVI como `y' = f(x,y)`, malla `x_n`, paso `h`; en EDP el segundo
  índice es `t_j = jk` (parabólicas) o `y_j = jh` (elípticas); iteraciones de
  Gauss–Seidel con superíndice `(m)`.
- Colores semánticos en ecuaciones y figuras: rojo `#d62728` = incógnita,
  azul `#1f77b4` = dato conocido.
- Terminología: "Euler adelante/atrás", "pendiente ponderada", "etapas" para
  los `k_i` de Runge–Kutta.
