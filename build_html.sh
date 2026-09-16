#!/bin/sh
# Genera el HTML de todos los capitulos en html/
# OJO: --mathjax necesita la URL explicita; el default de pandoc en Debian
# apunta a /usr/share/javascript/mathjax (MathJax 2 local, no instalado).
cd "$(dirname "$0")"

MATHJAX='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'
NAVDIR=$(mktemp -d)
trap 'rm -rf "$NAVDIR"' EXIT

titulo() { sed -n 's/^# //p' "$1" | head -n 1; }

# Nav anterior / indice / siguiente por capitulo (ventana de 3 sobre 0*.md)
ant2=; ant=
for f in 0*.md END; do
    if [ -n "$ant" ]; then
        {
            printf '<nav class="chapnav">'
            if [ -n "$ant2" ]; then
                printf '<a class="prev" href="%s.html">&larr; %s</a>' \
                    "${ant2%.md}" "$(titulo "$ant2")"
            else
                printf '<span></span>'
            fi
            printf '<a class="index" href="index.html">&Iacute;ndice</a>'
            if [ "$f" != END ]; then
                printf '<a class="next" href="%s.html">%s &rarr;</a>' \
                    "${f%.md}" "$(titulo "$f")"
            else
                printf '<span></span>'
            fi
            printf '</nav>'
        } > "$NAVDIR/${ant%.md}.html"
    fi
    ant2=$ant; ant=$f
done

cp style.css html/style.css

for f in 0*.md; do
    pandoc -s --mathjax="$MATHJAX" \
        --toc --toc-depth=2 \
        -M lang=es \
        -V pagetitle="$(titulo "$f")" \
        -c style.css \
        -B "$NAVDIR/${f%.md}.html" \
        -A "$NAVDIR/${f%.md}.html" \
        "$f" -o "html/${f%.md}.html"
done

# index.html desde el README (reescribe los links .md -> .html)
pandoc -s --mathjax="$MATHJAX" \
    -M lang=es \
    -M title="$(titulo README.md)" \
    -c style.css \
    README.md | sed 's/\.md"/.html"/g' > html/index.html
