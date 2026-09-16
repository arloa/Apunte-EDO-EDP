#!/bin/sh
# Arma docs/ listo para GitHub Pages: HTML + figs reales (Pages no sigue
# el symlink html/figs). index.html queda en la raiz del sitio.
cd "$(dirname "$0")"
./build_html.sh
rm -rf docs
mkdir docs
cp html/*.html html/style.css docs/
mkdir docs/figs
cp figs/*.png docs/figs/
touch docs/.nojekyll
echo "docs/ listo para publicar"
