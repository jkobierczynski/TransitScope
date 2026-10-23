#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Windows release (run in an MSYS2 UCRT64 shell): a zip with both programs, the DLLs they need
# from MSYS2 (Qt, CFITSIO, the GCC runtime, ...), the Qt platform and style plugins, and the data.
# usage: packaging/package-windows.sh BUILD_DIR VERSION OUT_DIR
set -euo pipefail
B=$(cd "$1" && pwd); V=$2; mkdir -p "$3"; OUT=$(cd "$3" && pwd)
NAME="transitscope-$V-windows-x86_64"
WORK=$(mktemp -d)
D="$WORK/$NAME"
cmake --install "$B" --prefix "$D" > /dev/null
prefix=$(cygpath -u "$MINGW_PREFIX")
QWIN=$(find "$prefix" -path '*plugins/platforms/qwindows.dll' | head -1)
[[ -n $QWIN ]] || { echo "qwindows.dll not found under $prefix" >&2; exit 1; }
PLUGINS=$(dirname "$(dirname "$QWIN")")
mkdir -p "$D/bin/platforms"
cp "$PLUGINS/platforms/qwindows.dll" "$D/bin/platforms/"
if [[ -d $PLUGINS/styles ]]; then mkdir -p "$D/bin/styles" && cp "$PLUGINS"/styles/*.dll "$D/bin/styles/"; fi

# every DLL from the MSYS2 prefix that the programs and plugins load, repeated until complete
for pass in 1 2 3 4; do
  find "$D/bin" -name '*.exe' -o -name '*.dll' | while read -r f; do
    ldd "$f" | awk '{print $3}' | grep -i "^$prefix/bin/" || true
  done | sort -u | while read -r dll; do
    [[ -f "$D/bin/$(basename "$dll")" ]] || cp "$dll" "$D/bin/"
  done
done
cat > "$D/README.txt" <<TXT
TransitScope $V for Windows x86-64

  bin\transitscope-gui.exe   desktop application
  bin\transitscope.exe       command line (transitscope project new myproject.toml, transitscope run ...)

Nothing has to be installed; keep the folder together. Data files (limb darkening tables,
passbands) are in share\transitscope\data; documentation and an example project in share\doc\transitscope.
TXT
(cd "$WORK" && zip -qr "$OUT/$NAME.zip" "$NAME")
rm -rf "$WORK"
ls -la "$OUT"
