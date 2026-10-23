#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Linux release: an AppImage of the desktop application, and a tarball with both programs and the
# libraries they need (bin/, lib/, plugins/, share/), made with linuxdeploy.
# usage: packaging/package-linux.sh BUILD_DIR VERSION OUT_DIR
set -euo pipefail
B=$(cd "$1" && pwd); V=$2; mkdir -p "$3"; OUT=$(cd "$3" && pwd)
WORK=$(mktemp -d)
TOOLS=${LINUXDEPLOY_TOOLS:-$WORK/tools}
mkdir -p "$TOOLS"
for t in linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage \
         linuxdeploy/linuxdeploy-plugin-qt/releases/download/continuous/linuxdeploy-plugin-qt-x86_64.AppImage; do
  f="$TOOLS/$(basename "$t")"
  [[ -x $f ]] || { curl -fsSL -o "$f" "https://github.com/$t"; chmod +x "$f"; }
done
export PATH="$TOOLS:$PATH" APPIMAGE_EXTRACT_AND_RUN=1 LINUXDEPLOY_OUTPUT_VERSION="$V"
export QMAKE=${QMAKE:-$(command -v qmake6 || echo /usr/lib/qt6/bin/qmake6)}
export EXTRA_PLATFORM_PLUGINS="libqoffscreen.so"   # headless use (--screenshot, CI)

APPDIR="$WORK/AppDir"
DESTDIR="$APPDIR" cmake --install "$B" --prefix /usr > /dev/null
(cd "$WORK" && linuxdeploy-x86_64.AppImage --appdir "$APPDIR" \
   --executable "$APPDIR/usr/bin/transitscope-gui" --executable "$APPDIR/usr/bin/transitscope" \
   --desktop-file "$APPDIR/usr/share/applications/transitscope.desktop" \
   --icon-file "$APPDIR/usr/share/icons/hicolor/256x256/apps/transitscope.png" \
   --plugin qt --output appimage)
mv "$WORK"/TransitScope-*.AppImage "$OUT/TransitScope-$V-x86_64.AppImage"

NAME="transitscope-$V-linux-x86_64"
cp -a "$APPDIR/usr" "$WORK/$NAME"
cat > "$WORK/$NAME/README.txt" <<TXT
TransitScope $V for Linux x86-64

  bin/transitscope        command line (transitscope project new myproject.toml, transitscope run ...)
  bin/transitscope-gui    desktop application

The libraries both programs need are in lib/; nothing has to be installed. Data files (limb
darkening tables, passbands) are in share/transitscope/data. Documentation and an example
project: share/doc/transitscope.
TXT
tar -C "$WORK" -czf "$OUT/$NAME.tar.gz" "$NAME"
rm -rf "$WORK"
ls -la "$OUT"
