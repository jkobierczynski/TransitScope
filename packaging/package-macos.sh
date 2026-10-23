#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# macOS release: TransitScope.app (desktop application, with the command-line program inside) in a
# .dmg, made self-contained with macdeployqt and signed ad hoc (not notarized).
# usage: packaging/package-macos.sh BUILD_DIR VERSION OUT_DIR
set -euo pipefail
B=$(cd "$1" && pwd); V=$2; mkdir -p "$3"; OUT=$(cd "$3" && pwd)
WORK=$(mktemp -d)
APP="$WORK/TransitScope.app"
cp -R "$B/transitscope-gui.app" "$APP"
cp "$B/transitscope" "$APP/Contents/MacOS/transitscope"
mkdir -p "$APP/Contents/Resources"
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cp -R "$ROOT/data" "$APP/Contents/Resources/data"
MACDEPLOYQT=${MACDEPLOYQT:-$(command -v macdeployqt || echo "$(brew --prefix qt)/bin/macdeployqt")}
"$MACDEPLOYQT" "$APP" -executable="$APP/Contents/MacOS/transitscope" -always-overwrite
codesign --force --deep --sign - "$APP"
codesign --verify --deep "$APP"

STAGE="$WORK/dmg"
mkdir -p "$STAGE"
mv "$APP" "$STAGE/"
ln -s /Applications "$STAGE/Applications"
cp "$ROOT/README.md" "$STAGE/README.md"
cp -R "$ROOT/examples" "$STAGE/examples"
ARCH=$(uname -m)
hdiutil create -volname "TransitScope $V" -srcfolder "$STAGE" -ov -format UDZO "$OUT/TransitScope-$V-macos-$ARCH.dmg"
rm -rf "$WORK"
ls -la "$OUT"
