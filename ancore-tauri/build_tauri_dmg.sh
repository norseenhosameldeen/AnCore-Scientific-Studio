#!/bin/bash
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR/src-tauri"
echo "Compiling Native macOS Rust/Tauri Application & DMG..."
cargo tauri build
echo "Done! Installer bundle generated in src-tauri/target/release/bundle/dmg/"
