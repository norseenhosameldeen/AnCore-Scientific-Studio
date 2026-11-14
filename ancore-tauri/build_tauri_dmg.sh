#!/bin/bash
set -e
echo "================================================================================"
echo " Building AnCore Studio Native Rust/Tauri macOS Application & DMG..."
echo "================================================================================"
export PATH="$HOME/.cargo/bin:/opt/homebrew/bin:/usr/local/bin:$PATH"

cd "/Users/ahmedsenary/Desktop/AnCore_project/ancore-tauri/src-tauri"
echo "Compiling Rust binary & packaging 'AnCore Studio.dmg'..."
cargo tauri build

DMG_OUT=$(find "/Users/ahmedsenary/Desktop/AnCore_project/ancore-tauri/src-tauri/target/release/bundle/dmg" -name "*.dmg" | head -n 1)
if [ -n "$DMG_OUT" ]; then
    cp "$DMG_OUT" "$HOME/Desktop/AnCore Studio.dmg"
    cp "$DMG_OUT" "$HOME/Downloads/AnCore Studio.dmg"
    echo "================================================================================"
    echo " [100%] SUCCESS! Standalone Tauri Rust DMG copied to:"
    echo "   - $HOME/Desktop/AnCore Studio.dmg"
    echo "   - $HOME/Downloads/AnCore Studio.dmg"
    echo "================================================================================"
    open -R "$HOME/Desktop/AnCore Studio.dmg"
fi
