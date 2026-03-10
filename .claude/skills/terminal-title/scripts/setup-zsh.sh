#!/bin/bash
# Setup zsh precmd hook to apply terminal titles set by Claude Code
# Usage: bash scripts/setup-zsh.sh

ZSHRC="$HOME/.zshrc"
TITLE_FILE="\$HOME/.claude/terminal_title"
MARKER="# claude-terminal-title"

# Check if already installed
if grep -q "$MARKER" "$ZSHRC" 2>/dev/null; then
    echo "Already installed in $ZSHRC"
    exit 0
fi

# Append precmd hook
cat >> "$ZSHRC" << 'EOF'

# claude-terminal-title: apply terminal title set by Claude Code
_claude_title_precmd() {
    if [[ -f "$HOME/.claude/terminal_title" ]]; then
        printf '\033]0;%s\007' "$(cat "$HOME/.claude/terminal_title")"
    fi
}
precmd_functions+=(_claude_title_precmd)
EOF

echo "Installed precmd hook in $ZSHRC"
echo "Restart your shell or run: source $ZSHRC"
