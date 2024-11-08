unsetopt AUTO_MENU       # Disable automatic menu on ambiguous completions
unsetopt MENU_COMPLETE   # Don’t use menu completion
unsetopt AUTO_LIST       # Don’t automatically list choices on ambiguous completion
unsetopt COMPLETE_IN_WORD # Disable completion within a word
bindkey "^I" self-insert
# Bind Shift + Tab to delete the last character (acts like "undoing" the tab)
bindkey "^[[Z" backward-delete-char

