setopt APPEND_HISTORY
setopt INC_APPEND_HISTORY

__dirhist_preexec() {
  [[ -o interactive ]] || return
  __dirhist_cmd="$1"
  __dirhist_pwd="$PWD"
}

__dirhist_precmd() {
  local exit_code=$?
  local histfile

  [[ -n "${__dirhist_cmd:-}" ]] || return

  histfile="$__dirhist_pwd/.local_zsh_history"
  print -r -- "[$(date '+%F %T')] [exit:$exit_code] $__dirhist_cmd" >> "$histfile" 2>/dev/null

  unset __dirhist_cmd
  unset __dirhist_pwd
}

autoload -Uz add-zsh-hook
add-zsh-hook preexec __dirhist_preexec
add-zsh-hook precmd __dirhist_precmd

