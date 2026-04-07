# Keep normal bash history behavior
shopt -s histappend

PROMPT_COMMAND="__dirhist_postexec${PROMPT_COMMAND:+;$PROMPT_COMMAND}"
PROMPT_COMMAND="history -a${PROMPT_COMMAND:+;$PROMPT_COMMAND}"

# Track only the next top-level interactive command
__dirhist_preexec() {
  # Only in interactive shells
  [[ $- == *i* ]] || return

  # Ignore nested DEBUG traps from functions, subshells, etc.
  [[ ${BASH_SUBSHELL:-0} -eq 0 ]] || return
  [[ ${#FUNCNAME[@]} -eq 1 ]] || return

  # Record once per prompt cycle
  [[ -n ${__dirhist_logged:-} ]] && return

  __dirhist_cmd=$BASH_COMMAND
  __dirhist_pwd=$PWD
  __dirhist_logged=1
}

# Append the command to a file in the directory it was run from
__dirhist_postexec() {
  local status=$?
  local histfile

  [[ -n ${__dirhist_cmd:-} ]] || return

  histfile="$__dirhist_pwd/.local_bash_history"
  printf '[%(%F %T)T] [exit:%s] %s\n' -1 "$status" "$__dirhist_cmd" >> "$histfile" 2>/dev/null

  unset __dirhist_cmd
  unset __dirhist_pwd
  unset __dirhist_logged
}

trap '__dirhist_preexec' DEBUG

