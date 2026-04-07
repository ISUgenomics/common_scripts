# Keep normal bash history behavior
shopt -s histappend

__dirhist_last_histnum=0
__dirhist_cmd_pwd=

# Capture the directory where the next top-level command starts
__dirhist_capture_pwd() {
  [[ $- == *i* ]] || return
  [[ -n ${__dirhist_in_postexec:-} ]] && return
  [[ ${BASH_SUBSHELL:-0} -eq 0 ]] || return
  [[ ${#FUNCNAME[@]} -eq 1 ]] || return
  [[ -n ${__dirhist_cmd_pwd:-} ]] && return

  __dirhist_cmd_pwd=$PWD
}

__dirhist_postexec() {
  local exit_code=$?
  local hist_line hist_num last_cmd histfile

  [[ $- == *i* ]] || return

  __dirhist_in_postexec=1
  history -a
  hist_line=$(HISTTIMEFORMAT= history 1)
  unset __dirhist_in_postexec

  [[ -n $hist_line ]] || {
    unset __dirhist_cmd_pwd
    return
  }

  hist_num=${hist_line%%[[:space:]]*}
  [[ $hist_num =~ ^[0-9]+$ ]] || {
    unset __dirhist_cmd_pwd
    return
  }

  (( hist_num > __dirhist_last_histnum )) || {
    unset __dirhist_cmd_pwd
    return
  }
  __dirhist_last_histnum=$hist_num

  last_cmd=${hist_line#*[0-9]}
  last_cmd=${last_cmd#"${last_cmd%%[![:space:]]*}"}
  [[ -n $last_cmd ]] || {
    unset __dirhist_cmd_pwd
    return
  }

  histfile="${__dirhist_cmd_pwd:-$PWD}/.local_bash_history"
  printf '[%s] [exit:%s] %s\n' \
    "$(date '+%F %T')" \
    "$exit_code" \
    "$last_cmd" >> "$histfile" 2>/dev/null

  unset __dirhist_cmd_pwd
}

trap '__dirhist_capture_pwd' DEBUG
PROMPT_COMMAND="__dirhist_postexec${PROMPT_COMMAND:+;$PROMPT_COMMAND}"

