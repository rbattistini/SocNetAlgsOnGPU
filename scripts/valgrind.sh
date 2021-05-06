#!/bin/sh

: <<'DOC'
Richiede compilazione con simboli di debug abilitati (flag -g in gcc)
DOC

valgrind  --tool=memcheck \
          --leak-check=yes \
          --show-reachable=yes \
          --num-callers=20 \
          --track-fds=yes \
          "$@"
