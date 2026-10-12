#!/bin/bash
# drivers B'' and C (23:14): split the remaining GHS scan into two disjoint queues (1 thread each).
# usage: colla_drive_bc.sh WAITPID JOB...   (waits for the Windows process WAITPID to end first)
cd "$(dirname "$0")"
w=$1; shift
while tasklist 2>/dev/null | grep -q " $w "; do sleep 10; done
COLLA_THREADS=1 ./colla_drive.sh "$@"
