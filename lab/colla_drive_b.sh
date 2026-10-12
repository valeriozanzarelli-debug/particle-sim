#!/bin/bash
# driver B (second version, 22:25): waits for the running ghs95_B2 chunk, then the GHS scan in priority order.
cd "$(dirname "$0")"
while tasklist 2>/dev/null | grep -q " 13004 "; do sleep 10; done
COLLA_THREADS=1 ./colla_drive.sh ghs95_B2 ghs95_B4 ghs80_B1 ghs80_B2 ghs80_B4 ghs95_B3 ghs80_B3 ghs0_B1 ghs0_B2 ghs0_B4 ghs50_B1 ghs50_B2 ghs50_B4 ghs95_B2pair ghs95_B4tet
