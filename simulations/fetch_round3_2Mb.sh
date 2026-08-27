#!/bin/bash
# Mirror the O2 archive into local simulation_data/.
#
#   bash simulations/fetch_round3_2Mb.sh              # report + mirror
#   DRY=1 bash simulations/fetch_round3_2Mb.sh        # report only
#
# TWO O2 CONNECTIONS, NOT ONE PER REPLICATE
#
# This used to loop over ~130 replicates calling `ssh` and `rsync -e ssh` once
# each. Every one of those is an unmultiplexed connection to O2, and O2 re-runs
# Duo for any connection that does not reuse an existing master -- so a single
# run could fire a hundred-odd 2FA pushes at the user. That is the whole reason
# the o2-ssh / o2-rsync wrappers exist: they hold a lock while opening ONE
# shared connection, and everything afterwards rides on it.
#
# So: one `o2-ssh` to list the archive, one `o2-rsync` to mirror it. The archive
# already has the exact layout simulation_data/ uses (<arm>/<replicate>/, flat),
# which is what makes a single recursive copy sufficient.
#
# Run archive_round3_to_data2.sh on O2 first; this only mirrors what is there,
# and that script only puts a replicate there once all four stages have landed.
#
# NEVER add a bare `ssh`, `scp`, or `rsync -e ssh` to this file.
#
# AND DO NOT CALL `o2-ssh --check` HERE. It looks like a harmless read-only
# probe and is not: --check runs the DEEP path, which bounds a no-op at 10
# seconds and, on timeout, runs `ssh -O exit` and tears the master down
# (~/.local/bin/o2-ssh:68-78). A loaded O2 login node answers a no-op in more
# than 10 seconds often enough that a pre-flight --check regularly destroys the
# connection the user just authenticated -- costing the very Duo push it was
# meant to save. Plain `o2-ssh '<cmd>'` uses the cheap socket check instead and
# opens the master on demand, so just issue the command.

set -uo pipefail

ARCHIVE="${ARCHIVE:-/n/data2/hms/dbmi/sunyaev/lab/nconnally/new_simulation_set}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCAL="$(cd "$HERE/.." && pwd)/simulation_data"
DRY="${DRY:-}"

command -v o2-ssh >/dev/null || { echo "ERROR: o2-ssh not on PATH -- see the remote-o2 skill." >&2; exit 1; }

echo "Listing the archive (1 connection) ..."
# No pre-flight check: o2-ssh opens the master on demand via its cheap path, and
# probing first can tear down a healthy connection (see the header).
REMOTE=$(o2-ssh "cd '$ARCHIVE' && for a in */; do a=\${a%/}; for r in \"\$a\"/*/; do echo \"\$a/\$(basename \$r)\"; done; done") \
    || { echo "ERROR: could not list $ARCHIVE" >&2; exit 1; }

n_remote=$(printf '%s\n' "$REMOTE" | grep -c .)
missing=""; present=0
while read -r rel; do
    [ -n "$rel" ] || continue
    if compgen -G "$LOCAL/$rel/*.enloc.sig.out" > /dev/null 2>&1; then
        present=$((present + 1))
    else
        missing="$missing $rel"
    fi
done <<< "$REMOTE"

echo
echo "archive: $n_remote replicates   local: $present   missing locally:$(
    [ -n "$missing" ] && echo "$missing" || echo ' none')"

if [ -n "$DRY" ]; then
    echo
    echo "DRY=1 -- nothing transferred."
    exit 0
fi

echo
echo "Mirroring (1 connection, via the transfer node) ..."
# -rlt not -a: -a implies -pgo, and preserving the archive's setgid group
# permissions onto local APFS fails with "fchmodat: Operation not permitted",
# which makes rsync exit non-zero on transfers that actually succeeded.
# transfer: rather than the login node -- this can move gigabytes.
# NO --info=stats2 and no --stats: macOS ships openrsync (protocol 29), which
# has neither, and an unknown option there is fatal -- it prints its whole usage
# block and exits non-zero, the same class of wart as the --ignore-missing-args
# failure the header describes. The transfer then only happens on the retry.
o2-rsync -rlt "transfer:$ARCHIVE/" "$LOCAL/" \
    || o2-rsync -rlt "transfer:$ARCHIVE/" "$LOCAL/"

echo
still=0
while read -r rel; do
    [ -n "$rel" ] || continue
    compgen -G "$LOCAL/$rel/*.enloc.sig.out" > /dev/null 2>&1 || { echo "  STILL MISSING: $rel"; still=$((still + 1)); }
done <<< "$REMOTE"
echo "================================================================"
echo "archive=$n_remote replicates   local now=$((n_remote - still))   still missing=$still"
exit 0
