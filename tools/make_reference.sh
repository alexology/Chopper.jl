#!/usr/bin/env bash
# Regenerates test/reference/ with the original chopper (Linux binary of the
# release in UPSTREAM.toml), one thread so that the output order is defined:
#
#     bash tools/make_reference.sh /path/to/chopper-linux
#
# On Windows run it inside WSL, from the package folder (as /mnt/d/...).
# For every line of test/reference/cases.tsv it writes <name>.fastq (stdout),
# <name>.stderr and <name>.status (exit status).
set -u
chopper=${1:?usage: make_reference.sh /path/to/chopper}
here=$(cd "$(dirname "$0")/.." && pwd)
ref="$here/test/reference"
"$chopper" --version > "$ref/VERSION"
cd "$here/test/data"
tail -n +2 "$ref/cases.tsv" | while IFS=$'\t' read -r name input arguments; do
    # shellcheck disable=SC2086  # arguments are split on purpose
    "$chopper" -t 1 -i "$input" $arguments > "$ref/$name.fastq" 2> "$ref/$name.stderr"
    echo $? > "$ref/$name.status"
    printf '%-22s exit %s, %s records\n' "$name" "$(cat "$ref/$name.status")" \
        "$(grep -c '' "$ref/$name.fastq" | awk '{print $1 / 4}')"
done
