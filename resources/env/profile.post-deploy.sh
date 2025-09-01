#!/usr/bin/env bash
set -o pipefail

python_version="$(python -V | awk '{print $2}' | cut -d '.' -f 1-2)"
humann_file="${CONDA_PREFIX}/lib/python${python_version}/site-packages/humann/search/prescreen.py"

# This edit makes it so humann works with metaphlan profiles that have a single-column
sed -i \
    -e '92c\
        data=line.rstrip().split("\\t")' \
        "${humann_file}"
