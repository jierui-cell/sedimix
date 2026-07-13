#!/usr/bin/env bash
set -euo pipefail

kraken2_dir="${KRAKEN2_DIR:-/opt/kraken2}"
export PATH="${kraken2_dir}:${PATH}"

for executable in \
  kraken2 kraken2-build kraken2-inspect \
  classify build_db k2 estimate_capacity; do
  test -x "${kraken2_dir}/${executable}"
done

kraken2 --version 2>&1 | grep -q "Kraken version 2.1.3"

workdir="$(mktemp -d)"
trap 'rm -rf "${workdir}"' EXIT
mkdir -p "${workdir}/db/taxonomy"

cat > "${workdir}/db/taxonomy/nodes.dmp" <<'EOF'
1	|	1	|	no rank	|
9606	|	1	|	species	|
EOF
cat > "${workdir}/db/taxonomy/names.dmp" <<'EOF'
1	|	root	|		|	scientific name	|
9606	|	Homo sapiens	|		|	scientific name	|
EOF
python3 - "${workdir}/library.fa" <<'PY'
import random
import sys

random.seed(42)
sequence = "".join(random.choice("ACGT") for _ in range(500))
with open(sys.argv[1], "w") as fasta:
    fasta.write(">test_sequence|kraken:taxid|9606\n")
    fasta.write(sequence + "\n")
PY

kraken2-build --add-to-library "${workdir}/library.fa" --db "${workdir}/db"
kraken2-build --build --db "${workdir}/db" --threads 1
kraken2 \
  --db "${workdir}/db" \
  --output "${workdir}/classification.tsv" \
  --report "${workdir}/report.tsv" \
  "${workdir}/library.fa"

grep -q $'^C\t' "${workdir}/classification.tsv"
grep -q $'\t9606\t' "${workdir}/classification.tsv"

echo "Kraken2 functional smoke test passed"
