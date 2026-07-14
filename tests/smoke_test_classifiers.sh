#!/usr/bin/env bash
set -euo pipefail

export PATH="/opt/kraken2:/usr/local/bin:${PATH}"

workdir="$(mktemp -d)"
trap 'rm -rf "${workdir}"' EXIT
mkdir -p \
  "${workdir}/scripts" \
  "${workdir}/kraken2_db/taxonomy" \
  "${workdir}/kraken2_run/0_data" \
  "${workdir}/centrifuge_run/0_data"
cp /opt/sedimix/scripts/select_reads_from_kraken2.py "${workdir}/scripts/"
cp /opt/sedimix/scripts/select_reads_from_centrifuge.py "${workdir}/scripts/"

python3 - "${workdir}" <<'PY'
import random
import sys
from pathlib import Path

root = Path(sys.argv[1])
random.seed(42)
sequence = "".join(random.choice("ACGT") for _ in range(500))
quality = "I" * len(sequence)

(root / "kraken2_db/taxonomy/nodes.dmp").write_text(
    "1\t|\t1\t|\tno rank\t|\n9606\t|\t1\t|\tspecies\t|\n"
)
(root / "kraken2_db/taxonomy/names.dmp").write_text(
    "1\t|\troot\t|\t\t|\tscientific name\t|\n"
    "9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|\n"
)
(root / "kraken2_reference.fa").write_text(
    ">tiny_reference|kraken:taxid|9606\n" + sequence + "\n"
)
(root / "centrifuge_reference.fa").write_text(
    ">tiny_reference\n" + sequence + "\n"
)
(root / "seqid2taxid.map").write_text("tiny_reference\t9606\n")
(root / "nodes.dmp").write_text(
    "1\t|\t1\t|\tno rank\t|\n9606\t|\t1\t|\tspecies\t|\n"
)
(root / "names.dmp").write_text(
    "1\t|\troot\t|\t\t|\tscientific name\t|\n"
    "9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|\n"
)
fastq = "@test_read\n" + sequence + "\n+\n" + quality + "\n"
(root / "kraken2_run/0_data/sample.fq").write_text(fastq)
(root / "centrifuge_run/0_data/sample.fq").write_text(fastq)
PY

# Build the smallest useful databases locally. This avoids multi-gigabyte
# production downloads while exercising the real database builders and readers.
kraken2-build --add-to-library "${workdir}/kraken2_reference.fa" --db "${workdir}/kraken2_db"
kraken2-build --build --db "${workdir}/kraken2_db" --threads 1
centrifuge-build \
  --conversion-table "${workdir}/seqid2taxid.map" \
  --taxonomy-tree "${workdir}/nodes.dmp" \
  --name-table "${workdir}/names.dmp" \
  "${workdir}/centrifuge_reference.fa" \
  "${workdir}/centrifuge_index"

cat > "${workdir}/kraken2_run/config.yaml" <<EOF
classification_software: kraken2
classification_index: ${workdir}/kraken2_db
taxID: 9606
min_length: 35
threads: 1
memory_mb: 1024
ref_genome: unused.fa
use_snp_panel: false
snp_panel_bed: unused.bed
min_quality: 0
calculate_from_mapdamage: false
lineage_sites: unused.tsv
types: Neanderthal
keep_non_hominin_reads: false
to_clean: false
EOF

cat > "${workdir}/centrifuge_run/config.yaml" <<EOF
classification_software: centrifuge
classification_index: ${workdir}/centrifuge_index
taxID: 9606
min_length: 35
threads: 1
memory_mb: 1024
ref_genome: unused.fa
use_snp_panel: false
snp_panel_bed: unused.bed
min_quality: 0
calculate_from_mapdamage: false
lineage_sites: unused.tsv
types: Neanderthal
keep_non_hominin_reads: false
to_clean: false
EOF

(
  cd "${workdir}/kraken2_run"
  snakemake \
    --snakefile /opt/sedimix/rules/snakefile_sedimix \
    --cores 1 \
    temp/sample_length_filtered_classified.fq
)
(
  cd "${workdir}/centrifuge_run"
  snakemake \
    --snakefile /opt/sedimix/rules/snakefile_sedimix \
    --cores 1 \
    temp/sample_length_filtered_classified.fq
)

grep -q '^@test_read$' "${workdir}/kraken2_run/temp/sample_length_filtered_classified.fq"
grep -q $'^C\ttest_read\t9606\t' "${workdir}/kraken2_run/1_classification/sample.kraken2"
grep -q '^@test_read$' "${workdir}/centrifuge_run/temp/sample_length_filtered_classified.fq"
grep -q $'^test_read\t.*\t9606\t' "${workdir}/centrifuge_run/1_classification/sample.centrifuge"
test -s "${workdir}/centrifuge_run/1_classification/sample.k2report"

echo "Kraken2 and Centrifuge workflow smoke tests passed"
