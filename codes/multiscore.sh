#!/usr/bin/env bash
#PGS_IDS=(PGS002724 PGS000013 PGS002237 PGS001782 PGS004034)
PGS_IDS=(PGS002802 PGS000027)

set -euo pipefail

SCORES_DIR="scores"
mkdir -p logs "$SCORES_DIR"

for id in "${PGS_IDS[@]}"; do
  sb="pgs_${id}.sbatch"

  cat > "$sb" <<EOF
#!/bin/bash
#SBATCH -J pgs_${id}
#SBATCH -o logs/pgs_${id}.%j.out
#SBATCH -e logs/pgs_${id}.%j.err
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH -t 02:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=haxiao@unmc.edu

set -euo pipefail
umask 0002
cd "\$SLURM_SUBMIT_DIR"

PGS_ID="${id}"
BASE_DIR="${SCORES_DIR}/\${PGS_ID}"
RUN_DIR="\${BASE_DIR}/run_\${SLURM_JOB_ID}"
FINAL_DIR="${SCORES_DIR}/all_scores"
mkdir -p "\${RUN_DIR}" "\${RUN_DIR}/temp" "\${FINAL_DIR}"

shopt -s nullglob
FILES=(ALL.chr*.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz)
if (( \${#FILES[@]} == 0 )); then
  echo "No VCF：ALL.chr*.genotypes.vcf.gz" >&2
  exit 2
fi
shopt -u nullglob

conda run -n new ./pgs-calc apply --ref "\${PGS_ID}" \
  \${FILES[@]} \
  --genotypes GT \
  --out "\${RUN_DIR}/\${PGS_ID}_scores.txt"

mv "\${RUN_DIR}/\${PGS_ID}_scores.txt" "\${FINAL_DIR}/\${PGS_ID}_scores.txt"

echo "Finished：\${PGS_ID} → \${FINAL_DIR}/\${PGS_ID}_scores.txt"
EOF

  jid=$(sbatch --parsable "$sb")
  echo "Submitted ${id} → JobID=${jid}"
done
