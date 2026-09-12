#!/bin/bash
# ==========================================================
# SoMoSeq preprocessing - step 4: run zUMIs (filter, align, count)
# ==========================================================
#SBATCH -t 10-00:00:00
#SBATCH -n 8
#SBATCH --mem=400G

# usage: sbatch 04_run_zumis.sh <path-to-zUMIs-yaml>
# Runs zUMIs against the given config (see configs/zUMIs_template.run.yaml and paths.example.yaml). 
# Requires ZUMIS_INSTALL_DIR and R_LIBS_USER to be set for your environment.

ZUMIS_YAML=$1
if [[ -z "$ZUMIS_YAML" ]]; then
    echo "Usage: $0 <path-to-zUMIs-yaml>"
    exit 1
fi

echo 'R_LIBS_USER="${R_LIBS_USER}"' > $HOME/.Renviron

${ZUMIS_INSTALL_DIR}/zUMIs.sh -c -y $ZUMIS_YAML
