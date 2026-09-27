# UVA
---
```
#!/bin/bash
#SBATCH --job-name=nextflow_bulknaseq_star_rsem
#SBATCH --nodes=1
#SBATCH --output=autput_%x_%A_%a.out
#SBATCH --partition=standard
#SBATCH --ntasks=24
#SBATCH --mem=256gb
#SBATCH --time=72:00:00
#SBATCH --account=my_account
#SBATCH --array=1-1

APPTAINER_STORAGE=<where_apptainer_image_stored>
export NXF_SINGULARITY_CACHEDIR=$APPTAINER_STORAGE
export APPTAINER_CACHEDIR=$APPTAINER_STORAGE

module load nextflow
module load apptainer

nextflow run nf-core/rnaseq \
-profile apptainer \
-w ./work_bulkrna \
--input samplesheet_bulk_bulk.csv \
--fasta GRCh38.primary_assembly.genome.fa \
--gtf gencode.v49.primary_assembly.annotation.gtf \
--gencode \
--trimmer fastp \
--aligner star_rsem \
--outdir ./results_dlfpc_bulk_bulk_rsem \
--skip_bigwig \
--skip_stringtie \
--skip_dupradar \
--skip_qualimap \
--skip_rseqc \
--skip_preseq \
--skip_markduplicates
-resume
```
# M

```
#!/bin/bash
#SBATCH --job-name=nf_rna_star_rsem_avatar5_batch
#SBATCH --nodes=1
#SBATCH --output=output_%x_%A_%a.out
#SBATCH --ntasks=32
#SBATCH --mem=128gb
#SBATCH --time=72:00:00
#SBATCH --array=1-10
#SBATCH --qos=small

BATCH_DIR="batch_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$BATCH_DIR"

SAMPLESHEET=/home/myaccount/RNA_pipelines/nextflow_avatar5_simple_action1/action1_samplesheet.csv
HEADER=$(head -1 "$SAMPLESHEET")

N=4
START=$(( (SLURM_ARRAY_TASK_ID - 1) * N + 2 ))
END=$(( START + N - 1 ))

{ echo "$HEADER"; sed -n "${START},${END}p" "$SAMPLESHEET"; } > "$BATCH_DIR/samplesheet_batch${SLURM_ARRAY_TASK_ID}.csv"

sed "s|^input:.*|input: samplesheet_batch${SLURM_ARRAY_TASK_ID}.csv|" params.yaml > "$BATCH_DIR/params.yaml"

cd "$BATCH_DIR"


# 1. Properly initialize conda for non-interactive shell environments
eval "$(conda shell.bash hook)"
conda activate bulkRNA

# 2. CLEAR the module's default JAVA_TOOL_OPTIONS
unset JAVA_TOOL_OPTIONS

# 3. SET your working paths and custom truststore configuration
export JAVA_TOOL_OPTIONS="-Djavax.net.ssl.trustStore=/home/myaccount/cacerts -Djavax.net.ssl.trustStorePassword=changeit"

# 4. FORCE Nextflow to use your custom truststore and ignore certificate blocks / network lookups
export NXF_OPTS="-Djavax.net.ssl.trustStore=/home/myaccount/cacerts -Djavax.net.ssl.trustStorePassword=changeit -Dnextflow.disable.http.certificates.check=true"
export NXF_OFFLINE=true

# 5. Apptainer specifics
export APPTAINER_CACHEDIR=/home/myaccount/apptainer_cache
export APPTAINER_CACERT=/home/myaccount/SectigoPublicServerAuthenticationRootE46.crt
export NXF_APPTAINER_CACHEDIR=/home/myaccount/apptainer_cache

6. Run Nextflow pipeline
nextflow run nf-core/rnaseq \
    -params-file params.yaml \
    -profile apptainer \
    -offline \
    -w "$TMPDIR/work" \
    -resume

find . -type f -name "*.bam" -delete
rm -rf $TMPDIR/work
```
