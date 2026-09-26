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
#SBATCH --account=iprime
#SBATCH --array=1-1

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
# -resume
```
