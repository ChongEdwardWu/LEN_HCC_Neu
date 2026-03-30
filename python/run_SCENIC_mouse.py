#!/usr/bin/env python3
import os
import subprocess
import sys
import logging
import pandas as pd
import loompy as lp
from pyscenic.binarization import binarize

# Soft limits to avoid over-threading inside BLAS/MKL.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

def ensure_dir(dir_path):
    """Ensure that a directory exists; if not, create it."""
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def run_command(cmd, step_name):
    """Run a shell command using subprocess and handle errors."""
    logging.info(f"Running {step_name}: {cmd}")
    try:
        result = subprocess.run(
            cmd, shell=True, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True
        )
        logging.info(f"{step_name} completed successfully:\n{result.stdout}")
        if result.stderr.strip():
            logging.info(f"{step_name} stderr:\n{result.stderr}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in {step_name}:\n{e.stderr}")
        sys.exit(1)

def main():
    # Configure logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    # ===== User-configurable =====
    DATASET_ID = "01_pySCENIC_counts_merged"   # without the .loom suffix
    WORKDIR = "path_to_project_root"
    AUXILIARIES_DIR = "path_to_scenic_reference"  # mouse or human SCENIC resources
    N_WORKERS_GRN = 32
    N_WORKERS_OTHERS = 4
    # =============================

    # Paths
    R_FOLDERNAME = os.path.join(WORKDIR, "04_R", "results")
    RESULTS_FOLDERNAME = os.path.join(WORKDIR, "05_pySCENIC", "results")
    FIGURES_FOLDERNAME = os.path.join(WORKDIR, "05_pySCENIC", "figures")
    ensure_dir(RESULTS_FOLDERNAME)
    ensure_dir(FIGURES_FOLDERNAME)

    INPUTLOOM_FNAME = os.path.join(R_FOLDERNAME, f"{DATASET_ID}.loom")

    # Ranking databases in a fixed order.
    RANKING_DBS = [
        os.path.join(AUXILIARIES_DIR, "mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"),
        os.path.join(AUXILIARIES_DIR, "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"),
    ]
    RANKING_DBS_FNAMES = " ".join(RANKING_DBS)

    # Auto-detect the motif annotation file and TF list from the resource folder.
    motif_files = sorted([f for f in os.listdir(AUXILIARIES_DIR)
                          if f.startswith("motifs") and f.endswith(".tbl")])
    if not motif_files:
        raise FileNotFoundError(f"No motif annotation file (*.tbl) found in {AUXILIARIES_DIR}")
    MOTIF_ANNOTATIONS_FNAME = os.path.join(AUXILIARIES_DIR, motif_files[0])
    logging.info(f"Motif annotations: {MOTIF_ANNOTATIONS_FNAME}")

    tfs_files = sorted([f for f in os.listdir(AUXILIARIES_DIR)
                        if f.startswith("allTFs") and f.endswith(".txt")])
    if not tfs_files:
        raise FileNotFoundError(f"No TF file (allTFs*.txt) found in {AUXILIARIES_DIR}")
    MM_TFS_FNAME = os.path.join(AUXILIARIES_DIR, tfs_files[0])
    logging.info(f"TF list: {MM_TFS_FNAME}")

    # Outputs
    ADJACENCIES_FNAME = os.path.join(RESULTS_FOLDERNAME, f"{DATASET_ID}.adjacencies.tsv")
    MOTIFS_FNAME = os.path.join(RESULTS_FOLDERNAME, f"{DATASET_ID}.motifs.csv")
    AUCELL_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, f"{DATASET_ID}.auc.csv")
    BIN_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, f"{DATASET_ID}.bin.csv")
    THR_FNAME = os.path.join(RESULTS_FOLDERNAME, f"{DATASET_ID}.thresholds.csv")
    LOOM_FNAME = os.path.join(RESULTS_FOLDERNAME, f"{DATASET_ID}.scenic.loom")

    # ---- Preflight checks ----
    for p in [INPUTLOOM_FNAME, MOTIF_ANNOTATIONS_FNAME, MM_TFS_FNAME] + RANKING_DBS:
        if not os.path.exists(p):
            logging.error(f"Missing file: {p}")
            sys.exit(1)

    # Quick sanity check on gene naming and TF overlap.
    try:
        with lp.connect(INPUTLOOM_FNAME, mode="r") as ds:
            # Gene row attribute keys may be stored as 'Gene' or 'gene'.
            if "Gene" in ds.ra:
                genes = pd.Index(ds.ra["Gene"][:])
            elif "gene" in ds.ra:
                genes = pd.Index(ds.ra["gene"][:])
            else:
                logging.warning("Gene names not found in loom row attributes (Gene/gene). Skipping gene-format check.")
                genes = pd.Index([])
    except Exception as e:
        logging.warning(f"Could not open loom for gene check: {e}")
        genes = pd.Index([])

    if len(genes) > 0:
        frac_ensembl = (genes.str.upper().str.startswith("ENSMUSG")).mean()
        if frac_ensembl > 0.05:
            logging.error(
                f"Your loom seems to use ENSEMBL IDs ({frac_ensembl:.1%}). "
                f"Current ranking DBs & motif annotations are MGI-based. Please convert to MGI symbols."
            )
            sys.exit(1)
        try:
            tfs = pd.read_csv(MM_TFS_FNAME, header=None)[0].astype(str).str.strip()
            overlap = len(set(genes).intersection(set(tfs)))
            if overlap < 100:
                logging.warning(f"Only {overlap} TFs found in expression genes. "
                                f"Check symbol convention and TF list.")
        except Exception as e:
            logging.warning(f"TF overlap check skipped: {e}")

    # ===== STEP 1: GRN inference =====
    # Explicitly set the method and random seed for reproducibility.
    cmd1 = (
        f"pyscenic grn {INPUTLOOM_FNAME} {MM_TFS_FNAME} "
        f"-o {ADJACENCIES_FNAME} --num_workers {N_WORKERS_GRN} "
        f"--method grnboost2 --seed 123"
    )
    run_command(cmd1, "SCENIC STEP 1: GRN inference")

    # ===== STEP 2: cisTarget / regulon prediction =====
    cmd2 = (
        f"pyscenic ctx {ADJACENCIES_FNAME} {RANKING_DBS_FNAMES} "
        f"--annotations_fname {MOTIF_ANNOTATIONS_FNAME} "
        f"--expression_mtx_fname {INPUTLOOM_FNAME} "
        f"--output {MOTIFS_FNAME} --auc_threshold 0.05 "
        f"--num_workers {N_WORKERS_OTHERS}"
    )
    run_command(cmd2, "SCENIC STEP 2: Regulon prediction")

    # ===== STEP 3a: AUCell -> CSV (used for binarization) =====
    cmd3a = (
        f"pyscenic aucell {INPUTLOOM_FNAME} {MOTIFS_FNAME} "
        f"--output {AUCELL_MTX_FNAME} --num_workers {N_WORKERS_OTHERS}"
    )
    run_command(cmd3a, "SCENIC STEP 3a: AUCell -> CSV")

    # ===== STEP 3b: AUCell -> LOOM (for visualization and downstream analysis) =====
    cmd3b = (
        f"pyscenic aucell {INPUTLOOM_FNAME} {MOTIFS_FNAME} "
        f"--output {LOOM_FNAME} --num_workers {N_WORKERS_OTHERS}"
    )
    run_command(cmd3b, "SCENIC STEP 3b: AUCell -> LOOM")

    # ===== STEP 4: Regulon activity binarization =====
    try:
        auc_mtx = pd.read_csv(AUCELL_MTX_FNAME, index_col=0)
    except Exception as e:
        logging.error(f"Failed to read AUC matrix from {AUCELL_MTX_FNAME}: {e}")
        sys.exit(1)

    bin_mtx, thresholds = binarize(auc_mtx, seed=123, num_workers=N_WORKERS_OTHERS)
    bin_mtx.to_csv(BIN_MTX_FNAME)
    thresholds.to_frame().rename(columns={0: "threshold"}).to_csv(THR_FNAME)
    logging.info("SCENIC STEP 4: Regulon activity binarization is done!")

if __name__ == "__main__":
    main()
