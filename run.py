"""
Perform diapause signature scoring on scRNA-seq data
"""

import os
import sys
import yaml

import anndata as ad
import pandas as pd
import numpy as np


def signature_score(
    adata: ad.AnnData, geneset_path: os.PathLike, geneset_id: str, inplace: bool = True
) -> None:
    """
    Calculate gene signature scores from average expression of weighted genes.

    Args:
        adata (ad.AnnData): AnnData object with scRNA-seq read counts.
        geneset_path (os.PathLike): Path to geneset CSV file.
        geneset_id (str): Gene set identifier (e.g. diapause).
        inplace (bool, optional): Whether signature scores will be amended to
        the AnnData .obs layer or returned as a DataFrame. Defaults to True.

    Returns:
        None: By default, AnnData .obs layer is automatically updated.
    """

    try:
        geneset = pd.read_csv(geneset_path)
    except AttributeError:
        "Cannot open path to geneset file"

    # Z-score scaling allows for read counts to be centered around the mean, so that
    # final scores can be maximized regardless of expression directionality.
    try:
        raw = adata.X.todense()
    except AttributeError:
        "Cannot find raw layer in AnnData object"
    stdev = np.std(raw, axis=0)
    znorm = raw - np.mean(raw, axis=0) / stdev
    reads = (
        pd.DataFrame(znorm, index=adata.obs_names, columns=adata.var_names)
        .fillna(0)
        .replace([np.inf, -np.inf], 0)
    )

    # Weights correspond to upregulated/downregulated genes from the gene set.
    # Downregulated genes are multiplied by -1 to ensure that the final score
    # represents the magnitude of the gene signature.
    gene_up = [gene for gene in geneset[geneset["weight"] >= 0]["gene"]]
    gene_down = [gene for gene in geneset[geneset["weight"] < 0]["gene"]]
    adata_up = reads[reads.columns.intersection(gene_up)]
    adata_down = reads[reads.columns.intersection(gene_down)] * -1
    adata_all = pd.concat([adata_up, adata_down], axis=1)
    gene_signature_score = np.round(adata_all.mean(axis=1), 3)
    if inplace == True:
        adata.obs.loc[:, f"{geneset_id}_score"] = gene_signature_score
    else:
        return gene_signature_score


def main():
    with open("./config.yaml", "r") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    adata = ad.read_h5ad(config['adata_path'])
    signature_score(adata, geneset_path=config['geneset_path'], geneset_id="diapause")
    adata.obs.to_csv(config['output_path'])


if __name__ == "__main__":
    sys.exit(main())

