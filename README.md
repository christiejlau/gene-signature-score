# Gene signature scoring
Calculates gene signature scores for single-cell RNA sequencing data, using methods described in Rehman et al _Cell_ 2021.

1. Read count matrix is z-score normalized gene-wise. 
2. From the defined gene set and associated weights (i.e. -1 for genes anticipated to be downregulated, +1 for upregulated genes), normalized read count values belonging to the defined gene set are multiplied by the associated weight. 
3. Weighted read count values belonging to the gene set are averaged to derive the final signature score. 
