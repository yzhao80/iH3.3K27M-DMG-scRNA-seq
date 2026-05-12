# iH3.3K27M DMG scRNA-seq Analysis
### Khairkhah et al. — R analysis pipeline
### Author: Yue Zhao, yuedz@umich.edu

This repository contains all R scripts for reproducing the single-cell RNA-seq analysis in **Khairkhah et al.** (manuscript in preparation), studying the transcriptional consequences of inducible and reversible H3.3K27M oncohistone expression in a mouse diffuse midline glioma (DMG) model.

---

## Repository structure

```
.
├── 01_QC_doublet_removal.R                    # QC, DecontX, scDblFinder
├── 02_normalization_integration_clustering.R  # Seurat pipeline, RPCA, UMAP
├── 03_differential_expression.R               # FindMarkers ON vs OFF, heatmaps, violins
├── 04_cellchat_analysis.R                     # CellChat construction & core analysis
├── 05_cellchat_circle_plots.R                 # Circle plots per pathway
├── 06_cellchat_heatmap_pathways.R             # LR heatmaps & diverging bar charts
├── 07_immune_subclustering.R                  # Immune cell re-clustering, UMAP, barplot
├── utilities/
│   ├── rankNet_patched.R                      # Patched CellChat rankNet (see below)
│   ├── netAnalysis_signalingRole_heatmap2.R   # Extended signaling role heatmap
│   └── violin_plot_celltype.R                 # Annotated violin plot panels
├── data/                                      # Input data (not tracked – see below)
└── results/                                   # Output figures & tables (generated)
```

---

## Script → Figure mapping

| Script | Figure panels |
|--------|--------------|
| `01_QC_doublet_removal.R` | — (QC filtering) |
| `02_normalization_integration_clustering.R` | **Fig 5E** (annotated UMAP), **Fig 5F** (cell-type proportion bar) |
| `03_differential_expression.R` | **Sup Fig 7B** (H3f3a violin), **Sup Fig 7D** (excitatory neuron markers), **Sup Fig 7E** (DMG glutamate/AMPAR markers: Slc1a3, Glud1, Got1, Got2, Gria3, Gria4), **Sup Fig 7F** (stemness/diff heatmap), **Sup Fig 7G** (inhibitory neuron markers) |
| `04_cellchat_analysis.R` | **Fig 5B** (rankNet bar), **Fig 5D** (overall signaling heatmap), **Fig 6A** (circle – N interactions), **Fig 6B** (circle – interaction strength), **Fig 6G** (outgoing/incoming heatmaps), **Sup Fig 7H** (differential circles), **Sup Fig 7I** (total interaction bars) |
| `05_cellchat_circle_plots.R` | **Fig 5D** (per-pathway circles: NCAM, NRXN, PTN, CADM), **Fig 6** (EPHA, NT, OPIOID, SOMATOSTATIN circles) |
| `06_cellchat_heatmap_pathways.R` | **Fig 5C** (diverging bar – DMG signaling), **Fig 5G** (LR heatmap), **Fig 5H** (LR bubble plot), **Fig 6F** (OFF-specific LR bubble) |
| `07_immune_subclustering.R` | **Fig 6D** (immune cell UMAP, split ON/OFF), **Fig 6E** (immune subtype composition barplot) |

---

## Pipeline overview (Methods order)

```
Raw Cell Ranger output (mm10, FLEX)
       │
       ▼
01  DecontX (ambient RNA) → scDblFinder (doublets) → QC filter
       │
       ▼
02  LogNormalize → VST variable features → PCA (50 PCs)
    → RPCA integration (4 samples) → SNN graph → Louvain clustering (res=0.4)
    → UMAP (12 dims) → Manual cell-type annotation
       │
       ▼
03  FindAllMarkers (cluster markers, non-integrated RNA)
    FindMarkers per cell type: ON vs OFF (Wilcoxon, Bonferroni)
    → violin_plot_celltype.R: annotated violin panels
      • Excitatory neurons (Sup Fig 7D)
      • DMG tumor cells – glutamate/AMPAR markers (Sup Fig 7E)
      • Inhibitory neurons (Sup Fig 7G)
       │
       ├──────────────────────────────────────────────────────────┐
       ▼                                                          ▼
04  CellChat (v1.6.1, CellChatDB.mouse)                  07  Immune subclustering
    → Build ON / OFF objects separately                       → Subset "Immune Cells"
    → Communication probabilities, centrality scores,         → Re-normalize, re-cluster
      patterns                                                   (PCA 15 dims, res=0.3)
    → Merge for comparative analysis                          → Annotate 9 immune subtypes
    → rankNet_patched.R: information flow bar (Fig 5B)        → UMAP split ON/OFF (Fig 6D)
    → netAnalysis_signalingRole_heatmap2.R: role              → Composition barplot (Fig 6E)
      heatmaps (Fig 5D, 6G)                                   → Fisher's exact test
       │
    ┌──┴──────────────────┐
    ▼                     ▼
05  Circle plots       06  Pathway heatmaps & LR bubble plots
                           → rankNet_patched.R (data extraction)
                           → netAnalysis_signalingRole_heatmap2.R
```

---

## Immune cell subtypes (script 07)

The 9 immune subtypes identified after re-clustering cluster 18 (Louvain res = 0.3, 15 PCs):

| Short label | Full name | Key markers |
|-------------|-----------|-------------|
| ITRM | Immunosuppressive tissue-resident macrophages | Cd68, Mrc1, Tgfb1 |
| Tregs | Regulatory T cells | Foxp3, Ctla4, Il2ra |
| CD8+ T | CD8+ cytotoxic T cells | Cd8a, Gzma, Prf1 |
| DC | Dendritic cells | Itgax, Cd74, H2-Ab1 |
| Act. DC | Activated dendritic cells | Cd86, Cd80, Il12b |
| B cells | B cells | Cd79a, Ms4a1, Pax5 |
| pDC | Plasmacytoid dendritic cells | Siglech, Bst2, Irf7 |
| Mast cells | Mast cells | Kit, Cpa3, Fcer1a |


---

## Software versions

| Package | Version |
|---------|---------|
| R | 4.4.0 |
| Seurat | 5.1.0 |
| celda (DecontX) | 1.20.0 |
| scDblFinder | 1.18.0 |
| CellChat | 1.6.1 |
| ComplexHeatmap | 2.20.0 |
| patchwork | 1.2.x |
| ggplot2 | 3.5.x |
| ggrepel | 0.9.x |

---

## Data availability

Raw FASTQ files are deposited at GEO under accession **[GSE######]** (to be added upon publication).

The processed Seurat object (`seurat_merged_clustered.rds`) is available at **[link]**.

Place downloaded data files in the `data/` directory before running scripts.

---

## Running the pipeline

Scripts are numbered in execution order. Edit the `input_dir` / `output_dir` paths at the top of each script to match your local environment, then run sequentially:

```bash
Rscript 01_QC_doublet_removal.R
Rscript 02_normalization_integration_clustering.R
Rscript 03_differential_expression.R
Rscript 04_cellchat_analysis.R
Rscript 05_cellchat_circle_plots.R
Rscript 06_cellchat_heatmap_pathways.R
Rscript 07_immune_subclustering.R
```

---

## Utility functions

Three helper scripts in `utilities/` are sourced by the main analysis scripts. They should not be run directly.

### `utilities/rankNet_patched.R`
A patched drop-in replacement for `CellChat::rankNet()` targeting `mode = "comparison"`.

**Key differences from the original:**
- `return.data = TRUE` returns `list(signaling.contribution = df, gg.obj = gg)`, making the underlying data frame directly accessible for export or custom plotting
- ON/OFF bars are colored gold (`#F1C40F`) and grey (`#A7A9AC`) by default to match the manuscript palette
- Non-significant pathways are removed cleanly rather than rendered in black

**Used in:** `04_cellchat_analysis.R` (Fig 5B), `06_cellchat_heatmap_pathways.R` (Fig 5C data extraction)

---

### `utilities/netAnalysis_signalingRole_heatmap2.R`
An extended version of `CellChat::netAnalysis_signalingRole_heatmap()`.

**Key differences from the original:**
- `color.use` accepts a **named vector** keyed by cell-type name (not just positional), so colors stay consistent when the heatmap column order changes between ON and OFF
- `top.bar.ylim` and `right.bar.ylim` parameters allow fixed axis scales across the ON and OFF heatmap panels, making cross-condition comparison valid

**Used in:** `04_cellchat_analysis.R` (Fig 5D, Fig 6G), `06_cellchat_heatmap_pathways.R`

---

### `utilities/violin_plot_celltype.R`
A self-contained violin-plot helper built on top of Seurat's `VlnPlot()` and `patchwork`.

**Exported functions:**

| Function | Purpose |
|----------|---------|
| `plot_violin_panels(seu, genes, de_table, ...)` | One-call combined figure: runs VlnPlot per gene, annotates each panel with `adj.P` and `log2FC` from a FindMarkers table, and assembles into a patchwork grid |
| `make_violin_plot_list(seu, genes, de_table, ...)` | Same logic but returns a plain list of ggplot objects for custom layouts |
| `make_de_label(gene, de_table, ...)` | Internal helper that formats the annotation string |
| `fmt_adj_p(p, ...)` | Internal helper that formats the adjusted p-value label |

**Key parameters for `plot_violin_panels()` / `make_violin_plot_list()`:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `seu` | — | Seurat object (already subsetted to the cell type of interest) |
| `genes` | — | Character vector of gene names to plot |
| `de_table` | — | `FindMarkers()` result with genes as rownames |
| `group.by` | `"treatment"` | Metadata column for x-axis grouping |
| `cols` | `c(ON="#F1C40F", OFF="#A7A9AC")` | Fill colors per group |
| `y_lim` | `c(0, 7)` | y-axis limits (passed to `coord_cartesian`) |
| `ncol` | `4` | Number of columns in the patchwork grid |
| `include_logfc` | `TRUE` | Whether to include log2FC in the annotation |

**Used in:** `03_differential_expression.R` (Sup Fig 7D, Sup Fig 7E, Sup Fig 7G)

---

## Contact

Questions: [yuedz@umich.edu](mailto:yuedz@umich.edu).
