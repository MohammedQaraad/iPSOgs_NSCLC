An Optimized Computational Framework for Non-Small Cell Lung Cancer Subtype Classification and Biomarker Discovery, Artificial Intelligence Review Journal, Springer, 2025

Mohammed Qaraad 1   , Luke H. Hoeppner 1,2  , Bushra Shakir 3 , and David Guinovart 1*    
1 The Hormel Institute, University of Minnesota, 801 16th Ave NE, Austin, 55912, MN, USA, qaraa001@umn.edu,   hoepp005@umn.edu, and guino001@umn.edu 
2 Masonic Cancer Center, University of Minnesota, Minneapolis, MN, United States, hoepp005@umn.edu. 
3 Department of Medical Instrumentation Engineering Techniques, Imam Ja’afar Al-Sadiq University, Iraq, bushra.shaker@ijsu.edu.iq
* Corresponding author: David Guinovart, guino001@umn.edu



Abstract
Non-small cell lung cancer (NSCLC), primarily consisting of lung squamous cell carcinoma (LUSC) and lung adenocarcinoma (LUAD), is a significant cause of cancer-related death globally. Accurate subtyping and identifying reliable biomarkers for NSCLC are needed to design personalized therapy; however, the molecular heterogeneity and high dimensionality of gene expression data pose significant challenges. This work proposes iPSOgs, an enhanced particle swarm optimization algorithm with a golden section search refinement strategy and an adaptive crossover-based solution generation method to improve search efficiency and convergence stability. The XGBoost classifier's hyperparameter optimization and gene selection are done simultaneously using iPSOgs to predict NSCLC subtypes. On transcriptomic datasets (TCGA-LUAD, TCGA-LUSC, and GSE81089), the developed framework performs exceptionally well, achieving an accuracy of 0.9580, a receiver operating characteristic area under the curve of 0.9879, an F1 score of 0.9560, a recall of 0.9456, and a precision of 0.9668. The robustness of iPSOgs in managing complex optimization landscapes is further confirmed by comparison with cutting-edge metaheuristics on the CEC2017 test suite. The role of the biologically significant genes DSG3, KRT5, and SPRR2E, which were confirmed as involved in NSCLC pathogenesis by enrichment and protein-protein analysis, is further highlighted by the SHAP-based explanation. This work demonstrates the efficacy of iPSOgs as a diagnostic and biomarker discovery tool, offering a scalable and explainable solution for precision oncology in lung cancer.

It provides complete, step-by-step scripts for reproducing all results in the paper, including data preprocessing, DESeq2 differential analysis, iPSOgs optimization, model training, external validation, and final interpretation.
All code in this repository follows the same workflow used in the manuscript.
## 📁 Repository Structure

```text
iPSOgs-NSCLC/
│
├──Step 1-(a) process_and_mapped RAW LUSC.R
├──Step 1-(b) process_and_mapped RAW LUAD.R
│      → Raw TCGA STAR-count preprocessing, gene mapping, QC, clinical matching
│
├──Step 2-merged.R
│      → Merges LUAD + LUSC datasets, aligns shared genes, merges labels
│
├──Step 3-DSeq2.R
│      → DESeq2 normalization, variance stabilization, QC, DGE results
│
├──Step 4-iPSOgs_NSCLC.py
│      → Nested CV, iPSOgs hyperparameter optimization, gene selection,
│        XGBoost model training, performance output, saved model artifacts
│
├──Step 5-Independant Dataset.py
│      → External validation (GSE81089) using the optimized model
│
├──optimized_xgboost_model.pkl
├──label_encoder.pkl
│
└──README.md
```


📦 Data Sources
TCGA-LUAD & TCGA-LUSC
As described in the manuscript, raw STAR-aligned read counts were downloaded from:
cBioPortal for Cancer Genomics (https://www.cbioportal.org/). 
The RNA-seq TCGA dataset, downloaded from (https://gdc.xenahubs.net) as TCGA-LUAD and TCGA-LUSC gene expression RNAseq - STAR - Counts
External Validation Dataset: GSE81089 Used for independent evaluation, obtained from NCBI GEO.
⚙️ Installation
R environment
Required for Steps 1–3.
install.packages(c("tidyverse", "reshape2", "ggplot2"))
BiocManager::install(c("DESeq2", "biomaRt"))

Python environment
Required for Steps 4–5.
conda create -n ipsogs python=3.10 -y
conda activate ipsogs
pip install numpy pandas scikit-learn xgboost shap matplotlib


▶️ Full Pipeline Execution
1. Process raw LUSC
   Rscript "Step 1-(a) process_and_mapped RAW LUSC.R"
2. Process raw LUAD
   Rscript "Step 1-(b) process_and_mapped RAW LUAD.R"
Outputs include:
TCGA-LUSC.star_counts_updated.csv
TCGA-LUAD.star_counts_updated.csv
TCGA-LUSC_matched_samples.csv
TCGA-LUAD_matched_samples.csv

3. Merge LUAD & LUSC
   Rscript Step 2-merged.R
Outputs:
TCGA-LUAD_LUSC_expression_combined.csv
TCGA-LUAD_LUSC_labels_combined.csv

4. Perform DESeq2 normalization
   Rscript Step 3-DSeq2.R
   
5. Train iPSOgs-optimized model
   python Step 4-iPSOgs_NSCLC.py

6. Independent external validation
   python Step 5-Independant Dataset.py

📊 Reproducing Figures & Tables
| Manuscript Item          | Script                        | Output               |
| ------------------------ | ----------------------------- | -------------------- |
| QC, normalization plots  | Step 3-DSeq2.R                | qc_plots/            |
| Convergence curves       | Step 4-iPSOgs_NSCLC.py        | figures/convergence/ |
| Optimizer comparison     | Step 4-iPSOgs_NSCLC.py        | figures/optimizer/   |
| Performance curves       | Step 4-iPSOgs_NSCLC.py        | figures/performance/ |
| External dataset results | Step 5-Independant Dataset.py | figures/external/    |
| Final biomarker panel    | Step 4-iPSOgs_NSCLC.py        | gene_panel.csv       |

⏱ Representative Runtimes
| Stage                           | Runtime       |
| ------------------------------- | ------------- |
| Raw preprocessing (LUAD/LUSC)   | 3–5 min each  |
| Merge + QC                      | < 1 min       |
| DESeq2 normalization            | < 1 min        |
| Nested CV + iPSOgs optimization | **10–15 min** |
| External validation             | < 1 min      |

🧠 Notes on Reproducibility
- All random seeds & CV folds are deterministic.

- DESeq2 processing is fully scripted.

- iPSOgs logs contain convergence curves & best-fitness records per fold.

- The final trained model and encoder are included in this repo.
