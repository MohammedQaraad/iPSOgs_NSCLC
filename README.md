An Optimized Computational Framework for Non-Small Cell Lung Cancer Subtype Classification and Biomarker Discovery, Artificial Intelligence Review Journal, Springer, 2025

Mohammed Qaraad 1   , Luke H. Hoeppner 1,2  , Bushra Shakir 3 , and David Guinovart 1*    
1 The Hormel Institute, University of Minnesota, 801 16th Ave NE, Austin, 55912, MN, USA, qaraa001@umn.edu,   hoepp005@umn.edu, and guino001@umn.edu 
2 Masonic Cancer Center, University of Minnesota, Minneapolis, MN, United States, hoepp005@umn.edu. 
3 Department of Medical Instrumentation Engineering Techniques, Imam Ja’afar Al-Sadiq University, Iraq, bushra.shaker@ijsu.edu.iq
* Corresponding author: David Guinovart, guino001@umn.edu



Abstract
Non-small cell lung cancer (NSCLC), primarily consisting of lung squamous cell carcinoma (LUSC) and lung adenocarcinoma (LUAD), is a significant cause of cancer-related death globally. Accurate subtyping and identifying reliable biomarkers for NSCLC are needed to design personalized therapy; however, the molecular heterogeneity and high dimensionality of gene expression data pose significant challenges. This work proposes iPSOgs, an enhanced particle swarm optimization algorithm with a golden section search refinement strategy and an adaptive crossover-based solution generation method to improve search efficiency and convergence stability. The XGBoost classifier's hyperparameter optimization and gene selection are done simultaneously using iPSOgs to predict NSCLC subtypes. On transcriptomic datasets (TCGA-LUAD, TCGA-LUSC, and GSE81089), the developed framework performs exceptionally well, achieving an accuracy of 0.9580, a receiver operating characteristic area under the curve of 0.9879, an F1 score of 0.9560, a recall of 0.9456, and a precision of 0.9668. The robustness of iPSOgs in managing complex optimization landscapes is further confirmed by comparison with cutting-edge metaheuristics on the CEC2017 test suite. The role of the biologically significant genes DSG3, KRT5, and SPRR2E, which were confirmed as involved in NSCLC pathogenesis by enrichment and protein-protein analysis, is further highlighted by the SHAP-based explanation. This work demonstrates the efficacy of iPSOgs as a diagnostic and biomarker discovery tool, offering a scalable and explainable solution for precision oncology in lung cancer.
