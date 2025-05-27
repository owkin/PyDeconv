# Implemented methods

Name |  Link paper | Year | Code | Implemented | Verified*
--- | --- | --- | --- | --- | ---
OLS | [Wikipedia](https://en.wikipedia.org/wiki/Ordinary_least_squares) | | [R](https://github.com/dtsoucas/DWLS) | `OLS` | ✅
RLR | [Wikipedia](https://en.wikipedia.org/wiki/Robust_regression) | | [R](https://github.com/cran/MASS/blob/master/R/rlm.R) | `RLR` |
NNLS | [Deconvolution of blood microarray data identifies cellular activation patterns in systemic lupus erythematosus](https://pubmed.ncbi.nlm.nih.gov/19568420/) | 2009 | [R](https://github.com/dtsoucas/DWLS) | `NNLS` | ✅
DSA | [Digital sorting of complex tissues for cell type-specific gene expression profiles](https://pubmed.ncbi.nlm.nih.gov/23497278/) | 2013 | [R](https://github.com/zhandong/DSA) | |
CIBERSORT | [Robust enumeration of cell subsets from tissue expression profiles](https://pubmed.ncbi.nlm.nih.gov/25822800/) | 2015 | | `NuSVR` |
Bseq-SC | [A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure](https://pmc.ncbi.nlm.nih.gov/articles/PMC5228327/) | 2016 | [R](http://github.com/shenorrlab/bseq-sc) | |
EPIC | [Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data](https://elifesciences.org/articles/26476) | 2017 | [R](https://github.com/GfellerLab/EPIC) | |
DWLS | [Accurate estimation of cell-type composition from gene expression data](https://www.nature.com/articles/s41467-019-10802-z) | 2019 | [R](https://github.com/dtsoucas/DWLS) | `DWLS` | ✅
MuSiC | [Bulk tissue cell type deconvolution with multi-subject single-cell expression reference](https://pubmed.ncbi.nlm.nih.gov/30670690/) | 2019 | [R](https://github.com/xuranw/MuSiC) | `WNNLS` |
FARDEEP | [Fast and robust deconvolution of tumor infiltrating lymphocyte from expression profiles using least trimmed squares](https://pubmed.ncbi.nlm.nih.gov/31059559/) | 2019 | [R](https://github.com/cran/FARDEEP) |  |
CIBERSORTx | [Determining cell-type abundance and expression from bulk tissues with digital cytometry](https://pmc.ncbi.nlm.nih.gov/articles/PMC6610714/) | 2019 | | |
Bisque | [Accurate estimation of cell composition in bulk expression through robust integration of single-cell information](https://pmc.ncbi.nlm.nih.gov/articles/PMC7181686/) | 2020 | [R](https://github.com/cozygene/bisque) | |
SCDC | [SCDC: bulk gene expression deconvolution by multiple single-cell RNA sequencing references](https://academic.oup.com/bib/article/22/1/416/5699815) | 2020 | [R](https://github.com/meichendong/SCDC) |  |
Scaden | [Deep learning–based cell composition analysis from tissue expression profiles](https://www.science.org/doi/10.1126/sciadv.aba2619) | 2020 | [Python](https://github.com/KevinMenden/scaden) | `Scaden` | ✅
TAPE | [Deep autoencoder for interpretable tissue-adaptive deconvolution and cell-type-specific gene analysis](https://www.nature.com/articles/s41467-022-34550-9) | 2022 | [Python](https://github.com/poseidonchan/TAPE) | `Tape` | ✅
Kassandra | [Precise reconstruction of the TME using bulk RNA-seq and a machine learning algorithm trained on artificial transcriptomes](https://www.sciencedirect.com/science/article/pii/S1535610822003191) | 2022 | | |
ENIGMA | [Approximate estimation of cell-type resolution transcriptome in bulk tissue through matrix completion](https://pubmed.ncbi.nlm.nih.gov/37529921/) | 2023 | [R](https://github.com/WWXkenmo/ENIGMA?tab=readme-ov-file) | |

(*) verified: a test exists that compares the current package implementation with the original one.
