This repository was created to regenerate the Figures and reproduce the results from Akana, Yoe et al., under review.

**Instructions**

To repeat the analyses and regenerate the figures:
1. Download and unzip the zip files from the Figshare repository. Downloading the code should take a few minutes.
2. Change the setup.R code to include the path to the data, input, and output directories
3. Source setup.R with the required R pacakges
4. Run main.R

The figures will be saved in the output directory.
The expected output for figures will also be provided here upon manuscript publication.

**Run time:** Expected run time if reanalyzing the data to reprocude the results is 3-4 hours. Figures can be reproduced within minutes.

**Hardware requirements**
Running the code requires only a standard computer with enough RAM to support the in-memory operations.

**Software requirements**
**OS Requirements**
The code was developed and tested on:
macOS: Sequoia (15.6.1)

**Dependencies (see setup.R):** scater (v1.26.1); edgeR (v3.40.2); beanplot (v1.3.1); cowplot (v1.1.3); Seurat (v5.0.1); EBImage (v4.40.1); survival (v3.5-7); rms (v6.7-1); mixtools (v2.0.0); MASS (v7.3-60); ggplot2 (v3.5.1); nnet (v7.3-19); ppcor (v1.1); ROCR (v1.0-11); tsne (v0.1-3.1); gplots (v3.2.0); ggpubr (v0.6.0); EnhancedVolcano (v1.16.0); plyr (v1.8.9); reshape2 (v1.4.4); plotrix (v3.8-4); stats (v4.2.0); Matrix (v1.6-4); Rtsne (v0.17); lmerTest (v3.1-3); devtools (v2.4.5); gplots (v3.2.0); heatmap3 (v1.1.9); e1071 (v1.7-14); openxlsx (v4.2.5.2); RColorBrewer (v1.1-3)

**License** 

BSD 3-Clause License provided ([here](https://github.com/Jerby-Lab/Immunomodulators_PertrubSeq/LICENSE)).
