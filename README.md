# Transfer learning under large-scale low-rank regression models

## Directory
- [Function](https://github.com/ishspsy/transfer_learning/tree/main/Function): includes functions that implement transfer learning algorithms and generate synthetic data used in Section 3.
- [Data](https://github.com/ishspsy/transfer_learning/tree/main/Data): includes CCLE dataset.
- [Codeexample](https://github.com/ishspsy/transfer_learning/tree/main/Codeexample): includes R codes demonstrating the use of R functions in [Function](https://github.com/ishspsy/transfer_learning/tree/main/Function) directory.


## Data
The Cancer Cell Line Encyclopedia (CCLE) dataset consists of the expression levels of 19,221 genes and drug response data for 24 compounds across a range of human cancer cell lines, where the drug response of each cell line is measured as the area under the dose-response curve. Cancer cell line datasets, such as CCLE, have been widely used to construct predictive models of drug response. We aim to predict drug response in non-small cell lung cancer (NSCLC) cell lines with mutations in the Kirsten rat sarcoma viral oncogene homolog (KRAS) gene. We focus on a subset of the CCLE data. Specifically, we consider 8 drugs for which drug response is available for all cell lines. These 8 drugs are ``AZD0530", ``Crizotinib", ``Dovitinib", ``Lapatinib", ``Nutlin-3", ``PD0325901", ``TAE684", and ``Topotecan". In addition, we select the top 100 genes with the highest variance across KRAS-mutant NSCLC cell lines. Finally, we have 28 KRAS-mutant NSCLC cell lines on 100 predictors and 8 responses.


### Contact
ishspsy@yonsei.ac.kr
