The dataset was downloaded from Qiita (https://qiita.ucsd.edu/) with study ID 1011, which consists of 98 samples and 6674 OTUs. We selected the
caffeine intake as the outcome of interest since caffeine intake was found to have a significant impact on the
gut microbiota. We aimed to predict the caffeine intake based on the OTU abundances.
Before applying the prediction methods, we implemented a series of preprocessing steps designed in  to make the microbiome data more amenable to predictive modeling. First, we removed outlier
samples based on an outlier index defined on Bray-Curtis distance and removed rare OTUs with prevalence
less than 10% to reduce the dimensionality of OTUs, leaving 98 samples and 499 OTUs. Second, we
normalized OTU raw read counts using GMPR (Chen et al., 2018) followed by a replacement of outlier
counts using winsorization at 97% quantile. Third, we transformed the normalized OTU abundance data
using square-root transformation to reduce the influence of highly abundant observation. Finally, we applied
quantile transformation to the caffeine intake to make it approximately normally distributed.
