# Complex networks for analyzing the urban acoustic environment

Repsoitory to store the Preprint, analyses and code for Haselhoff et al. 2023.

## Files
:page_facing_up: **Complex networks for analyzing the urban acoustic environment.pdf** -> Preprint Manucript.

:arrow_forward: **AnalysesCode.py** -> Code used for analyses in Haselhoff et al. 2023. The data in this study is available on request from the corresponding author. The data is not publicly available due to potential privacy issues.

:arrow_forward: **Appendix C.py** -> Refurbished code and functions to apply Frequency Correlation Matrices based Complex Networks for data from the urban acoustic environment (with demo on test_data.csv).

:clipboard:	 **test_data.csv** -> Example Dataset of normalized and de-noised frequency spectra (s. Haselhoff et al. 2023: 2.1 & 3.1) for one day of recordings from one recording device. Use this for a straight forward application for functions in Appendix C.


## Abstract

The urban acoustic environment (AE) provides comprehensive acoustic information related to the diverse systems of urban areas, such as traffic, the built environment, or biodiversity. The decreasing cost of acoustic sensors and rapid growth of storage space and computational power have fostered the collection of large amounts of acoustical data to be processed. However, despite the extensive information that is recorded by modern acoustic sensors, few approaches are established to capture the rich complex dynamics embedded in the time-frequency domain of the urban AE. Quantitative methods need to account for this complexity, while effectively reducing the high dimensionality of acoustic features within the data. Therefore, we introduce complex networks as a tool for analyzing the complex structure of large-scale urban AE data. We present a framework to construct networks based on frequency correlation matrices (FCMs). FCMs have shown to be a promising tool to depict environment specific interrelationships between consecutive power spectra. Accordingly, we show the capabilities of complex networks for the quantification of these interrelationships and thus, to characterize different urban AEs. 
We demonstrate the scope of the proposed method, using one of the world’s most extensive longitudinal audio datasets, considering 3-min audio recordings (n = 319,385 ≙ 665 days) from 23 sites. We construct networks from hour-of-day specific audio recordings for each site. We show that the average shortest path length (ASPL) as an indicator for dominance of sound sources in the urban AE exhibits spatial- and temporal-specific patterns between the sites, which allows us to identify four to seven clusters of distinct urban AEs. To validate our findings, we use the land use mix around each site as a proxy for the AE and compare those between and within the clusters. The identified clusters show high intra- and low inter-cluster correlations of ASPL diel cycles as well as strong intra-similarities in land use mix. Our results indicate that complex networks are a promising approach to analyze large-scale audio data, expanding our understanding of the time-frequency domain of the urban AE.

### Literature
_Haselhoff, T., Braun, T., Fiebig, A., Hornberg, J., Lawrence, B. T., Marwan, N., & Moebus, S. (2023). Complex networks for analyzing the urban acoustic environment. Earth ArXiv. https://doi.org/10.31223/X5J08D._
