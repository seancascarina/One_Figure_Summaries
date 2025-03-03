# Discovering new SR/SR-related proteins using a serine and arginine composition grid search
## Background
- Many proteins involved in RNA processing contain regions with a high fraction of **arginine (R)** and **serine (S)* amino acids. 
- These types of protein regions are called **"RS domains"**, and they often help with RNA processing.
- Previous methods for finding new proteins with RS domains were limited and incomplete.

## Key Discoveries
- I developed a simple computational approach to find new proteins with RS domains.
- RS domains could be found with high sensitivity and specificity using a minimum threshold for combined S and R percentage.
- Analysis of existing data indicates that the newly identify proteins with RS domains have typical functions of RS proteins.

## Analysis Techniques
- Recently developed a Python program called "LCD-Composer", which could easily look for RS domains in protein sequences at a large scale.
- Ran a naive grid search over various R and S minimum thresholds, then examined what fraction of resulting proteins are known to bind to RNA.
- Validated newly identified RS proteins using Gene Ontology (GO) statistical analysis, analysis of existing data, and prediction models for RNA binding.

![RS_GridSearchFigure](https://github.com/seancascarina/One_Figure_Summaries/blob/main/2022_RNA/RS_GridSearchFigure.png)
**Figure credit:** Cascarina and Ross (2022) *RNA*.

## Brief Figure Summary
- **Panel A:** Progressively increasing the minimum composition threshold for R and S results in a decreasing number of proteins with RS domains.
- **Panel B:** Progressively increasing the minimum composition threshold for R and S leads to a higher fraction of RNA-binding proteins among the proteins with RS domains.
- Panels A and B together show that LCD-Composer searches using this naive grid search approach for R and S composition thresholds is extremely effective at discovering new proteins with RS domains with high sensitivity and specificity. 
- **Python code for analyses and data visualizations are included in this repository.**