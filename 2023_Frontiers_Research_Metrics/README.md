# Quantifying self-referencing rates in journal publications for biological disciplines
## Background
- Authors often cite their own previous publications (**"self-referencing"**) to communicate work that served as a foundation for a new study.
- No clear guidelines exist regarding an appropriate or typical self-referencing rate.
- The number of citations to an author's publications is often one of the metrics used during hiring decisions, grant funding decisions, journal peer-review, and performance evaluations. This can lead to an incentive for researchers to cite their own work more frequently.

## Key Discoveries
- Median self-referencing rates in biological disciplines are between ~8-13% across many journals.
- Right-skewed distributions indicate that some publications have relatively high self-referencing rates, and this occurs for all journals evaluated.
- Self-referencing rates depend on other factors such as the total number of references in the publication, the number of authors on the publication, the author status/rank, author position in the author list, and total number of publications for each author. These effects were all quantified in more detail in the study.

## Analysis Techniques
- Used the Entrez API to retrieve publication metadata for >94,000 publications from PubMed.
- Coerced seaborn's violinplot() function to make what appear to be a series of histograms with percentile markers to show self-referencing rate distributions for each journal. In reality, these are rotated violin plots with only half of each violin shown.
- Performed ordinary least-squares regression to confirm accuracy of the automated quantification method for self-referencing rates, and *t*-tests to statistically compare groups of data.

![self-referencing figure](https://github.com/seancascarina/One_Figure_Summaries/blob/main/2023_Frontiers_Research_Metrics/SelfReferencing_Figure.png)
**Figure credit:** From Cascarina (2023) *Frontiers in Research Metrics and Analytics*.

## Brief Figure Summary
- **Panel A:** Median self-referencing rates and the shapes of self-referencing rate distributions vary across the PLOS journals. PLOS journals were initially chosen because they were open-access from their inception, so their entire publication history could be analyzed.
- **Panel B:** Self-referencing rates are relatively stable over time for the PLOS journals.
- **Panel C:** Self-referencing rates vary more substantially across a broader range of journals (non-PLOS journals).
- **Panel D:** Median self-referencing rates for non-PLOS have tended to converge over time toward an industry-standard ~8-13%.
- **Python code for analyses and data visualizations are included in this repository.** Note that the QueryPubmed.py is included for demonstration purposes and will not run as-is because it requires a user email and a user-specific API key for making queries to Pubmed.