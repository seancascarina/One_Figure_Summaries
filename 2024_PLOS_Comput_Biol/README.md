# How common are low-complexity domains in different types of organisms?
## Background
- **"Low-complexity domains" (LCDs)** are regions of protein sequences made up of only a few amino acids (e.g., "NNQNNQYQNQNYYQNQQQNNQNN" is an LCD).
- There are at least 400 types of LCDs, depending on which amino acids are enriched in the LCD sequence.
- LCDs have unusual biochemical properties, which affects what they do in living organisms and how common they are in each organism. *However, they prevalence and functions of LCDs across organisms had not been studied on a large scale.*

## Key Discoveries
- LCDs differ in their abundance across organisms, potentially explaining how organisms adapt to different environments.
- Statistically significant LCD enrichment differs between organisms and depends on LCD class. Statistically significant LCD enrichment can occur even when the amino acids are extremely common in the proteome (e.g., N, D, and E LCDs in the malaria proteome).
- Many types of LCDs are often associated with unique functions that depend on the amino acids enriched in the LCD.

## Analysis Techniques
- I recently developed a Python program called "LCD-Composer" to look for LCDs in complete proteomes specifically by LCD type.
- These searches required massive computational parallelization: >21k organisms (representing all known organisms), each with between 5k and 100k proteins, each analyzed 800 times with LCD-composer.
- Generated and analyzed randomized proteomes to evaluate the frequencies of LCDs expected by random change ("null" expectation), then compared the observed LCD frequencies to the null expectation to determine statistically significant LCD enrichment using Fisher's exact test. This included a multiple-test correction (Sidak) to control the Type I error rate.

![LCD figure](https://github.com/seancascarina/One_Figure_Summaries/blob/main/2024_PLOS_Comput_Biol/LCD_Frequencies_Malaria_and_Humans.png)
**Figure credit:** Cascarina and Ross (2024) *PLOS Comput Biol*.

## Brief Summary
- **Panel A:** Certain types of LCDs (bright colors) are common in the malaria proteome. This heatmap shows the natural log of the odds ratio for each type of LCD compared to a randomized version of the Malaria proteome. The "Primary Amino Acid" (x-axis) essentially represents the most prevalent amino acid in the LCD, whereas the "Secondary Amino Acid" (y-axis) represents the next most prevalent amino acid in the LCD. Together, these two amino acids comprise the "LCD class".
- **Panel B:** Indication of statistical significance for LCD enrichment according to each LCD class. Red squares indicate that the observed LCD enrichment is statistically significant.
- **Panel C:** Same as Panel A, but for the human proteome. Notice that different classes of LCDs are enriched in the human proteome compared to the malaria proteome.
- **Panel D:** Same as Panel B, but for the human proteome. Notice that more classes and different classes of LCDs are statistically significant compared to the malaria proteome.
- **Python code for analyses and data visualizations for all Panels are included in this repository.** Note that the "compare_Observed_vs_Scrambled_Frequencies.py" script creates a file that is required to run the "plot_IndividualOrganism_lnORs_and_Pvals.py" script, so they must be run sequentially.
