# Discovering protein features that affect their stress-response behavior in cells
## Background
- Certain proteins form clumps called **"stress granules" (SGs)** in response to hot temperatures.
- Our lab studies **"prion-like domains" (PrLDs)**, which are a certain type of protein that often goes to SGs.
- No one really knows what about the protein sequences causes only specific PrLDs to go to SGs.

## Key Discoveries
- Hydrophobic amino acids (abbreviated I, L, and V) cause PrLDs to go to SGs.
- Progressively removing I, L, or V from PrLDs causes a step-wise reduction in the ability to go to SGs.

## Analysis Techniques
- Developed local, Python-based webserver to automate the quantification of SG localization in yeast microscopy images.
- Integrated an open-source neural network for initial identification of yeast cells in microscopy images.
- Performed one-way ANOVA with Tukey's HSD post-hoc test (scipy.stats) to determine statistically significant differences in SG localization between proteins.
- Used scipy.optimize.curve_fit to fit sigmoidal objective function to median SG Enrichment Scores.

![SG figure](https://github.com/seancascarina/One_Figure_Summaries/blob/main/2024_JMB/SG_figure.png)
**Figure credit:** Modified from Baer, Cascarina, and Paul *et al.* (2024) *JMB*.

## Brief Figure Summary
- **Panel A:** Progressively removing I, L, and V amino acids from the green "sPrLD2" protein prevents it from going to SGs during stress.
- **Panel B:** I developed a Python-based webserver that automatically quantifies the ability of sPrLD to go to SGs, using pixel intensity data extracted directly from the yeast microscopy images.
- **Panel C:** A statistical test ("Tukey's HSD" post-hoc test) indicates which groups significantly differ from the others, while controlling the family-wise error rate across multiple hypothesis tests.
- **Python code for analyses and data visualizations in Panels B and C are included in this repository.** Note that the plotting script creates a file that is required to run the ANOVA script, so they must be run sequentially.

## Additional Details
**Detailed Background:**
Cells can form clumps of protein and RNA called "stress granules" (SGs) in response to hot temperatures. Only certain proteins go to SGs during stress, but it is not clear why some proteins go to SGs and others do not. All proteins are made of building blocks called "amino acids". There are 20 different types of amino acids, and their order and type in a protein sequence determines how each protein behaves in cells. Our lab discovered that specific amino acids - in this case, the amino acids abbreviated I, L, and V - control the ability of certain proteins called "prion-like domains" (PrLDs) to join SGs in cells. The figure shows microscopy images of yeast cells producing these PrLDs with progressively decreasing amounts of I, L, and V in their sequences, along with corresponding single-cell quantification of the amount of the PrLD in SGs.

**Figure Panel Descriptions:**
**Panel A** shows images of yeast cells with two different proteins shown in different colors. The green color shows us where our protein of interest is in cells. The red color shows us where SGs are located in cells. Notice that the red protein always forms clumps that appear as bright spots. Compare that with the green "WT" protein, which starts out forming bright green spots at the same locations as the red bright spots. However, as more I, L, and V are removed (-3 ILV, -6 ILV, ... , -19 ILV), the green spots gradually become dimmer and dimmer before disappearing entirely around -12 ILV, even though the bright red spots still occur in the same cells.

**Panel B** shows a quantification of the behavior of the green protein using data extracted from the yeast microscopy images. I developed a Python-based webserver that automatically quantifies single-cell "SG enrichment scores" for our yeast microscopy images. The SG enrichment score on the y-axis measures how strongly our protein of interest forms bright spots that coincide with SGs: the better the protein is at forming these bright spots, the higher the SG score will be. The "Construct" on the x-axis indicates how many I, L, and V have been removed from the sequence and corresponds to the images show in Panel A. In simple terms, the SG enrichment score is basically a ratio of the green pixel intensities ***within SG regions*** to the green pixel intensities ***outside of the SG regions***. Notice that as the green spots become dimmer and dimmer, the SG enrichment scores decrease until -12 ILV, where the scores level out at ~1.3 (the score at which no detectable bright spots occur). The panel includes a curve fit using scipy.optimize.curve_fit with a sigmoidal objective function.

**Panel C** shows the results of statistical comparisons of each construct to all other constructs using Tukey's "Honestly Significant Difference" (HSD) statistical test. The heatmap shows -log(p-value) for each comparison, which means that larger values indicate more statistically significant differences. The traditional p<0.05 threshold corresponds to a -log(p-value)>1.30, meaning that any value in the heatmap greater than 1.30 would be considered statistically significant using traditional p-value thresholds. Notice that constructs with very different SG enrichment scores in Panel B have very high values in this panel, indicating that the comparison is very statistically significant. For example, comparing WT and -19ILV has the highest value in the heatmap and these groups are the furthest apart in Panel B. Comparing -6ILV and -9ILV has a smaller value in the heatmap and are closer together in Panel B: there is still a statistically significant difference between these groups, but the difference is smaller than the WT versus -19ILV comparison. Finally a comparison such as WT versus -3ILV has a small difference in Panel B and a correspondingly small value in the heatmap, indicating that these differences are not statistically significantly different.
Note that a one-way ANOVA test determined that SG enrichment scores differed significantly between constructs prior to running Tukey's HSD post-hoc test, but the pairwise group comparisons from Tukey's HSD were more important for our purpose than the omnibus ANOVA test results.
