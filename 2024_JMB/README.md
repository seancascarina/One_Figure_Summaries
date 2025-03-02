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
- **Panel A:** Progressively removing I, L, and V amino acids from the green "sPrLD2" protein prevents it from going to SGs during stress. Notice that the bright green spots progressively disappear as more I, L, and V are removed from the PrLD.
- **Panel B:** I developed a Python-based webserver that automatically quantifies the ability of sPrLD to go to SGs, using pixel intensity data extracted directly from the yeast microscopy images. A high SG enrichment score indicates that the bright green spots in Panel A overlap with the bright red spots that show were SGs are located in the cell.
- **Panel C:** A statistical test ("Tukey's HSD" post-hoc test) indicates which groups significantly differ from the others, while controlling the family-wise error rate across multiple hypothesis tests.
- **Python code for analyses and data visualizations in Panels B and C are included in this repository.** Note that the plotting script creates a file that is required to run the ANOVA script, so they must be run sequentially.
