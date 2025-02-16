
from scipy.stats import f_oneway, tukey_hsd
import math
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
constructs = ['WT', '-3ILV', '-6ILV', '-9ILV', '-12ILV', '-15ILV', '-19ILV']

def main():

    prld = 'sPrLD2'
    output = open('sPrLD2_HydrophobicDeletionSeries_TukeyHSD_Posthoc_TestResults.tsv', 'w')
    output.write('\t'.join(['PrLD', 'Construct'] + constructs) + '\n')
    file = 'sPrLD2_AliphaticDeletionSeries_SummaryStatistics.tsv'
    medians_df = get_data(file)
    
    #CONSTRUCT matrix FROM DATA==================#
    matrix = []
    for construct in constructs:
        matrix.append(medians_df[construct])
    #============================================#

    # RUN STATISTICS USING matrix================#
    fstat, pval = f_oneway(*matrix)
    tukey_result = tukey_hsd(*matrix)
    #============================================#

    #OUTPUT P-VALUES FROM TUKEY'S HSD============#
    for index, row in enumerate(tukey_result.pvalue):
        construct = constructs[index]
        output.write('\t'.join([prld, construct] + [str(x) for x in row]) + '\n')
    output.close()
    #============================================#
    
    # CONVERT P-VALUES AND PLOT HEATMAP==========#
    log_pvals_matrix = [[math.log10(x) if x == 1.0 else -math.log10(x) for x in row] for row in tukey_result.pvalue]
    plot_pvalues_heatmap(log_pvals_matrix, constructs, prld)
    #============================================#
    
    
def plot_pvalues_heatmap(log_pvals_matrix, constructs, prld):
    """Plot heatmap of -log(p-value) for each pairwise comparison from Tukey's HSD test.
    Returns:
        None
    """
    
    mask = np.triu(log_pvals_matrix)
    ax = sns.heatmap(log_pvals_matrix, annot=True, mask=mask, fmt='.3f', annot_kws={'family':'Arial'})
    ax.set_facecolor('darkgrey')
    plt.xticks([x+0.5 for x in range(len(constructs))], labels=constructs)
    plt.yticks([x+0.5 for x in range(len(constructs))], labels=constructs, rotation=0)
    plt.xticks(fontname='Arial')
    plt.yticks(fontname='Arial')
    plt.xlabel('Comparison Group 1', fontname='Arial')
    plt.ylabel('Comparison Group 2', fontname='Arial')
    plt.title(prld)
    plt.savefig(prld + '_HydrophobicDeletions_ANOVA_TukeyHSD_Pvalues_LogScale.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()
    
    
def get_data(file):
    """Gather summary statistics data from file output by "plot_sPrLD_SGenrichmentScores.py"
    Returns:
        df = dictionary with construct as keys and list of median SG enrichment scores for each replicate as values.
    """
    
    h = open(file)
    header = h.readline()
    df = {}
    
    for line in h:
        construct, median, mean_of_medians, standard_error = line.rstrip().split('\t')
        construct, rep = construct.split('_')
        df[construct] = df.get(construct, [])
        df[construct].append( float(median) )
        
    h.close()
    
    return df


if __name__ == '__main__':
    main()