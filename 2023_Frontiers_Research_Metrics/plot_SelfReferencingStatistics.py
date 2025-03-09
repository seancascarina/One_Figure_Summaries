
import matplotlib.pyplot as plt
import seaborn as sns
import statistics
from scipy import stats
import pandas as pd
from matplotlib.lines import Line2D
import random

journal_sets = [['PlosBiol', 'PlosGenet', 'PlosPath', 'PlosMed', 'PlosCompBiol'], ['CurrBiol', 'Genetics', 'mBio', 'BMCmed', 'Bioinformatics', 'Cell', 'Nature', 'Science', 'MCB', 'JCellBiol']]
journal_label_sets = [['PLOS Biology', 'PLOS Genetics', 'PLOS Pathogens', 'PLOS Medicine', 'PLOS Computational\nBiology'], ['Current Biology', 'Genetics', 'mBio', 'BMC Medicine','Bioinformatics', 'Cell', 'Nature', 'Science', 'Molecular and\nCellular Biology', 'J Cell Biology']]
figure_panels = [('Fig1A', 'Fig1B'), ('Fig1C', 'Fig1D')]
fig_heights = [5, 6]
  

def main():

    output = open('TableS3_Summarized_SelfReferencingRate_DistributionStatistics_by_Journal.tsv', 'w')
    output.write('\t'.join(['Journal', 'Median Self-Reference Percentage', 'Average Self-Reference Percentage'] + ['Score at ' + str(x) + 'th Percentile' for x in (10, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90, 100)] + ['Median in '+str(pub_year) for pub_year in range(2003, 2023)]) + '\n')
        
    minimum_references = 20
    pmid_to_year = get_pub_year_df()
    
    for i, journals in enumerate(journal_sets):
    
        journal_labels = journal_label_sets[i]
        j_to_jlabel = {journal:journal_labels[i] for i, journal in enumerate(journals)}
        fig_height = fig_heights[i]
        fig_panel_tup = figure_panels[i]
        
        df = {'Publication Year':[],
            'Journal':[],
            'Self-Reference Percentage':[]}
                    
        df = get_selfref_rates(minimum_references, pmid_to_year, df, journals, j_to_jlabel)

        percentile_scores = []
        means = []
        boxplot_locs = []

        for xval, journal in enumerate(journal_labels):
            vals = [df['Self-Reference Percentage'][i] for i in range(len(df['Self-Reference Percentage'])) if df['Journal'][i] == journal]
            scores_at_percentiles = stats.scoreatpercentile(vals, [x for x in range(10, 110, 10)])
            boxplot_percentiles = stats.scoreatpercentile(vals, [25, 75])
            boxplot_locs.append(boxplot_percentiles)
            
            ave = statistics.mean(vals)
            means.append(ave)
            max_val = max(vals)
            
            pub_year_medians = []
            for pub_year in range(2003, 2023):
                temp_df = pd.DataFrame.from_dict(df)
                temp_df = temp_df[temp_df['Publication Year'] == pub_year]
                temp_df = temp_df[temp_df['Journal'] == journal]
                year_vals = temp_df['Self-Reference Percentage']
                if len(year_vals) > 0:
                    med = stats.scoreatpercentile(temp_df['Self-Reference Percentage'], 50)
                else:
                    med = 'N/A'
                pub_year_medians.append(med)
            
            output.write('\t'.join([journal.replace('\n', ' '), str(scores_at_percentiles[4]), str(ave)] + [str(x) for x in sorted(list(scores_at_percentiles) + list(boxplot_percentiles))] + [str(x) for x in pub_year_medians]) + '\n')

            percentile_scores.append(scores_at_percentiles[:-1])    # LEAVE OUT 100 PERCENTILE MARK FROM PLOTS, SINCE VIOLIN EXTENDS TO MAX VALUE.

        plot_maindistributions(df, journal_labels, percentile_scores, means, boxplot_locs, fig_panel_tup[0], fig_height)

        df = pd.DataFrame.from_dict(df)
        years = sorted(list(set(df['Publication Year'])))

        median_df = make_median_df(df, years, journal_labels)
        median_df = pd.DataFrame.from_dict(median_df)
        
        plot_year_lineplot(median_df, years, journals, journal_labels, fig_panel_tup[1], fig_height)
        
    # RUN CALCULATIONS FOR ALL PUBLICATIONS COMBINED
    scores_at_percentiles = stats.scoreatpercentile(df['Self-Reference Percentage'], [x for x in range(10, 110, 10)])
    boxplot_percentiles = stats.scoreatpercentile(df['Self-Reference Percentage'], [25, 75])
    pub_year_medians = []
    for pub_year in range(2003, 2023):
        temp_df = df[df['Publication Year'] == pub_year]
        year_vals = temp_df['Self-Reference Percentage']
        if len(year_vals) > 0:
            med = stats.scoreatpercentile(temp_df['Self-Reference Percentage'], 50)
        else:
            med = 'N/A'
        pub_year_medians.append(med)
    
    output.write('\t'.join(['All Journals Combined', str(scores_at_percentiles[4]), str(ave)] + [str(x) for x in sorted(list(scores_at_percentiles) + list(boxplot_percentiles))] + [str(x) for x in pub_year_medians]) + '\n')
    
    output.close()

    
def plot_maindistributions(plotting_df, journal_labels, percentile_scores, means, boxplot_locs, figure_panel, fig_height):

    plotting_df['hue'] = [0 for x in range(len(plotting_df['Journal']))]
    
    for xval, journal in enumerate(journal_labels):
        plotting_df['Journal'].append(journal)
        plotting_df['hue'].append(1)
        plotting_df['Self-Reference Percentage'].append(-1000)
        plotting_df['Publication Year'].append(1950)
    
    max_val = max(plotting_df['Self-Reference Percentage'])
    ax = plt.gca()

    colors = sns.color_palette('pastel')
    
    sns.violinplot(x='Self-Reference Percentage', y='Journal', hue='hue', data=plotting_df, inner=None, palette=[colors[2]], width=1.0, split=True, cut=0, orient='h', order=journal_labels)
    for xval, journal in enumerate(journal_labels):
        for j, score in enumerate(percentile_scores[xval]):
            if j == 4:
                color = '#d62728'
                offset = 0.15
                linewidth=2
            else:
                color = 'black'
                offset=0.10
                linewidth=1
            plt.plot((score, score), (xval-offset, xval+offset), color=color, linewidth=linewidth)
        plt.plot((means[xval], means[xval]), (xval-0.15, xval+0.15), color='#1f77b4', linewidth=2)

        boxplot_loc = boxplot_locs[xval]
        coords = [(boxplot_loc[0], xval-0.05), (boxplot_loc[1], xval-0.05), (boxplot_loc[1], xval+0.05), (boxplot_loc[0], xval+0.05)]
        ax = plt.gca()
        ax.add_patch(plt.Polygon(coords, color ='grey'))
    
    ax.get_legend().remove()
    
    plt.xticks(fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('Percent Self-References', fontname='Arial', fontsize=12)
    plt.ylabel('Journal', fontname='Arial', fontsize=12)
    plt.xlim(0-5, max_val+5)
    plt.ylim(len(journal_labels)-0.5, 0-1)
    fig = plt.gcf()
    fig.set_size_inches(4, fig_height)
    plt.savefig(figure_panel + '_SelfReferencing_Distributions_AllJournals.tif', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def plot_year_lineplot(df, years, journals, journal_labels, fig_panel, fig_height):

    # FILTER OUT DUMMY VALUES THAT GET ADDED IN plot_maindistributions() (THESE GET ADDED GLOBALLY EVEN THOUGH THE ADDING IS PERFORMED WITHIN A FUNCTION)
    df = df[df['Median Self-Reference Percentage'] > -1000]
    
    colors = sns.color_palette()
    sns.lineplot(x='Publication Year', y='Median Self-Reference Percentage', data=df, hue='Journal', palette=colors[:len(journal_labels)], hue_order=journal_labels)
    sns.scatterplot(x='Publication Year', y='Median Self-Reference Percentage', data=df, hue='Journal', palette=colors[:len(journal_labels)], hue_order=journal_labels)
    plt.xticks([x for x in years], labels=list(years), fontname='Arial', fontsize=12, rotation=45)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('Publication Year', fontname='Arial', fontsize=14)
    plt.ylabel('Median Percent Self-References', fontname='Arial', fontsize=14)
    
    ax = plt.gca()

    legend_elements = [Line2D([0], [0], marker='o', color=colors[i], label=journal_labels[i], markerfacecolor=colors[i], markersize=6) for i in range(len(journals))]
    ax.legend(handles=legend_elements, loc='upper right', prop={'family':'Arial', 'size':9})
    plt.xlim(2002.5, 2022.5)
    plt.ylim(4, 21)
    fig = plt.gcf()
    fig.set_size_inches(8, fig_height)
    plt.savefig(fig_panel + '_MedianSelfReference_Percentage_ByYear_Linegraph_AllJournals.tif', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def make_median_df(df, years, journal_labels):

    median_df = {'Publication Year':[],
        'Journal':[],
        'Median Self-Reference Percentage':[]}
    for year in years:
        for journal in journal_labels:
            df2 = df[df['Publication Year'] == year]
            df3 = df2[df2['Journal'] == journal]
            vals = list(df3['Self-Reference Percentage'])
            if len(vals) == 0:
                continue
            med_selfref = stats.scoreatpercentile(vals, 50)
            median_df['Publication Year'].append(year)
            median_df['Median Self-Reference Percentage'].append(med_selfref)
            median_df['Journal'].append(journal)
            
    return median_df
    
    
def get_selfref_rates(minimum_references, pmid_to_year, df, journals, j_to_jlabel):

    h = open('TableS1_SelfReferencingRate_Estimates.tsv')
    header = h.readline()

    for line in h:
        journal, pmid, authors, all_author_num_selfrefs, all_author_perc_selfrefs, max_author, num_selfrefs, total_refs, perc_selfrefs = line.rstrip().split('\t')
        
        if pmid == '2269344':   # ONE PMID WAS MIS-ATTRIBUTED TO Nature, WHEN IT IS REALLY A FEBS Letters PUBLICATION
            continue
        
        if journal not in journals:
            continue
            
        journal_label = j_to_jlabel[journal]
        num_selfrefs, total_refs = int(num_selfrefs), int(total_refs)
        if total_refs < minimum_references:
            continue

        perc = num_selfrefs / total_refs * 100
        df['Publication Year'].append( int(pmid_to_year[pmid]) )
        df['Self-Reference Percentage'].append( perc )
        df['Journal'].append(journal_label)
        
    h.close()
    
    return df
    

def get_pub_year_df():

    h = open('Filtered_PMIDs_AllJournals.csv')
    header = h.readline()
    
    df = {}
    for line in h:
        items = line.rstrip().split(',')
        year = items[3]
        pmid = items[9]
        if pmid == '':
            continue
        df[pmid] = year
        
    h.close()
    
    return df


if __name__ == '__main__':
    main()
    