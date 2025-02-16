
import matplotlib.pyplot as plt
from matplotlib import patches
import seaborn as sns
from scipy.stats import sem
import statistics
import pandas as pd
from scipy.optimize import curve_fit
import numpy as np

def main():
    
    # GATHER RAW DATA==============================================#
    data_label = 'sPrLD2_AliphaticDeletionSeries'
    label_order = ['WT', '-3ILV', '-6ILV', '-9ILV', '-12ILV', '-15ILV', '-19ILV']
    file = 'sPrLD2_HydrophobicDeletionSeries_CombinedData.tsv'
    df, df2 = get_data(file)      
    #==============================================================#

    # CALCULATE MEANS AND STANDARD ERRORS==========================#
    sem_output = open(data_label + '_SummaryStatistics.tsv', 'w')
    sem_output.write('\t'.join(['Data Label', 'Median SG Enrichment Score', 'Mean of Medians', 'Standard Error of the Mean']) + '\n')
    sem_df, sem_output = calc_sems(df2, label_order, sem_output)
    sem_output.close()
    #==============================================================#
    
    #PLOTTING PARAMETERS==========================================#
    include_title = False
    apply_xtickrotation = 0
    width_scaling = 0.90
    ymin = 0
    ymax = None
    ybreakpoint = 6
    dotsize = 2.0
    dotopacity = 0.3
    grouped_dotsize = 2.5
    grouped_dotopacity = 0.5
    include_curvefit = True
    gfp_mean_of_medians = 1.313476665626228     # CALCULATED FROM THE GFP REPLICATES
    
    curvefit_estimates = []
    if include_curvefit:
        curvefit_estimates = get_curvefit_estimates(sem_df)
    #==============================================================#
    
    # CALL PLOTTING FUNCTIONS==========================================#
    plotting(df, df2, sem_df, include_curvefit, curvefit_estimates, data_label, include_title, label_order, apply_xtickrotation, width_scaling, ymin, ymax, ybreakpoint, dotsize, dotopacity, gfp_mean_of_medians)
    boxplot(df, df2, sem_df, data_label, include_title, label_order, apply_xtickrotation, width_scaling, ymin, ymax, ybreakpoint, grouped_dotsize, grouped_dotopacity, gfp_mean_of_medians)
    if ybreakpoint:
        boxplot_breakpoint(df, df2, sem_df, data_label, include_title, label_order, apply_xtickrotation, width_scaling, ymin, ymax, ybreakpoint, grouped_dotsize, grouped_dotopacity, gfp_mean_of_medians)
        plotting_breakpoint(df, df2, sem_df, include_curvefit, curvefit_estimates, data_label, include_title, label_order, apply_xtickrotation, width_scaling, ymin, ymax, ybreakpoint, dotsize, dotopacity, gfp_mean_of_medians)
    #==============================================================#
    
    
def get_curvefit_estimates(sem_df):
    """Calculate initial curvefit estimates.
    Returns:
        curvefit_estimates = a list of initial curvefit estimates to be used in scipy.optimize.curve_fit
    """
    
    max_yval = max(sem_df['Mean of Medians'])
    min_yval = min(sem_df['Mean of Medians'])
    midpoint_xval_est = len(set(sem_df['Construct'])) / 2
    L_est = sem_df['Mean of Medians'][-1] - sem_df['Mean of Medians'][0]
    if L_est < 0:
        b0_est = max(sem_df['Mean of Medians'])
    else:
        b0_est = min(sem_df['Mean of Medians'])
    curvefit_estimates = [L_est, midpoint_xval_est, 1, b0_est]
    
    return curvefit_estimates
    
    
def get_data(file):
    """Gather raw SG Enrichment Score data.
    Returns:
        df = long-form dictionary with raw, single-cell SG enrichment scores for each construct and replicate.
        df2 = long-form summary statistics dictionary with median SG enrichment scores for each construct and replicate.
    """
    
    df = {'Construct':[],
            'SG Enrichment Score':[], 
            'Replicate':[]
            }
    h = open(file)
    header = h.readline()
    
    for line in h:
        construct, capture, cellnum, sg_score = line.rstrip().split('\t')
        construct, replicate = construct.split('_')
        sg_score = float(sg_score)
        df['Construct'].append(construct)
        df['SG Enrichment Score'].append(sg_score)
        df['Replicate'].append(int(replicate))
        
    h.close()
    
    temp = pd.DataFrame.from_dict(df)
    df2 = temp.groupby(['Construct', 'Replicate']).median().reset_index()
    
    return df, df2
    
    
def sigmoidal_objective(x, L ,x0, k, b):
    """Sigmoidal objective function."""
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return y
    

def calc_sems(df2, label_order, sem_output):
    """Calculate means and standard errors for each construct.
    Returns:
        sem_df = long-form dictionary with means of medians and standard errors of the means for each construct.
    """
    
    sem_df = {'Construct':[],
            'Mean of Medians':[], 
            'Standard Error of Mean' :[]
            }

    for label in label_order:
        temp = df2[df2['Construct'] == label]
        medians = []
        for repnum in sorted(list(set(temp['Replicate']))):
            temp2 = temp[temp['Replicate'] == repnum]
            val = list(temp2['SG Enrichment Score'])[0]
            medians.append(val)

        ave = statistics.mean(medians)
        standard_error = sem(medians)
        for i in range(len(medians)):
            sem_output.write('\t'.join([label + '_' + str(i+1), str(medians[i]), str(ave), str(standard_error)]) + '\n')
        sem_df['Construct'].append(label)
        sem_df['Mean of Medians'].append(ave)
        sem_df['Standard Error of Mean'].append(standard_error)

    return sem_df, sem_output
            

def plotting(df, df2, sem_df, include_curvefit, curvefit_estimates, data_label, include_title, label_order, apply_xtickrotation, width_scaling, ymin, ymax, ybreakpoint, dotsize, dotopacity, gfp_mean_of_medians):
    """Make superplot of SG enrichment scores.
    Returns:
        None
    """
    
    num_constructs = len(set(df['Construct']))

    pal = sns.color_palette()
    max_repnum = max(df['Replicate'])
    pal = pal[:max_repnum]
    
    ax = sns.barplot(x='Construct', y='Mean of Medians', data=sem_df, color='white', linewidth=0, fill=False, yerr=sem_df['Standard Error of Mean'], zorder=101)
    
    #=========================
    # SETTING OF DIFFERENT EDGE COLORS FOR BARS IN BARPLOT WAS ADAPTED FROM:
    # https://stackoverflow.com/questions/65781401/is-it-possible-to-set-different-edgecolors-for-left-and-right-edge-of-matplotlib
    
    bars = [rect for rect in ax.get_children() if isinstance(rect, patches.Rectangle)]
    for bar in bars[:-1]:
        x, y = bar.get_xy()
        w, h = bar.get_width(), bar.get_height()

        ax.plot([x, x + w], [y + h, y + h], color='black', lw=1, zorder=100)

    #=========================
    
    sns.stripplot(x='Construct', y='SG Enrichment Score', data=df, hue='Replicate', palette=pal, alpha=dotopacity, size=dotsize, linewidth=0)
    
    sns.swarmplot(x='Construct', y='SG Enrichment Score', data=df2, hue='Replicate', palette=pal, linewidth=0.5)
    
    if include_curvefit:
        params_output = open(data_label + '_CurvefitParameters.tsv', 'w')
        params_output.write('Parameter Type\tEstimated Value\n')
        
        # RUN CURVE FITTING USING SIMULATED DATA
        popt, pcov = curve_fit(sigmoidal_objective, [x for x in range(len(sem_df['Mean of Medians']))], sem_df['Mean of Medians'], curvefit_estimates, method='lm')
        
        # UNPACK PARAMETER ESTIMATES
        est_L, est_x0, est_k, est_b = popt
        param_labels = ['L: difference between y-intercept and other y-limit', 'x0: x-axis value at curve midpoint', 'k: controls steepness of curve', 'b: y-intercept (one of the y-limits)']
        for i, param in enumerate([est_L, est_x0, est_k, est_b]):
            param_label = param_labels[i]
            params_output.write(param_label + '\t' + str(param) + '\n')
        params_output.close()
            
        sim_xvals = np.linspace(0, len(set(sem_df['Construct'])), num=1000, endpoint=True)
        sim_yvals = [sigmoidal_objective(x, est_L, est_x0, est_k, est_b) for x in sim_xvals]

        plt.plot(sim_xvals, sim_yvals, color='0.5')
        
    plt.xticks([x for x in range(len(label_order))], labels=label_order, fontname='Arial', fontsize=14, rotation=apply_xtickrotation)
    plt.plot((-0.5, len(label_order)), (gfp_mean_of_medians, gfp_mean_of_medians), linestyle='--', color='0.5', zorder=0)
    plt.yticks(fontname='Arial', fontsize=14)
    plt.ylabel('SG Enrichment Score', fontname='Arial', fontsize=16)
    plt.xlabel('Construct', fontname='Arial', fontsize=16)
    if include_title:
        plt.title(data_label, fontname='Arial', fontsize=18)

    ax.get_legend().remove()
    
    if ymax:
        plt.ylim(ymin, ymax)
    else:
        plt.gca().set_ylim(bottom=ymin)
        
    if ybreakpoint:
        plt.ylim(ymin, ybreakpoint)
        ax.spines.top.set_visible(False)
        ax.set_title('')
    else:
        plt.gca().set_ylim(bottom=ymin)
        
    fig = plt.gcf()
    fig.set_size_inches(width_scaling*num_constructs, 4)
    plt.savefig(data_label + '_PartitionCoefficients_Superplot.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()

    
def plotting_breakpoint(df, df2, sem_df, include_curvefit, curvefit_estimates, data_label, include_title, label_order, apply_xtickrotation, width_scaling, ymin, ymax, ybreakpoint, dotsize, dotopacity, gfp_mean_of_medians):
    """Plot the upper half of the plot if a ybreakpoint is used. Some plots have extreme outliers,
        so it may be necessary to condense the top half of the y-axis to visualize the 
        central tendency of the data.
    Returns:
        None
    """
    num_constructs = len(set(df['Construct']))

    pal = sns.color_palette()
    max_repnum = max(df['Replicate'])
    pal = pal[:max_repnum]
    
    df = pd.DataFrame.from_dict(df)
    high_partco_df = df[df['SG Enrichment Score'] >= ybreakpoint]

    ax = sns.stripplot(x='Construct', y='SG Enrichment Score', data=df, hue='Replicate', palette=pal, alpha=dotopacity, size=dotsize, linewidth=0)

    plt.xticks([x for x in range(len(label_order))], labels=label_order, fontname='Arial', fontsize=14, rotation=apply_xtickrotation)
    plt.xlabel('')
    plt.yticks(fontname='Arial', fontsize=14)
    plt.ylabel('SG Enrichment Score', fontname='Arial', fontsize=16, color='white')
    if include_title:
        plt.title(data_label, fontname='Arial', fontsize=18)

    ax.get_legend().remove()
    ax.spines.bottom.set_visible(False)
    ax.set_xticks([])

    plt.ylim(ybreakpoint, max(df['SG Enrichment Score']) + (0.07 * (max(df['SG Enrichment Score']) - ybreakpoint)))
    fig = plt.gcf()
    fig.set_size_inches(width_scaling*num_constructs, 1)
    plt.savefig(data_label + '_PartitionCoefficients_Superplot_TopHalfOnly_REVISED.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()
    
    
def boxplot(df, df2, sem_df, data_label, include_title, label_order, apply_xtickrotation, width_scaling, ymin, ymax, ybreakpoint, grouped_dotsize, grouped_dotopacity, gfp_mean_of_medians):
    """Make boxplots of SG enrichment scores. This differs from a superplot in that each replicate
        gets a separate boxplot to better visualize their distributions.
    Returns:
        None
    """
    
    num_constructs = len(set(df['Construct']))

    pal = sns.color_palette()
    max_repnum = max(df['Replicate'])
    pal = pal[:max_repnum]
    
    ax = sns.boxplot(x='Construct', y='SG Enrichment Score', data=df, hue='Replicate', palette=['0.8'], showfliers=False)
    sns.stripplot(x='Construct', y='SG Enrichment Score', data=df, hue='Replicate', palette=pal, alpha=grouped_dotopacity, size=grouped_dotsize, linewidth=0, dodge=True)
    
    plt.plot((-0.5, len(label_order)), (gfp_mean_of_medians, gfp_mean_of_medians), linestyle='--', color='0.5', zorder=0)

    plt.xticks([x for x in range(len(label_order))], labels=label_order, fontname='Arial', fontsize=14, rotation=apply_xtickrotation)

    ax.get_legend().remove()
    
    plt.yticks(fontname='Arial', fontsize=14)
    plt.ylabel('SG Enrichment Score', fontname='Arial', fontsize=16)
    if include_title:
        plt.title(data_label, fontname='Arial', fontsize=18)
    
    if ymax:
        plt.ylim(ymin, ymax)
    else:
        plt.gca().set_ylim(bottom=ymin)
        
    if ybreakpoint:
        plt.ylim(ymin, ybreakpoint)
        ax.spines.top.set_visible(False)
        ax.set_title('')
    else:
        plt.gca().set_ylim(bottom=ymin)
        
    fig = plt.gcf()
    fig.set_size_inches(width_scaling*num_constructs, 4)
    plt.savefig(data_label + '_PartitionCoefficients_GroupedBoxplots_SeparatedBoxes.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()
    
    
def boxplot_breakpoint(df, df2, sem_df, data_label, include_title, label_order, apply_xtickrotation, width_scaling, ymin, ymax, ybreakpoint, grouped_dotsize, grouped_dotopacity, gfp_mean_of_medians):
    """Plot the upper half of the boxplot if a ybreakpoint is used. Some plots have extreme outliers,
        so it may be necessary to condense the top half of the y-axis to visualize the 
        central tendency of the data.
    Returns:
        None
    """
    
    num_constructs = len(set(df['Construct']))
    
    pal = sns.color_palette()
    max_repnum = max(df['Replicate'])
    pal = pal[:max_repnum]

    df = pd.DataFrame.from_dict(df)
    high_partco_df = df[df['SG Enrichment Score'] >= ybreakpoint]

    ax = sns.stripplot(x='Construct', y='SG Enrichment Score', data=df, hue='Replicate', palette=pal, alpha=grouped_dotopacity, size=grouped_dotsize, linewidth=0, dodge=True)
    
    plt.xlabel('')

    ax.get_legend().remove()
    
    ax.spines.bottom.set_visible(False)
    ax.set_xticks([])
        
    plt.yticks(fontname='Arial', fontsize=14)
    plt.ylabel('P', fontname='Arial', fontsize=16, color='white')   # ENSURES IDENTICAL FIGURE SIZING FOR FINAL SAVED BOXPLOTS
    if include_title:
        plt.title(data_label, fontname='Arial', fontsize=18)
        
    plt.ylim(ybreakpoint, max(df['SG Enrichment Score']) + (0.07 * (max(df['SG Enrichment Score']) - ybreakpoint)))
        
    fig = plt.gcf()
    fig.set_size_inches(width_scaling*num_constructs, 1)
    plt.savefig(data_label + '_PartitionCoefficients_GroupedBoxplots_SeparatedBoxes_TopHalfOnly.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()
    
    
if __name__ == '__main__':
    main()