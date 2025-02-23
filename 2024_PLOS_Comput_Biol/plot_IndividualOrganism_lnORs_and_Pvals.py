
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import math

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_strings = [aa for aa in amino_acids]
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)

domains = ['Archaea','Bacteria','Eukaryota','Viruses']
proteomes_of_interest = ['UP000005640_9606', 'UP000001450_36329']
fig_labels = {('UP000001450_36329', 'lnOR'):'Fig5A',
            ('UP000001450_36329', 'Statistical Significance'):'Fig5B',
            ('UP000005640_9606', 'lnOR'):'Fig5C',
            ('UP000005640_9606', 'Statistical Significance'):'Fig5D'
            }

def main():

    domain = 'Eukaryota'
    aa_freqs_df, ordered_aas_df = determine_aa_order()

    for proteome_of_interest in proteomes_of_interest:
        pvals_df, lnOR_df, proteome_pvals_df, datalines_df, header_line, pvals_list, lcd_class_list = get_data(proteome_of_interest)

        # THIS IS USING THE PRE-CALCULATED HOLM-SIDAK CORRECTED P-VALUES
        corrected_pvals_df = {lcd_class:pvals_list[i] for i, lcd_class in enumerate(lcd_class_list)}

        matrix = []
        ordered_aas = ordered_aas_df[proteome_of_interest]
        df = {aa:[] for aa in ordered_aas}
        lnORs_master = {aa:[] for aa in ordered_aas}
        for res1 in ordered_aas:
            row = []
            for res2 in ordered_aas:
                if res1 == res2:
                    pvals = pvals_df[domain][res1]
                    corrected_pval = corrected_pvals_df[res1]
                else:
                    pvals = pvals_df[domain][res1+res2]
                    corrected_pval = corrected_pvals_df[res1+res2]

                pvals = [float(x) for x in pvals if x !='N/A']

                # USE THIS WHEN USING THE lnOR FOR THE HEATMAP============
                if res1 == res2:
                    lnOR = lnOR_df[proteome_of_interest][res1]
                else:
                    lnOR = lnOR_df[proteome_of_interest][res1+res2]

                # USE FOR CORRECTED LOG PVALS
                row.append(corrected_pval)
                df[res1].append(corrected_pval)
                
                lnORs_master[res1].append(lnOR)

            matrix.append(row)
            
        binary_df = {}
        for i, res in enumerate(list(df.keys())):
            binary_df[res] = []
            for j in range(len(df[res])):
                val = df[res][j]

                if val == 'N/A':
                    binary_df[res].append(1)    # MASKED VALUE: NO LCDs
                elif val < 0.05:
                    binary_df[res].append(3)
                else:
                    binary_df[res].append(2)
                    
        df = {lcd_class:[-math.log10(df[lcd_class][i]) if df[lcd_class][i] != 'N/A' else 0.0 for i in range(len(df[lcd_class]))] for lcd_class in df}
        
        df = pd.DataFrame.from_dict(df)
        df['Secondary Amino Acid'] = list(ordered_aas)
        df.set_index('Secondary Amino Acid', inplace=True)
        df.to_csv(proteome_of_interest + '_SecondaryLCDs_Matrix_Log10correctedPvals.tsv', sep='\t')

        binary_df = pd.DataFrame.from_dict(binary_df)
        binary_df['Secondary Amino Acid'] = list(ordered_aas)
        binary_df.set_index('Secondary Amino Acid', inplace=True)
        fig_label = fig_labels[ (proteome_of_interest, 'Statistical Significance') ]
        plot_binary_heatmap(binary_df, proteome_of_interest, fig_label)
        binary_df.to_csv(proteome_of_interest + '_SecondaryLCDs_HeatmapMatrix_Log10correctedPvals_BINARY_TABLE.tsv', sep='\t')

        lnORs_master = pd.DataFrame.from_dict(lnORs_master)
        lnORs_master['Secondary Amino Acid'] = list(ordered_aas)
        lnORs_master.set_index('Secondary Amino Acid', inplace=True)
        fig_label = fig_labels[ (proteome_of_interest, 'lnOR') ]
        plot_heatmap(lnORs_master, proteome_of_interest, fig_label)
        lnORs_master.to_csv(proteome_of_interest + '_SecondaryLCDs_HeatmapMatrix_lnORs.tsv', sep='\t')
        

def get_data(proteome_of_interest):

    h = open('Observed_vs_Scrambled_FisherExact_Results.tsv')
    header = h.readline()
    df = {domain:{lcd_class:[] for lcd_class in aa_strings} for domain in domains}
    lnOR_df = {}
    proteome_pvals_df = {}
    lcd_freqs_df = {}
    datalines_df = {}
    pvals_list = []
    lcd_class_list = []
    
    for line in h:
        items = line.rstrip().split('\t')
        proteome = items[0]
        if proteome != proteome_of_interest:
            continue
        lnOR_df[proteome] = lnOR_df.get(proteome, {})
        domain = items[1]
        lcd_class = items[2]
        pval = items[12]    # SIDAK-HOLM CORRECTED P-VAL
        obs_num = int(items[3])
        scr_num = int(items[4])
        total_prots = int(items[5])
        
        datalines_df[lcd_class] = line

        lnOR = items[7]
        if lnOR == 'N/A':
            if items[14] == 'N/A':
                print('\nSome log odds were not calculable even with biased estimates. Analysis could not be completed for this organism. This only occurs when all of the proteins in a proteome contain an LCD of one type, and is believed to occur only for small viruses with a total of 1-2 proteins in the proteome. Exiting program.\n')
                exit()
            lnOR = float(items[14])     # USE BIASED lnOR IF THE ACTUAL lnOR IS NOT AVAILABLE
        else:
            lnOR = float(lnOR)
        if pval != 'N/A':
            pval = float(pval)
        lnOR_df[proteome][lcd_class] = lnOR

        if pval == 0.0:
            pval = 10**-300

        df[domain][lcd_class].append(pval)
        pvals_list.append(pval)
        lcd_class_list.append(lcd_class)
        
    h.close()
    
    return df, lnOR_df, proteome_pvals_df, datalines_df, header, pvals_list, lcd_class_list
    

def plot_heatmap(df, proteome_of_interest, fig_label):

    pal = sns.color_palette('gist_heat', as_cmap=True)
    
    cg = sns.heatmap(df, cmap=pal)
    ax = plt.gca()
    ax.set_facecolor("grey")
    cb=ax.collections[0].colorbar
    cb.set_ticklabels(cb.ax.get_yticklabels(), fontsize=10, fontname='Arial') 
    
    plt.xticks(rotation=0, fontname='Arial', fontsize=10)
    plt.yticks(rotation=0, fontname='Arial', fontsize=10)
    plt.xlabel('Primary Amino Acid', fontname='Arial', fontsize=12)
    plt.ylabel('Secondary Amino Acid', fontname='Arial', fontsize=12)
    
    plt.ylim(20, 0)
    fig = plt.gcf()
    plt.savefig(fig_label + '_' + proteome_of_interest + '_Original-vs-Scrambled-Proteome_Heatmap_lnORvalues.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()
    
    
def plot_binary_heatmap(df, proteome_of_interest, fig_label):

    pal = ['grey', 'black', 'red']

    cg = sns.heatmap(df, cmap=pal)
    ax = plt.gca()
    ax.set_facecolor("white")
    cb=ax.collections[0].colorbar

    cb.set_ticklabels(cb.ax.get_yticklabels(), fontsize=10, fontname='Arial') 
    
    plt.xticks(rotation=0, fontname='Arial', fontsize=10)
    plt.yticks(rotation=0, fontname='Arial', fontsize=10)
    plt.xlabel('Primary Amino Acid', fontname='Arial', fontsize=12)
    plt.ylabel('Secondary Amino Acid', fontname='Arial', fontsize=12)
    
    plt.ylim(20, 0)
    fig = plt.gcf()
    plt.savefig(fig_label + '_' + proteome_of_interest + '_Original-vs-Scrambled-Proteome_Heatmap_Pvalues_BinaryStatSigClassification.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()
    
    
def determine_aa_order():

    h = open('Background_AAfrequencies_AllProteomes.tsv')
    header = h.readline()
    df = {}
    aa_orders_df = {}
    for line in h:
        domain, proteome, *aa_freqs = line.rstrip().split('\t')
        aa_freqs = [int(x) for x in aa_freqs]
        df[proteome] = aa_freqs
        
        aas = list(amino_acids)
        aa_freqs, aas = zip(*sorted(zip(aa_freqs, aas), reverse=True))
        aa_orders_df[proteome] = aas
        
    h.close()
    
    return df, aa_orders_df
        
if __name__ == '__main__':
    main()
    