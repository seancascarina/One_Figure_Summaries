
from scipy.stats import fisher_exact
import math
import mpmath

# CREATE LIST OF STRINGS FOR ALL 400 LCD CLASSES
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_strings = [aa for aa in amino_acids]
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)
        
def main():

    # GET LCD FREQUENCY DATA
    obs_df, org_to_domain = get_frequencies('TableS1_LCDfrequency_NumberOfProtsWithLCDs_Humans_Malaria_Only.tsv')
    scr_df, org_to_domain_scr = get_frequencies('TableS1_LCDfrequency_NumberOfProtsWithLCDs_SCRAMBLED_Humans_Malaria_Only.tsv')

    # PREP OUTPUT FILE
    output = open('Observed_vs_Scrambled_FisherExact_Results.tsv', 'w')
    output.write('\t'.join(['Proteome', 'Domain of Life', 'LCD Class', '# of Proteins with LCDs, Actual Proteome ("Observed")', '# of Proteins with LCDs, Scrambled Proteome', 'Total Proteins in Proteome', 'OddsRatio', 'lnOR', 'Fold Change (# in actual proteome / # in scrambled proteome)', '95% Confidence Interval (lower bound, upper bound)', '95% Confidence Interval for Odds Ratio Excludes 1?', 'p-value', 'Sidak-Holm Corrected p-value', 'Biased OddsRatio', 'Biased lnOR (when necessary)', 'Biased Fold Change [(# in actual proteome + 1) / (# in scrambled proteome + 1)]', 'Biased 95% Confidence Interval (lower bound, upper bound)', 'Biased Raw p-value (when necessary)']) + '\n')
    
    # LOOP OVER PROTEOMES TO MAKE PLOTS FOR
    for proteome in obs_df:
        scr_label = proteome + '_SCRAMBLED'
        if scr_label in scr_df:
            pvals = []
            lcd_classes = []

            data_df = {}
            biased_df = {}
            domain = org_to_domain[proteome]
            
            # LOOP OVER 400 LCD CLASSES
            for lcd_class in aa_strings:
                obs = obs_df[proteome][lcd_class]
                scr = scr_df[scr_label][lcd_class]
                total_prots_obs = obs_df[proteome]['Total Proteins']
                total_prots_scr = scr_df[scr_label]['Total Proteins']
                
                # INITIALIZE BIAS VARIABLES AS N/A (CASES IN WHICH BIASED ESTIMATES ARE NOT NECESSARY)
                biased_oddsratio, biased_lnOR, biased_upper_CI, biased_lower_CI, biased_pval = ['N/A']*5
                
                if obs == 0 and scr == 0:   # NEED BIASED ESTIMATES DUE TO ZEROS
                    data_line = [str(x) for x in (proteome, domain, lcd_class, obs, scr, total_prots_obs)] + ['N/A']*6
                    data_df[lcd_class] = data_line
                    biased_oddsratio, biased_lnOR, biased_upper_CI, biased_lower_CI, biased_pval = calc_lnOR(obs+1, total_prots_obs+1, scr+1, total_prots_scr+1, proteome)
                    biased_fold_change = (obs+1) / (scr+1)
                    if biased_oddsratio == 'N/A':
                        biased_ci = 'N/A'
                    else:
                        biased_ci = (biased_lower_CI, biased_upper_CI)
                
                    biased_df[lcd_class] = [str(x) for x in (biased_oddsratio, biased_lnOR, biased_fold_change, biased_ci, biased_pval)]

                    continue

                # CALCULATE ODDS RATIO AND CONFIDENCE INTERVALS
                oddsratio, lnOR, upper_CI, lower_CI, pval = calc_lnOR(obs, total_prots_obs, scr, total_prots_scr, proteome)

                if upper_CI != 'N/A' and lower_CI != 'N/A':
                    is_in_CI = 0
                    if lower_CI <= 0 <= upper_CI:
                        is_in_CI = 1
                else:
                    is_in_CI = 'N/A'

                if obs == 0 or scr == 0:    # FOR BIASED ESTIMATES
                    fold_change = 'N/A'
                    biased_oddsratio, biased_lnOR, biased_upper_CI, biased_lower_CI, biased_pval = calc_lnOR(obs+1, total_prots_obs+1, scr+1, total_prots_scr+1, proteome)
                    biased_fold_change = (obs+1) / (scr+1)
                else:   # FOR NON-BIASED ESTIMATES
                    fold_change = obs / scr
                    biased_fold_change = 'N/A'

                if oddsratio == 'N/A':
                    ci = 'N/A'
                else:
                    ci = (lower_CI, upper_CI)
                if biased_oddsratio == 'N/A':
                    biased_ci = 'N/A'
                else:
                    biased_ci = (biased_lower_CI, biased_upper_CI)
                    
                # CREATE DATA STRING FOR OUTPUT
                data_line = [str(x) for x in (proteome, domain, lcd_class, obs, scr, total_prots_obs, oddsratio, lnOR, fold_change, ci, is_in_CI, pval)]

                data_df[lcd_class] = data_line
                biased_df[lcd_class] = [str(x) for x in (biased_oddsratio, biased_lnOR, biased_fold_change, biased_ci, biased_pval)]
                
                if pval != 'N/A':
                    pvals.append(pval)
                    lcd_classes.append(lcd_class)

            # CORRECT P-VALUES FOR MULTIPLE HYPOTHESIS TESTS
            if len(pvals) == 0:
                corrected_pvals = []
            else:
                corrected_pvals = sidak_correction(pvals)
            
            # PREPARE DATA FOR OUTPUT
            for lcd_class in aa_strings:
                if lcd_class in lcd_classes:
                    index = lcd_classes.index(lcd_class)
                    corrected_pval = corrected_pvals[index]
                    data_df[lcd_class].append(str(corrected_pval))
                else:
                    data_df[lcd_class].append('N/A')
                    
                data_df[lcd_class] += biased_df[lcd_class]
                output.write('\t'.join(data_df[lcd_class]) + '\n')
                
    output.close()


def sidak_correction(pvals):
    """Performs multiple test correction using Sidak-Holm method.
    Returns:
        final_pvals = list of corrected p-values, sorted in ascending order.
    """
    
    positions = [i for i in range(len(pvals))]
    sorted_pvals, pval_positions = zip(*sorted(zip(pvals, positions)))
    m = len(pvals)

    mpmath.mp.dps = 300

    corrected_pvals = [1-(1-mpmath.mp.mpf(sorted_pvals[len(pvals)-i]))**(i) for i in range(m, 0, -1)]

    final_pvals = [corrected_pvals[0]]
    for i in range(1, len(corrected_pvals)):
        pval = max(final_pvals[-1], corrected_pvals[i])
        final_pvals.append(pval)
        
    corrected_pvals = [corrected_pvals[0]] + [max(corrected_pvals[i-1], corrected_pvals[i]) for i in range(1, len(corrected_pvals))]
    corrected_pvals = [mpmath.nstr(x, 16) for x in corrected_pvals]
    final_pvals = [mpmath.nstr(x, 16) for x in final_pvals]

    pval_positions, final_pvals = zip(*sorted(zip(pval_positions, final_pvals)))

    return final_pvals
    
         
def calc_lnOR(obs, total_prots_obs, scr, total_prots_scr, proteome):
    """Calculate the natural log of the odds ratio for the LCD frequency for a given
    LCD class in the original proteome versus scrambled proteome.
    Returns:
        oddsratio = float
        lnOR = natural log of the odds ratio, float
        upper_CI = upper bound on the 95% confidence interval, float
        lower_CI = lower bound on the 95% confidence interval, float
        pval = p-value from Fisher's exact test, float
    """
    contingency = [ [obs, (total_prots_obs-obs)], [scr, (total_prots_scr-scr)] ]

    oddsratio, pval = fisher_exact(contingency, alternative='two-sided')

    # lnOR, oddsratio, AND CONFIDENCE INTERVALS CAN'T BE CALCULATED IF 
    # OBSERVED OR SCRAMBLED FREQUENCIES ARE ZERO. REQUIRES BIASED ESTIMATION, 
    # WHICH IS PERFORMED LATER.
    if obs == 0 or scr == 0 or obs == total_prots_obs or scr == total_prots_scr:
        lnOR = 'N/A'
        upper_CI = 'N/A'
        lower_CI = 'N/A'
        oddsratio = 'N/A'

    # CASES WHERE BIASED ESTIMATION IS NOT NEEDED
    else:
        lnOR = math.log(oddsratio)
        
        upper_CI = lnOR + 1.96 * math.sqrt(1/obs + 1/total_prots_obs + 1/scr + 1/total_prots_scr)
        lower_CI = lnOR - 1.96 * math.sqrt(1/obs + 1/total_prots_obs + 1/scr + 1/total_prots_scr)

    return oddsratio, lnOR, upper_CI, lower_CI, pval
    

def get_frequencies(file):
    """Get LCD frequencies for each LCD class and organism.
    Returns:
        df = multi-dimensional dictionary
                1st dimension ---> proteome as keys, LCD classes as values
                2nd dimension ---> LCD classes as keys, # of proteins with that type of LCD as values
                    This dimension also contains "Total Proteins" key and corresponding # of proteins in the proteome as the value.
        org_to_domain = mapping dictionary with organism ID as keys and corresponding domain of life as values
    """
    
    h = open(file)
    header = h.readline()
    df = {}
    org_to_domain = {}
    
    for line in h:
        items = line.rstrip().split('\t')
        proteome = items[0]
        domain = items[2]
        total_prots = int(items[5])
        lcd_freqs = [int(x) for x in items[6:]]
        df[proteome] = df.get(proteome, {})
        for i, lcd_class in enumerate(aa_strings):
            df[proteome][lcd_class] = lcd_freqs[i]
        df[proteome]['Total Proteins'] = total_prots
        org_to_domain[proteome] = domain
        
    h.close()
        
    return df, org_to_domain


if __name__ == '__main__':
    main()
