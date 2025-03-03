
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    threshold = 70

    rbps = get_rbps()

    dataset = 'S-R'
    pos_aa = dataset[-1]
    df = {}
    prots = {}
    output = open(str(threshold) + 'perc_Combined_' + dataset + 'composition_RBPproportions.csv', 'w')
    output.write(','.join( [pos_aa + ' Comp Threshold -->'] + [str(x) for x in range(20, 85, 5)] ) + '\n')
    output.write('S Comp Threshold\n')
    matrix = []
    for s_comp in range(20, 105, 5):
        proportions = []
        col = []
        for pos_comp in range(20, 105, 5):
            if s_comp + pos_comp > 100:
                col.append(-2)  # DUMMY VALUE USED FOR MASKING CELLS
                continue
                
            comp_str = str(s_comp) + '_' + str(pos_comp)
            prots[comp_str] = set()
            
            h = open('Human_' + dataset + '_' + str(s_comp) + '-' + str(pos_comp) + '_RESULTS.tsv')
            for i in range(7):
                h.readline()
            header = h.readline()
            for line in h:
                items = line.rstrip().split('\t')
                junk, uniprot, *junk = items[0].split('|')
                seq = items[1]
                boundaries = items[2]
                prots[comp_str].add(uniprot)
                
                if s_comp + pos_comp >= threshold:
                    df[uniprot] = df.get(uniprot, [[], '0', [], []])  # value is a list, where the first item in the list is another list containing all of the S/R composition combinations in which the protein was identified. The second item in the list is a binary (0 or 1) string value, where 0 indicates that the protein is not an RBP and 1 indicates that the protein is an RBP. The third item in the list is another list containing each of the domain seqs identified for the corresponding composition criteria. The fourth item in the list is another list containing the domain boundaries for each of the domains identified at each composition criteria.
                    df[uniprot][0].append(comp_str)
                    if uniprot in rbps:
                        df[uniprot][1] = '1'
                    df[uniprot][2].append(seq)
                    df[uniprot][3].append(boundaries)
            h.close()
            
            if len(prots[comp_str]) == 0:
                col.append(-1)  # DUMMY VALUE TO DISTINGUISH BETWEEN CELLS THAT HAVE A TRUE RBP PROPOTION == 0.0 AND CELLS THAT DID NOT HAVE ANY IDENTIFIED SK OR SR PROTEINS
                continue
                
            proportion_rbp = sum( [1 for prot in prots[comp_str] if prot in rbps] ) / len(prots[comp_str])
            proportions.append(str(proportion_rbp))
            col.append( proportion_rbp )
            
        matrix.append(col)
        output.write(','.join( [str(s_comp)] + proportions ) + '\n')
        
    output.close()
    
    plot_heatmap(matrix, threshold, dataset, pos_aa)
    
    output = open('SR_proteins_with_Combined_' + dataset + '_Above_' + str(threshold) + '.tsv', 'w')
    output.write('\t'.join( ['Protein', 'Is RBP?', dataset + ' Composition Thresholds where Protein Identified', 'Domain Sequences Identified at Each Corresponding ' + dataset + ' Composition Threshold', 'Domain Boundaries at Each Corresponding ' + dataset + ' Composition Threshold'] ) + '\n')
    for prot in df:
        output.write('\t'.join( [prot, df[prot][1], ','.join(df[prot][0]), ','.join(df[prot][2]), ','.join(df[prot][3])] ) + '\n')
    output.close()
            

def plot_heatmap(matrix, threshold, dataset, pos_aa):

    matrix = np.matrix.transpose(np.array(matrix))
    colors = sns.color_palette("coolwarm", 1000)

    sns.heatmap(matrix, mask=matrix<-1.5, cmap=['0.5'], annot=False, cbar=False)   # MAKE SOLID GREY HEATMAP TO INDICATE CELLS WITH 0 IDENTIFIED SK OR SR PROTEINS
    sns.heatmap(matrix, mask=matrix<0, cmap=colors, fmt='.2f', annot=True, vmin=0, vmax=1)   # PLOT RBP PROPORTIONS FOR CELLS WITH AT LEAST 1 IDENTIFIED SK OR SR PROTEIN
    plt.xticks([x+0.5 for x in range(17)], labels=[str(x) for x in range(20, 105, 5)])
    plt.yticks([x+0.5 for x in range(17)], labels=[str(x) for x in range(20, 105, 5)], rotation=0)
    plt.xlabel('S Composition Threshold')
    plt.ylabel(pos_aa + ' Composition Threshold')
    plt.ylim(0, 13)
    plt.xlim(0, 13)
    fig = plt.gcf()
    fig.set_size_inches(8,5)
    if dataset == 'S-R':
        plt.savefig('Fig 1B - Human Pairwise ' + dataset +'comp Range_Proportion in RBPs.tiff', bbox_inches='tight', dpi=600)
    else:
        plt.savefig('Fig S1B - Human Pairwise ' + dataset +'comp Range_Proportion in RBPs.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def get_rbps():

    h = open('Complete_list_of_RBPs.txt')
    rbps = []
    for line in h:
        rbps.append(line.rstrip())
        
    return rbps

    
if __name__ == '__main__':
    main()