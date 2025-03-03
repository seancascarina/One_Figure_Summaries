
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def main():

    dataset = 'S-R'
    pos_aa = dataset[-1]
    matrix = []
    
    for s_comp in range(20, 105, 5):
        col = []
        for pos_comp in range(20, 105, 5):
            if s_comp + pos_comp > 100:
                col.append(0.01)    # DUMMY VALUE USED FOR MASKING CELLS
                continue

            total_prots = get_total_prots(s_comp, pos_comp, dataset)
            if total_prots == 0:
                col.append(0.1) # DUMMY VALUE FOR COLOR SCALE ONLY: LATER CONVERTED TO 0
                continue

            col.append(total_prots)
            
        matrix.append(col)

    plot_heatmap(matrix, dataset, pos_aa)
    
    
def plot_heatmap(matrix, dataset, pos_aa):

    matrix = np.matrix.transpose(np.array(matrix))
    mask = np.where(matrix<0.05, True, False)     # MASK CELLS WITH S+R COMPOSITION THRESHOLD EXCEEDING 100
    
    log_mat = np.log10(matrix)  # LOG TRANSFORM MATRIX FOR HEATMAP CELL COLORS
    log_mat = np.where(log_mat<0, 0, log_mat)   # CONVERT DUMMY VALUES TO 0
    
    colors = sns.color_palette("coolwarm", 1000)
    
    sns.heatmap(log_mat, cmap=colors, fmt='.0f', annot=matrix, mask=mask)
    plt.xticks([x+0.5 for x in range(17)], labels=[str(x) for x in range(20, 105, 5)])
    plt.yticks([x+0.5 for x in range(17)], labels=[str(x) for x in range(20, 105, 5)], rotation=0)
    plt.xlabel('S Composition Threshold')
    plt.ylabel(pos_aa + ' Composition Threshold')
    plt.ylim(0, 13)
    plt.xlim(0, 13)
    fig = plt.gcf()
    fig.set_size_inches(8,5)
    plt.savefig('Fig 1A - Human Pairwise ' + dataset + 'comp Range LCD Frequencies.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def get_total_prots(s_comp, pos_comp, dataset):

    h = open('Human_' + dataset + '_' + str(s_comp) + '-' + str(pos_comp) + '_RESULTS.tsv')
    prot_tracker = set()
    for i in range(7):
        h.readline()
    header = h.readline()
    total_prots = 0
    for line in h:
        items = line.rstrip().split('\t')
        junk, uniprot, *junk = items[0].split('|')
        if uniprot in prot_tracker:
            continue
            
        total_prots += 1
        prot_tracker.add(uniprot)
    h.close()

    return total_prots
    
    
if __name__ == '__main__':
    main()
    
