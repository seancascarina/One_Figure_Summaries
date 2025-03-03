

def main():

    output = open('RUN_LCD-Composer_Human_SRsearch_Batch.bat', 'w')
    
    for s_comp in range(20, 105, 5):
        for r_comp in range(20, 105, 5):
            if s_comp + r_comp > 100:
                continue
            
            # WRITE COMMANDS FOR S+R SEARCH
            output.write('python LCD-Composer.py UP000005640_9606.fasta Human_S-R_' + str(s_comp) + '-' + str(r_comp) + '_RESULTS -a S_R -c ' + str(s_comp) + '_' + str(r_comp) + '\n')

    output.close()

if __name__ == '__main__':
    main()