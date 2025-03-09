
from Bio import Entrez
from tqdm import tqdm

Entrez.email = ''    # USER WOULD NEED TO INPUT THEIR OWN EMAIL
my_key = ''         # USER WOULD NEED TO INPUT THEIR OWN PUBMED API KEY

def main():

    journals = ['PLoS Pathog', 'PLoS Biol', 'PLoS Genet', 'PLoS Comput Biol', 'PLoS Med', 'Cell', 'Nature', 'Science', 'Genetics', 'Bioinformatics', 'BMC Med', 'mBio', 'J Cell Biol', 'Curr Biol', 'Mol Cell Biol']
    journal_args = ['PlosPath', 'PlosBiol', 'PlosGenet', 'PlosCompBiol', 'PlosMed', 'Cell', 'Nature', 'Science', 'Genetics', 'Bioinformatics', 'BMCmed', 'mBio', 'JCellBiol', 'CurrBiol', 'MCB']
    journal_df = {journal_arg:journals[i] for i, journal_arg in enumerate(journal_args)}
    
    for journal_arg in journal_args:
        journal = journal_df[journal_arg]
        
        max_tries = 10

        pmids = get_pubmed_ids(journal)
        all_pmids = sorted([x for x in pmids])

        main_badrequests_output = open('FailedQueries_for_MainArticle_PubmedIDs_' + journal_arg + '.txt', 'w')
        
        refs_badrequests_output = open('FailedQueries_for_RefArticle_PubmedIDs_' + journal_arg + '.tsv', 'w')
        refs_badrequests_output.write('MainAritcle PubmedID\tRefArticle PubmedID\n')
        
        output = open('Pubmed_AuthorSearch_' + journal_arg + '.tsv', 'w')
        output.write('\t'.join(['Main Article PubmedID', 'Main Article Combinednames', 'Main Article Lastnames', 'Main Article Firstnames', 'Main Article Initials', 'Main Dictionary Contained Firstname Key?', 'Cited Article PubmedID', 'Cited Article Combinednames', 'Cited Article Lastnames', 'Cited Article Firstnames', 'Cited Article Initials', 'Cited Dictionary Contained Firstname Key?']) + '\n')
        for main_pmid in tqdm(all_pmids):
            successful_main = False
            tries = 0
            while not successful_main and tries < max_tries:
                try:
                    if tries > 0:
                        print(tries, main_pmid)
                    tries += 1
                    main_art_combinednames, main_art_lastnames, main_art_firstnames, main_art_initials, main_art_containedfirstnames = get_authorlist(main_pmid)
                    
                    results = Entrez.read(Entrez.elink(db="pubmed", LinkName="pubmed_pubmed_refs", from_uid=main_pmid, api_key=my_key))
                    id_list = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]
                    successful_main = True

                except:
                    continue
                    
            if tries == max_tries:
                main_badrequests_output.write(main_pmid + '\n')
                continue

            for ref_pmid in id_list:
                successful_ref = False
                iternum = 0
                while not successful_ref and iternum < max_tries:
                    try:
                        if iternum != 0:
                            print(iternum, main_pmid, ref_pmid)
                        iternum += 1
                        ref_art_combinednames, ref_art_lastnames, ref_art_firstnames, ref_art_initials, ref_art_containedfirstnames = get_authorlist(ref_pmid)
                        
                        output.write('\t'.join([main_pmid, ';'.join(main_art_combinednames), ';'.join(main_art_lastnames), ';'.join(main_art_firstnames), ';'.join(main_art_initials), main_art_containedfirstnames, ref_pmid, ';'.join(ref_art_combinednames), ';'.join(ref_art_lastnames), ';'.join(ref_art_firstnames), ';'.join(ref_art_initials), ref_art_containedfirstnames]) + '\n')
                        successful_ref = True

                    except:
                        continue
                        
                if iternum == max_tries:
                    refs_badrequests_output.write('\t'.join([main_pmid, ref_pmid]) + '\n')
                    continue

        main_badrequests_output.close()
        refs_badrequests_output.close()
        output.close()
        
        
def get_authorlist(pmid):

    handle = Entrez.efetch(db="pubmed", id=pmid, retmax='1', retmode="xml", api_key=my_key)
    results = Entrez.read(handle)
    authors = results['PubmedArticle'][0]['MedlineCitation']['Article']['AuthorList']

    lastnames = []
    firstnames = []
    initials = []
    combinednames = []
    containedfirstnames = '0'
    for author in authors:
        lastname = author['LastName']
        if 'ForeName' in author:
            firstname = author['ForeName']  # CAN INCLUDE A MIDDLE INITIAL. FORM WILL BE THE FIRSTNAME + SPACE + MIDDLE INITIAL.
            containedfirstnames = '1'
        elif 'Initials' in author:
            firstname = author['Initials']
        else:
            firstname = ''
            
        if 'Initials' in author:
            initial = author['Initials']
        else:
            initial = ''

        lastnames.append(lastname)
        firstnames.append(firstname)
        initials.append(initial)
        combinednames.append(firstname + ' ' + lastname)
        
    return combinednames, lastnames, firstnames, initials, containedfirstnames
    
    
def get_pubmed_ids(journal):

    h = open('Filtered_PMIDs_AllJournals.csv')
    df = {}
    
    total = 0
    total_to_analyze = 0
    header = h.readline()
    for line in h:
        items = line.rstrip().split(',')
        if items[0] != journal:
            continue
        pmid = items[9]
        if pmid == '':
            total += 1
            continue
        df[pmid] = items
        
        total += 1
        total_to_analyze += 1
        
    h.close()

    return df
    
        
if __name__ == '__main__':
    main()