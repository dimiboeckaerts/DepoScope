"""
Created on 03/03/23

@author: dimiboeckaerts

PhageDEPOdetection DATA PROCESSING TESTS
"""

def test_remove_sequences():
    domains_dict = {'ON337196_00055': ['IPR002925__53.3', 'IPR003925__24.9'], 
                    'JQ015307_00248': ['IPR001724__35.3'], 
                    'MF036692_00040': ['IPR011496__28.9']}
    to_remove = ['IPR002925']
    proteins_to_remove = []
    for proteinid in domains_dict.keys():
        domains_list = [domain.split('__')[0] for domain in domains_dict[proteinid]]
        if any([domain in domains_list for domain in to_remove]):
            proteins_to_remove.append(proteinid)
    assert proteins_to_remove == ['ON337196_00055']
