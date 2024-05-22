"""
DATA COLLECTION FOR BENCHMARKING

@date: 02/09/2023
@author: dimiboeckaerts
"""

import pandas as pd
from Bio import SeqIO, Entrez
absolute_path = '/Users/dimi/Documents/GitHub/PhageDEPOdetection/data/Benchmarking/'

# 1 - NCBI DATA COLLECTION, FASTAS, DATAFRAME
# ------------------------------------------
# collect the depolymerase proteins from NCBI
pires_depos = ['HM214492', 'YP_009056510', 'AGY46593', 'AGY46829', 'AGY46940', 'AGY46975', 'AGY46977', 'AGY46980', 'AGY47215', 'AGY47298', 'AGY47935', 'AGB62649',
               'YP_007349243', 'AGY48254', 'AGY48030', 'AGY48191', 'AGY48483', 'AGY48730', 'YP_007517493', 'Q37893', 'YP_006907302', 'YP_009031349', 'YP_009003154', 
               'YP_006907824', 'YP_009036975', 'YP_007517572', 'YP_007517417', 'YP_007517647', 'NP_073695.1', 'AGE60867', 'YP_009036526', 'YP_009035268', 'YP_009036150',
               'YP_008060120', 'ACH57080', 'NP_690774', 'ADF59152', 'ACE96035', 'YP_007678034', 'AAA88489', 'YP_007003341', 'YP_007003328', 'YP_002300374', 'NP_690719', 
               'NP_690860', 'AGE60945', 'YP_008430873', 'YP_008318429', 'YP_008318292', 'YP_007677121', 'AHB81145', 'AEY69698', 'YP_007236804', 'NP_944291', 'AGS80957', 
               'YP_008242013', 'AGO49035', 'AGO49697', 'AGO49696', 'AGO48665', 'YP_008240873', 'AGO49328', 'AGO49248', 'YP_008242178', 'AGO47740', 'YP_008241124', 'AGO49746',
               'YP_007675684', 'YP_007673458', 'ADR30483', 'AEW47853', 'YP_006383591', 'YP_006488624', 'YP_006488654', 'YP_006987059', 'YP_006987233', 'YP_007348361', 'ABQ88383',
               'YP_654148', 'YP_654147', 'CAJ29458', 'YP_338127', 'CAJ29390', 'CBY99572', 'CBY99579', 'YP_007007685', 'YP_009004869', 'YP_009007463', 'YP_001039683', 'YP_005098420', 
               'CBX45113', 'CBX44510', 'YP_009010167', 'AFQ96603', 'YP_004327331', 'AEJ81510', 'AFO10349', 'ADA82273', 'ADA82374', 'ADA82322', 'ADA82474', 'AGC35227', 'YP_007348539',
               'YP_007348546', 'YP_006987816', 'YP_003347555', 'ADA79896', 'ADA79897', 'YP_009044242', 'YP_007002901', 'YP_214529', 'AGN12301', 'AFH14728', 'CCE60816', 'YP_004286222',
               'YP_007237194', 'ADD80892', 'ADD80999', 'ABA54611', 'NP_112090', 'AGY47555', 'AGY47556', 'AGY47760', 'YP_004893855', 'AFU63680', 'YP_004327544', 'YP_007008117', 
               'YP_008130362', 'AFX93502', 'AFX93505', 'AFX93507', 'AAQ12204', 'YP_240092', 'AIA64067', 'YP_950618', 'YP_950681', 'AFM73722', 'AFM73788', 'AHG23941', 'CBW39235',
               'CBW39181', 'CBW39031', 'CBW38974', 'YP_009043178', 'YP_009043945', 'YP_009042723', 'ABL61072', 'AAL15086', 'AHN84645', 'CBW38917', 'YP_008050881', 'AFV51346', 'AFO10889',
               'YP_006906201', 'YP_003714746', 'YP_006990138', 'YP_008051427', 'YP_008060258', 'AEX65697']
Entrez.email= 'dimitri.boeckaerts@ugent.be'
records_gb = Entrez.efetch(db='protein', id=pires_depos, rettype='gbwithparts', retmode='text')

depo_ids = []
depo_seqs = []
for record in SeqIO.parse(records_gb, 'gb'):
    depo_seqs.append(str(record.seq))
    depo_ids.append(record.id)

# loop over the genome records to collect the other sequences
pires_accessions = ['NC_021856.1','NC_023693.1','NC_023693.1','NC_004167.1','NC_019487.1','NC_007636.1','NC_007636.1','NC_024792.1','HQ632825.1','NC_041856.1',
                    'NC_048631.1','NC_027993.1','NC_011048.1','NC_019923.1','NC_027994.1','NC_006883.2','JX262376.1','NC_019414.1','KC700557.1','KC700558.1',
                    'GU323708.1','NC_021802.1','NC_015296.1','NC_018836.1','NC_023610.1','NC_019526.1','KC700556.1','GU323708.1','NC_022768.1','NC_041897.1',
                    'NC_020083.1','KC821630.1','KC821629.1','KC821610.1','KC821634.1','NC_020860.1','NC_023579.1','GU196281.1','GQ413937.1','NC_016071.1',
                    'DQ831957.1','DQ834250.1','NC_018281.1','NC_018284.1','NC_017981.1','NC_022772.1','NC_007055.1','NC_004165.1','NC_021336.1','NC_049976.1',
                    'NC_001423.1','NC_022088.2','NC_019917.1','AY349011.3','NC_021804.1','NC_020842.1','NC_048629.1','NC_019523.1','NC_017980.1','NC_018083.1',
                    'NC_018084.1','NC_019400.1','CP000711.1','NC_008152.1','NC_007637.1','NC_009014.1','NC_019926.1','FQ482084.1','NC_015295.1','NC_027364.1',
                    'NC_019403.1','NC_024378.1','NC_017971.2','NC_015208.1','NC_011976.1','NC_002730.1','NC_027351.1','NC_005344.1','NC_024355.1','NC_023582.1',
                    'FR671411.1','FR671410.1','FR671407.1','FR671406.1','HG799497.1','HG799490.1','HG799496.1','KJ417497.1','FR671405.1','NC_022761.1','NC_002649.1',
                    'NC_011421.1','KC556894.1','NC_019446.1','KC821625.1','KC821621.1','NC_020078.1','NC_008152.1','NC_025446.1','NC_022773.1','KF669651.1','KF669659.1',
                    'NC_020478.1','NC_020479.1','NC_020477.1','KC330683.1','KC330681.1','NC_041858.1','KC821633.1','KC821608.1','NC_019510.1','NC_020079.1','NC_019454.1',
                    'NC_021563.1','NC_014229.1','NC_020873.1','KC821633.1','NC_022769.1','NC_022761.1','NC_004166.2','NC_021792.1','NC_023735.1','NC_022771.1',
                    'NC_022763.1','NC_022767.1','NC_020883.1','NC_023694.1','NC_020079.1','NC_022772.1','NC_022764.1','NC_022770.1','KF669657.1','NC_019530.1',
                    'NC_020083.1','NC_020081.2','NC_019929.1','NC_022761.1','NC_019401.1','NC_003157.5','NC_023557.1','NC_023501.1','NC_018857.1','NC_025422.1',
                    'NC_022761.1','NC_024213.1','NC_024211.1','NC_019487.1','NC_009819.1','NC_024216.1','NC_024205.1','NC_018856.1','NC_020083.1','NC_020081.2',
                    'NC_021856.1','NC_017972.1', 'NC_016767.1', 'NC_024137.1']
pires_accessions = list(set(pires_accessions)) # remove dups
Entrez.email= 'dimitri.boeckaerts@ugent.be'
records_gb = Entrez.efetch(db='nucleotide', id=pires_accessions, rettype='gbwithparts', retmode='text')

# loop over the records and save their CDS's in FASTA files
genome_ids = []
protein_ids = []
dna_seqs = []
protein_seqs = []
for record in SeqIO.parse(records_gb, 'gb'):
    file = open(absolute_path+'CDS_' + record.id + '.fasta', 'w')
    for feature in record.features:
        if feature.type == "CDS":
            feature_cds = str(feature.location.extract(record).seq)
            try:
                feature_id = feature.qualifiers['protein_id'][0]
                protein = feature.qualifiers['translation'][0]
                if protein in depo_seqs: # check sequences, even though perhaps ids dont match
                    feature_id = depo_ids[depo_seqs.index(protein)]
                    depo_seqs.remove(protein)
                    depo_ids.remove(feature_id)
                genome_ids.append(record.id)
                protein_ids.append(feature_id.split('.')[0])
                dna_seqs.append(feature_cds)
                protein_seqs.append(protein)
            except KeyError:
                try:
                    feature_id = feature.qualifiers['locus_tag'][0]
                    protein = feature.qualifiers['translation'][0]
                    if protein in depo_seqs: # check sequences, even though perhaps ids dont match
                        feature_id = depo_ids[depo_seqs.index(protein)]
                        depo_seqs.remove(protein)
                        depo_ids.remove(feature_id)
                    genome_ids.append(record.id)
                    protein_ids.append(feature_id.split('.')[0])
                    dna_seqs.append(feature_cds)
                    protein_seqs.append(protein)
                except KeyError:
                    pass
    file.close()
    print()
    print(record.id + ' done')
print('all done')

labels = [1 if x in pires_depos else 0 for x in protein_ids]
benchmark_df = pd.DataFrame({'genome_id': genome_ids, 'protein_id': protein_ids, 'dna_seq': dna_seqs, 'protein_seq': protein_seqs, 'label': labels})
benchmark_df.to_csv(absolute_path+'benchmark_dataframe.csv', index=False)
not_found = [x for x in pires_depos if x not in protein_ids]
"""
7 proteins were not found: ['HM214492', 'YP_009031349', 'NP_073695.1', 'NP_690774', 'NP_690860', 'CAJ29390', 'YP_005098420']
We will check these manually.
YP_005098420 -> in 'NC_016767.1', which was not in the list
CAJ29390 is actually a duplicate of YP_338127
NP_690860 is removed from NCBI
NP_690774 is not found although another protein from the same alleged genome is found; reality, there are two genomes: phi105 and phi-105. Both contain one of the proteins...
NP_073695 is in the protein_ids list but NP_073695.1 is not; strange. This is probably due to small sequence diffs. Should not affect the outcome.
YP_009031349 was not in list -> NC_024137.1 (was added now)
HM214492 is a 'colanic acid-degrading protein gene'; without associated genome. We will not include this one.
"""