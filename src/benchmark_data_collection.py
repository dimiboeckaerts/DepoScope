"""
DATA COLLECTION FOR BENCHMARKING

@date: 02/09/2023
@author: dimiboeckaerts
"""

from Bio import SeqIO, Entrez
absolute_path = '/Users/Dimi/Desktop/pires_dna_sequences/'

pires_accessions = ['NC_021856.1','NC_023693.1','NC_023693.1','NC_048631.1','NC_019487.1','NC_007636.1','NC_007636.1','NC_024792.1','HQ632825.1','NC_041856.1',
                    'NC_048631.1','NC_027993.1','NC_011048.1','NC_019923.1','NC_027994.1','NC_006883.2','JX262376.1','NC_019414.1','KC700557.1','KC700558.1',
                    'GU323708.1','NC_021802.1','NC_015296.1','NC_018836.1','NC_023610.1','NC_019526.1','KC700556.1','GU323708.1','NC_022768.1','NC_041897.1',
                    'NC_020083.1','KC821630.1','KC821629.1','KC821610.1','KC821634.1','NC_020860.1','NC_023579.1','GU196281.1','GQ413937.1','NC_016071.1',
                    'DQ831957.1','DQ834250.1','NC_018281.1','NC_018284.1','NC_017981.1','NC_022772.1','NC_007055.1','NC_004165.1','NC_021336.1','NC_049976.1',
                    'NC_001423.1','NC_022088.2','NC_019917.1','AY349011.3','NC_021804.1','NC_020842.1','NC_048629.1','NC_019523.1','NC_017980.1','NC_018083.1',
                    'NC_018084.1','NC_019400.1','CP000711.1','NC_008152.1','NC_007637.1','NC_009014.1','NC_019926.1','FQ482084.1','NC_015295.1','NC_027364.1',
                    'NC_019403.1','NC_024378.1','NC_017971.2','NC_015208.1','NC_011976.1','NC_002730.1','NC_027351.1','NC_005344.1','NC_024355.1','NC_023582.1',
                    'FR671411.1','FR671410.1','FR671407.1','FR671406.1','HG799497.1','HG799490.1','HG799496.1','KJ417497.1','FR671405.1','NC_022761.1','NC_002649.1',
                    'FJ230960.1','KC556894.1','NC_019446.1','KC821625.1','KC821621.1','NC_020078.1','NC_008152.1','NC_025446.1','NC_022773.1','KF669651.1','KF669659.1',
                    'NC_020478.1','NC_020479.1','NC_020477.1','KC330683.1','KC330681.1','NC_041858.1','KC821633.1','KC821608.1','NC_019510.1','NC_020079.1','NC_019454.1',
                    'NC_021563.1','NC_014229.1','NC_020873.1','KC821633.1','NC_022769.1','NC_022761.1','NC_004166.2','NC_021792.1','NC_023735.1','NC_022771.1',
                    'NC_022763.1','NC_022767.1','NC_020883.1','NC_023694.1','NC_020079.1','NC_022772.1','NC_022764.1','NC_022770.1','KF669657.1','NC_019530.1',
                    'NC_020083.1','NC_020081.2','NC_019929.1','NC_022761.1','NC_019401.1','NC_003157.5','NC_023557.1','NC_023501.1','NC_018857.1','NC_025422.1',
                    'NC_022761.1','NC_024213.1','NC_024211.1','NC_019487.1','NC_009819.1','NC_024216.1','NC_024205.1','NC_018856.1','NC_020083.1','NC_020081.2',
                    'NC_021856.1','NC_017972.1']
pires_accessions = list(set(pires_accessions)) # remove dups

# download the records from NCBI
Entrez.email= 'dimitri.boeckaerts@ugent.be'
records_gb = Entrez.efetch(db='nucleotide', id=pires_accessions, rettype='gbwithparts', retmode='text') 

# loop over the records and save their CDS's in FASTA files
for record in SeqIO.parse(records_gb, 'gb'):
    file = open(absolute_path+'CDS_' + record.id + '.fasta', 'w')
    if record.features:
        for feature in record.features:
            if feature.type == "CDS":
                try:
                    feature_cds = str(feature.location.extract(record).seq)
                    feature_id = feature.qualifiers['locus_tag'][0]
                    feature_product = feature.qualifiers['product'][0]
                    feature_translation = feature.qualifiers['translation'][0]
                    file.write('>' + feature_id + '_' + feature_product + '\n' + feature_cds + '\n')
                except KeyError:
                    pass
    file.close()
    print()
    print(record.id + ' done')