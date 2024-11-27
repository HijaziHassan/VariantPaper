#!/usr/bin/env python3

"""
Author: Hassan Hijazi (https://github.com/HijaziHassan)
This code generates peptidoforms of Lys27-Arg40 stretch of histone H3 and its varaints H33, H3mm7, and H3mm13.
The script was used to generate Figure 1 of the article: 
"Mind Your Spectra: Points to be Aware of when Validating the Identification of Isobaric Histone Peptidoforms": 
https://www.biorxiv.org/content/10.1101/2024.11.13.623238v1.full

For more information on the parameters please refer to the documentation of pyteomics library:
https://pyteomics.readthedocs.io/en/latest/mass.html
https://pyteomics.readthedocs.io/en/latest/api/parser.html#pyteomics.parser.isoforms

"""

import pyteomics.mass as mass
import pyteomics.parser as parser
import pandas as pd


output_file_name = "peptidoforms.csv"

#Define the list of sequences
sequences = ['KSAPATGGVKKPHR', 'KSAPSTGGVKKPHR', 'KSAPSIGGVKKPHR', 'KSVPSTGGVKKPHR']

charge= [0]

ion_type = 'M' #precursor ion

# link mod to amino acid(s) it modifies {mod:[a.a]}
variable_mods={'pr-': True, # N-term 
               'ac': ['K'],
               #'me1': ['K'],
               'me2': ['K'],
               'me3': ['K'],
               'bu': ['K'], # bu = me1 + pr
               'hb': ['K'], 
               'cr': ['K'],
               'la': ['K'],
               'pr': ['S', 'K', 'T']
               #'fo': ['S', 'K', 'T']
              }

# shorthand mod and its corresponding unimod "title" to extract composition 
modifications = {
    'ac': 'Acetyl',
    #'me1': 'Methyl',
    'me2': 'Dimethyl',
    'me3': 'Trimethyl',
    'bu': 'Butyryl',
    'hb': 'hydroxyisobutyryl',
    'cr': 'Crotonyl',
    'pr': 'Propionyl',
    'la': 'Lactylation'
    #'fo': 'Formyl',
}

# standard amino acid composition dictionary
aa_comp = dict(mass.std_aa_comp)

# get Unimod database
db = mass.Unimod()


# Add mod and their composition to the standard a.a dictionary 
for mod, mod_name in modifications.items():
    aa_comp[mod] = db.by_title(mod_name)['composition']

# Add N-term mod to the dictionary [ (mod + H) This how pyteomcis understands this N-term MOD]
aa_comp['pr-'] = aa_comp['pr'] + {'H': 1}






# generate all possible peptidoforms and calculate their m/z's
def get_peptidform(sequences: str, aa_comp: dict, charges: int, ion_type: str) -> pd.DataFrame:
    data: dict = {'sequence': [], 'peptidoform': [], 'z': [], 'm/z': []}
    for sequence in sequences:
        peptidoforms = parser.isoforms(sequence, variable_mods= variable_mods, show_unmodified_termini = True)
        for peptidoform in peptidoforms:
            for z in charges:
                moz = mass.calculate_mass(sequence = peptidoform, aa_comp=aa_comp, charge=z, ion_type= ion_type, max_mods=None)
                data['sequence'].append(sequence)
                data['peptidoform'].append(peptidoform)
                data['z'].append(z)
                data['m/z'].append(moz)
    return pd.DataFrame(data)





df = get_peptidform(sequences, aa_comp, charges= charge, ion_type = ion_type)


df.to_csv(output_file_name, index=False)

print(f'"{output_file_name}" is exported.')
input("Press ENTER to exit")




