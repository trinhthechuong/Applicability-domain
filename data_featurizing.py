#!/usr/bin/env python
# coding: utf-8

# In[14]:


import math
import numpy as np
import pandas as pd
import os
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import seaborn as sns

from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.rdBase import BlockLogs
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem import PandasTools, AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors


from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Avalon import pyAvalonTools as fpAvalon
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D,Generate
#from mol2vec.features import mol2alt_sentence, MolSentence, DfVec, sentences2vec
#from mol2vec.helpers import depict_identifier, plot_2D_vectors, IdentifierTable, mol_to_svg
#from gensim.models import word2vec
#from mol2vec import features
#from mol2vec import helpers
from mhfp.encoder import MHFPEncoder
#from mordred import Calculator, descriptors
from tqdm import tqdm # progress bar
tqdm.pandas()

from sklearn.preprocessing import MinMaxScaler
from pathlib import Path
import math
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, PandasTools
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.Draw import IPythonConsole
from tqdm.auto import tqdm
tqdm.pandas()

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit import DataStructs, Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.rdBase import BlockLogs
import warnings
warnings.filterwarnings(action='ignore')


class preprocess():
    def __init__(self, data, ro5):
        self.data = data
        self.ro5 = ro5
        #self.save_dir = save_dir
        self.ID = self.data.columns[0]
        
        
    def calculate_ro5_properties(self,smiles, fullfill = 4):
    # RDKit molecule from SMILES
        molecule = Chem.MolFromSmiles(smiles)
    # Calculate Ro5-relevant chemical properties
        molecular_weight = Descriptors.ExactMolWt(molecule)
        n_hba = Descriptors.NumHAcceptors(molecule)
        n_hbd = Descriptors.NumHDonors(molecule)
        logp = Descriptors.MolLogP(molecule)
    #tpsa = Descriptors.TPSA(molecule)
    # Check if Ro5 conditions fulfilled
        conditions = [molecular_weight <= 500, n_hba <= 10, n_hbd <= 5, logp <= 5]
        ro5_fulfilled = sum(conditions) >= fullfill
        return ro5_fulfilled

    def standardize(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        # removeHs, disconnect metal atoms, normalize the molecule, reionize the molecule
        clean_mol = rdMolStandardize.Cleanup(mol) 
        # if many fragments, get the "parent" (the actual mol we are interested in) 
        parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol)
        # try to neutralize molecule
        uncharger = rdMolStandardize.Uncharger() # annoying, but necessary as no convenience method exists
        uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)
        te = rdMolStandardize.TautomerEnumerator() # idem
        taut_uncharged_parent_clean_mol = te.Canonicalize(uncharged_parent_clean_mol)
        return taut_uncharged_parent_clean_mol 
        
        
    def filter_data(self):
        self.data['Canomicalsmiles'] = self.data.iloc[:,1].apply(Chem.CanonSmiles)
        self.index_0 = self.data.index.tolist()
        
        self.data_ro5 = self.data[self.data['Canomicalsmiles'].progress_apply(self.calculate_ro5_properties, fullfill = self.ro5)]
        self.index_ro5 = self.data_ro5.index.tolist()
        self.index_vi_ro5 = [x for x in self.index_0 if x not in self.index_ro5]
        #self.data_vi_ro5 = self.data.iloc[self.index_vi_ro5,:]
        
        #self.data = self.data_ro5.copy()
        self.data.reset_index(drop=True, inplace = True)
        block = BlockLogs()
        self.data['Molecule'] = self.data['Canomicalsmiles'].progress_apply(self.standardize)
        
        return self.data
    
    
    def featurizing(self):
        print("STANDARDIZING MOLECULES...")
        block = BlockLogs()
            #std = standardization(data=self.data,ID=self.ID, smiles_col=self.smile_col, active_col=self.activity_col, ro5 = self.ro5)
        self.data = self.filter_data()
        del block
            

        
        #. secfp
        self.secfp= self.data.copy()
        print("CALCULATING SECFP FINGERPRINTS...")
        self.secfp["secfp"] = self.secfp["Canomicalsmiles"].progress_apply(MHFPEncoder.secfp_from_smiles)
        X = np.stack(self.secfp["secfp"].values)
        d = pd.DataFrame(X)
        #self.secfp= pd.concat([self.secfp, d], axis = 1)
        self.secfp_ad= pd.concat([self.secfp, d], axis = 1).drop([self.secfp.columns[1],"secfp"], axis =1)
        self.secfp_visualize= pd.concat([self.secfp, d], axis = 1).drop(["Molecule",self.secfp.columns[1],"secfp"], axis =1)
        #self.secfp.columns[1]
        #self.secfp_visualize.to_csv(f"{self.save_dir}Secfp.csv", index= False)         
        print("FINISH CALCULATING!")
        #return self.secfp_visualize
        
    def fit(self):
        #self.folder = 'featurized_data'
        #isExist = os.path.exists(self.folder)
        #if not isExist:
            #os.makedirs(self.folder)
        #self.save_dir = self.save_dir+'/'+self.folder+"/"

        self.featurizing()
        #os.chdir(os.getcwd())
        #return self.secfp_visualize
    

