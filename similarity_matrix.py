import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit.Chem import  AllChem
from rdkit import Chem, DataStructs 
from mhfp.encoder import MHFPEncoder
from rdkit.Chem import rdMHFPFingerprint

import warnings
warnings.filterwarnings("ignore")

class similarity_matrix:
    def __init__(self, data_train, data_test):
        self.data_train = data_train
        self.data_test = data_test
        self.ID = self.data_train.columns[0]
        #self.mol_col = self.data_train.columns[2]
        self.smiles_col = self.data_train.columns[1]
        
        self.ID_test = self.data_test.columns[0]
        #self.mol_col_test = self.data_test.columns[2]
        self.smiles_col_test = self.data_test.columns[1]
    
    
    def tanimoto(self,vector1, vector2):
        a = np.where(vector1 == 1)[0]
        b = np.where(vector2 == 1)[0]
        return len(np.intersect1d(a, b)) / (float(len(a) + len(b)) - len(np.intersect1d(a, b)))
    
    def fp_similarity(self,fp1, fp2):
        return self.tanimoto(fp1, fp2)
    
    def train_process(self):
        self.list_training_name = list(self.data_train[self.ID].values)
        #for trainnames in self.data_train[self.ID]:
            #self.list_training_name.append(trainnames)       
        df_fp = self.data_train.drop([self.ID, self.smiles_col], axis = 1)
        self.list_training_fp = list(df_fp.values)
        
            
    def test_process(self):
        self.list_test_name = list(self.data_test[self.ID_test].values)
    
        df_fp_test = self.data_test.drop([self.ID_test, self.smiles_col_test], axis = 1)
        self.list_test_fp = list(df_fp_test.values)


    
    def fit(self):
        self.train_process()
        self.test_process()
        self.list_data_set=self.list_training_name+self.list_test_name #all data set-> training+test+external
        self.list_data_set_fp=self.list_training_fp+self.list_test_fp #all data set-> training+test+external
        
        
        size=len(self.list_data_set_fp)
        self.matrix=pd.DataFrame()
        for m, i in enumerate(self.list_data_set_fp):
            for n, j in enumerate(self.list_data_set_fp):
                #similarity=DataStructs.TanimotoSimilarity(i,j)
                similarity=self.fp_similarity(i,j)
                self.matrix.loc[self.list_data_set[m],self.list_data_set[n]]=similarity
