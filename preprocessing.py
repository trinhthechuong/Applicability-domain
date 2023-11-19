#import sys
#sys.path.append('ultility')
#from Featurizer import Featurize

import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt



class prepare_dataset:
    
    def __init__(self, data_train, data_test):
        self.data_train = data_train
        self.data_test = data_test
        self.smile_col = self.data_train.columns[1]
        self.mol_col = self.data_train.columns[2]
        #self.activity_col = activity_col
        self.ID = self.data_train.columns[0]
        #self.feature_col = feature_col
        #self.fp_type = fp_type # RDKFp;  ECFPs; MACCs
        # self.feature = Featurize(data = self.data_train, smile_col =self.smile_col,
        #                     activity_col=self.activity_col, ID = self.ID, save_dir = None, m2v_path = None)
        
    def fp_call(self):
        self.train = self.data_train.drop(self.data_train.columns[0:3], axis = 1)
        self.test = self.data_test.drop(self.data_test.columns[0:3], axis = 1)
            
#         if self.feature_col ==  None:
#             if self.fp_type == 'ECFPs':

#                 fp = self.data_train[self.mol_col].apply(self.feature.ECFPs)
#                 X = np.stack(fp.values)
#                 self.train = pd.DataFrame(X)

#                 fp = self.data_test[self.mol_col].apply(self.feature.ECFPs)
#                 X = np.stack(fp.values)
#                 self.test= pd.DataFrame(X)
#         else:
#             self.train = self.data_train[self.feature_col]
#             self.test = self.data_test[self.feature_col]
            
    def pca_reduce(self):
        pca = PCA(n_components=2)
        pca.fit(self.train)

        x_train = pca.transform(self.train)
        x_test = pca.transform(self.test)
        
        df_1 = pd.DataFrame(x_train)
        df_1['Data'] = 'Train'
        df_1["ID"] = self.data_train[self.ID]
        
        df_2 = pd.DataFrame(x_test)
        df_2['Data'] = 'Test'
        df_2["ID"] = self.data_test[self.data_test.columns[0]]
        
        self.df_pca = pd.concat([df_1,df_2], axis = 0).reset_index(drop=True)
        self.df_pca.columns = ['PC1','PC2','Data',"ID"]
        
    def fit(self):
        self.fp_call()
        self.pca_reduce()
