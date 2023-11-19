
import numpy as np
#import modin.pandas as pd
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import Isomap
from matplotlib.pyplot import figure
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from sklearn import manifold
import os


class similarity_convexhull:
    def __init__(self, data, list_training_fp, list_test_fp, figsize, save_fig=False):
        self.data = data
        self.list_training_fp = list_training_fp
        self.list_test_fp = list_test_fp
        self.figsize = figsize
        self.save_fig = save_fig
        
    def point_in_hull(self, point, hull, tolerance=1e-12):
        return all(
            (np.dot(eq[:-1], point) + eq[-1] <= tolerance)
            for eq in hull.equations)
    
    def isomap(self):
        mds = Isomap(n_components=2, n_jobs= -1)
        results = mds.fit(self.data)
        self.coords = results.embedding_
    
    def mds(self):
        mds = manifold.MDS(n_components=2, dissimilarity="euclidean", random_state=42)
        results = mds.fit(self.data)
        self.coords = results.embedding_
        
    def visualize(self, coords):
        coords= (coords - np.min(coords)) / (np.max(coords) - np.min(coords))#normalize
        coords_training=coords[:len(self.list_training_fp)]
        coords_test=coords[len(self.list_training_fp):len(self.list_test_fp)+len(self.list_training_fp)]
        
        hull = ConvexHull(coords_training)
        
        in_out = []
        for p in coords_test:
            point_is_in_hull = self.point_in_hull(p, hull)
            in_out.append(point_is_in_hull)

        self.df_test = pd.DataFrame(coords_test)
        self.df_test['convex'] = in_out
        




        coords_test_in = coords_test[in_out]
        coords_test_out = self.df_test[self.df_test['convex']==False].drop(['convex'], axis =1)
    
        
        idx = coords_test_out.index
        out = self.data.iloc[len(self.list_training_fp):len(self.list_test_fp)+len(self.list_training_fp)]
        out = out.iloc[idx,:]
        display(out.head())
        
        
        
        coords_test_out = coords_test_out.values
        plt.scatter(coords_training[:, 0], coords_training[:, 1], marker = 'o',label='Training')#training
     

        '''
        Visualize the convex hull
        '''
        #hull = ConvexHull(coords_training)

        for simplex in hull.simplices:
            plt.plot(coords_training[simplex, 0], coords_training[simplex, 1], 'k-')



        plt.scatter(coords_test_in[:, 0], coords_test_in[:, 1], marker = '^',label='Test_in', color ='g')#test
        plt.scatter(coords_test_out[:, 0], coords_test_out[:, 1], marker = 'd',label='Test_out', color = 'r')#test
        # for label, x, y in zip(list_data_set[len(list_training_fp):len(list_test_fp)+len(list_training_fp)], coords_test[:, 0], coords_test[:, 1]): #show molecule name 
        #     plt.annotate(label,xy = (x, y), xytext = (0, 0),textcoords = 'offset pixels', ha = 'center', va = 'bottom', fontsize=8) #show molecule name 

        #plt.scatter(coords_external[:, 0], coords_external[:, 1], marker = 'X', label='External')#external
        # for label, x, y in zip(list_data_set[len(list_test_fp)+len(list_training_fp):], coords_external[:, 0], coords_external[:, 1]): #show molecule name 
        #     plt.annotate(label,xy = (x, y), xytext = (0, 0),textcoords = 'offset pixels', ha = 'center', va = 'bottom', fontsize=8)#show molecule name 
        plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",mode="expand", borderaxespad=0, ncol=3,shadow=True, fontsize='12')
        plt.xlabel("isomap1",fontweight='bold', fontsize = 16)
        plt.ylabel("isomap2",fontweight='bold', fontsize = 16)
        self.index_out = out.index
        return self.index_out
    
    def fit(self):
        plt.figure(figsize=self.figsize, dpi=600)
        plt.subplot(121)
        self.isomap()
        self.visualize(coords =self.coords)
        self.index_out_isomap = self.index_out
        plt.subplot(122)
        self.mds()
        self.visualize(coords =self.coords)
        self.index_out_mds = self.index_out
        if self.save_fig == True:
            self.folder = 'AD_img'
            isExist = os.path.exists(self.folder)
            if not isExist:
                os.makedirs(self.folder)
            plt.savefig('AD_img/convexhull_similarity.png', transparent = False)
