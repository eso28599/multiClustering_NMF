import pandas as pd
import numpy as np
from iSSVD.functions import issvd, issvd_diagnostics, gen_sim_vec, gen_tmp
from numpy.random import seed

data_views = pandas.read_excel('/Users/ellaorme/GitHub/Data_integration/data_views.xlsx', sheet_name = ["Sheet 1","Sheet 2","Sheet 3"])
test1_data = [np.array(data) for sheet,data in data_views.items()]
test1 = issvd(test1_data, standr=False,pointwise=True,steps=100,size=0.5,
            vthr = 0.9,ssthr=[0.6,0.65],nbicluster=10,rows_nc=False,cols_nc=False,col_overlap=False
            ,row_overlap=False,pceru=0.9,pcerv=0.9,merr=0.0001,iters=100)
test1
#assess answers in the same way 
# what are we interested in 
# Bicluster samples identified by iSSVD
Rows = test1['Sample_index']
# Bicluster variables identified by iSSVD
Cols = test1['Variable_index']


import seaborn as sns
import matplotlib.pyplot as plt

df = test1_data
new_df = []
for d in range(2):
    cs = []
    d=0
    col = np.array([],dtype=int)
    for i in range(3):
        r1 = df[d][Rows[i],:]
        col = np.append(col,Cols[i][d])
        c1 = r1[:,col]
        c2 = np.delete(r1, col, axis=1)
        c3 = np.concatenate([c1,c2], axis=1)
        cs.append(c3)
    new_df.append(np.vstack(cs))
    
fign, axesn = plt.subplots(nrows=1, ncols=2, figsize=(15,3))
sns.heatmap(new_df[0],ax=axesn[0])
sns.heatmap(new_df[0],ax=axesn[1])
plt.show()
# True bicluster samples and variables - which indices
#each dataset with 200 samples, same clusters across views
#rowClusters1<-list(c(75, 75, 50),c(75, 75, 50),c(75, 75, 50))
#100, 50,250 features respectively
#colClusters1<-list(c(30,30,40),c(10, 20, 20),c(100,50,100))
row_ind = [[np.arange(75),np.arange(75),np.arange(75)],[np.arange(75,150),
             np.arange(75,150),np.arange(75,150)],[np.arange(150,200),
             np.arange(150,200),np.arange(150,200)]]
row_ind = [np.arange(75),np.arange(75,150),np.arange(150,200)]
col_ind = [[np.arange(30),np.arange(10),np.arange(100)],[np.arange(30,60),
             np.arange(10,30),np.arange(100,150)],[np.arange(60,100),
             np.arange(30,50),np.arange(150,250)]]



res1tmp, res2tmp = gen_tmp(Rows[0],Cols[0], row_ind, col_ind,n=200,p=250,D=3)

rev, rel, f, fp, fn = issvd_diagnostics(res1tmp,res2tmp,row_ind[0],col_ind[0])

#their guide 
seed(25)
data, rows, cols = gen_sim_vec(n=200,p=1000,D=3,rowsize=50, colsize=100, 
                               numbers=1,sigma=0.1,nbicluster=4, orthonm=False)


df = data[0]
test2 = issvd(X=df,standr=False,pointwise=True,steps=100,size=0.5,
            vthr = 0.9,ssthr=[0.6,0.65],nbicluster=10,rows_nc=False,cols_nc=False,col_overlap=False
            ,row_overlap=False,pceru=0.1,pcerv=0.1,merr=0.0001,iters=100)

Rows2 = test2['Sample_index']

# Bicluster variables identified by iSSVD
Cols2 = test2['Variable_index']
res1tmp2, res2tmp2 = gen_tmp(Rows2,Cols2, rows[0], cols[0],n=200,p=1000,D=2)


rev, rel, f, fp, fn = issvd_diagnostics(res1tmp2,res2tmp2,rows[0], cols[0])

def gen_tmp(Rows, Cols, true_Rows, true_Cols,n,p,D):
    
    # Rows
    row1 = []
    row2 = []
    
    for k in range(len(true_Rows)):
        tmp2 = np.zeros(n)
        tmp2[true_Rows[k]] = 1
        row2.append(tmp2.reshape(-1,1))
    rowf = np.concatenate(row2,axis=1)
    
    for k in range(len(Rows)):
        tmp1 = np.zeros(n)
        tmp1[Rows[k]] = 1
        row1.append(tmp1.reshape(-1,1))
    rowg = np.concatenate(row1, axis=1)
    
    # Cols
    colg = []
    colf = []
    for d in range(D):
        col1 = []
        col2 = []
        
        for k in range(len(Cols)):
            tmp1 = np.zeros(p)
            tmp1[Cols[k][d]] = 1
            col1.append(tmp1.reshape(-1,1))
        colg.append(np.concatenate(col1,axis=1))
        
        for k in range(len(true_Rows)):
            tmp2 = np.zeros(p)
            tmp2[true_Cols[d][k]] = 1
            col2.append(tmp2.reshape(-1,1))
        colf.append(np.concatenate(col2, axis=1))
    
    return [rowg,colg], [rowf,colf]