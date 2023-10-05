def issvd_diagnostics(res1,res2,true_rows,true_cols):
    
    #res1  - list of length n_views
    #res2  - list of length n_views
    #true_rows - list of arrays - one for each cluster,
    #  with the ith array containing the indices of the ith row cluster
    D = len(res1[1]) #no of data views  # res1[i] considers the ith cluster membership across all views
    k1 = res1[0].shape[1] #number of row clusters
    k2 = res2[0].shape[1] #number of col clusters - again assumes is equal to row - non overlapping
    mats = []
    rel = []
    rev = []
    fscore = []
    fps = []
    fns = []
    fpmats = []
    fnmats = []
    #for each view
    for d in range(D):
        #KxL matrix for false positives matrices for this view
        #The FP is defined as the ratio of number of falsely selected non-zero elements outside of 
        # true bicluster against number of elements in the true bicluster.  
        # the FN is defined as the ratio of number of non-zeros computed by the algorithm in the 
        # true bicluster against the number of elements in the true bicluster.
        fpmat = np.zeros((k1,k2))
        fnmat = np.zeros((k1,k2))
        mat = np.zeros((k1,k2))
        for i in range(k1):
            for j in range(k2):
                A1 = res1[0][:,i].reshape(-1,1)@res1[1][d].T[i,:].reshape(1,-1)
                A = np.where(A1.flatten()>0)[0]
                B1 = res2[0][:,j].reshape(-1,1)@res2[1][d].T[j,:].reshape(1,-1)
                B = np.where(B1.flatten()>0)[0]

                C = set(A).intersection(set(B))
                mat[i,j] = len(C)/(len(A)+len(B)-len(C))
                
                rows = true_rows[j].ravel()
                cols = true_cols[d][j].ravel()
                ele = rows.shape[0]*cols.shape[0]
                fnmat[i,j] = np.sum(A1[np.ix_(rows,cols)]<1) / ele
                A1[np.ix_(rows,cols)] = 0
                fpmat[i,j] = np.sum(A1>0) / ele
        
        fps.append(np.mean(np.min(fpmat,axis=1)))
        fns.append(np.mean(np.min(fnmat,axis=1)))
        mats.append(mat)
        fpmats.append(fpmat)
        fnmats.append(fnmat)
        reld = np.mean(np.max(mat,axis=0))
        revd = np.mean(np.max(mat,axis=1))
        rel.append(reld)
        rev.append(revd)
        fscore.append(2*reld*revd/(reld+revd))
        
    return np.mean(rev), np.mean(rel), np.mean(fscore), np.mean(fps),np.mean(fns)