import numpy as np
import math

def patient_correlation_measure(F, M):
    
    F = F.T
    num_genes = M.shape[1]
    dist = np.zeros((num_genes, num_genes))
        
    for i in range(num_genes):
        for j in range(i+1, num_genes):
            
            dist[i,j] = M[i,j] + math.sqrt(F[i] + F[j])/10 #smaller values are shrunk and larger values are inflated 
    
    dist = dist + dist.T + np.eye(num_genes)
    
    return dist
