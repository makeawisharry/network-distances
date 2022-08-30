import networkx as nx
import numpy as np
import time
import random

import linalg_tools as lt

def expected_rw_distance(G, node1, node2, t=False):
    '''
    Nodes should be named 0, 1, 2, ..., n-1
    
    Compare https://hal.inria.fr/hal-01944246/document, equations (10) and (17)
    '''
    n = G.number_of_nodes()
    
    # Net flow vector of G
    S = np.zeros(n)
    S[node1] = 1
    S[node2] = -1
    
    if t:
        start1 = time.time()
        print('Calculating Laplacian of G')
    
    # Laplacian matrix of G
    A = nx.adjacency_matrix(G).todense()
    D = np.diag([val for (node, val) in G.degree()])
    L = D-A # As it is now, L is singular -> use method in (Newman, "A measure of betweenness centrality based on random walks")
    
    if t:
        end1 = time.time()
        print(end1-start1)
        start2 = time.time()
        print('Calculating potential vector')
    
    v = random.sample([ele for ele in list(range(n)) if ele not in {node1, node2}], 1)[0] # choose random index
    L = np.delete(L, v, axis=0) # delete v-th row of L 
    L = np.delete(L, v, axis=1) # delete v-th column of L 
    S = np.delete(S, v) # delete v-th entry of S 
    
    # If matrix is positive definite, solve system iteratively to calculate the potential vector - else, solve directly
    if lt.is_pd(L):
        #print('Matrix is PD, solving iteratively')
        V = lt.gauss_seidel(L, S, tolerance=1e-10, max_iterations=1000)
    else:
        #print('Matrix is not PD, solving directly')
        V = np.linalg.solve(L, S) 
        
    V = np.insert(V, v, 0) # insert v-th entry of V as 0
    
    if t:
        end2 = time.time()
        print(end2-start2)
        start3 = time.time()
        print('Calculating distance')
    
    # Calculate distance with equation (10) from link above
    l = [abs(V[edge[0]]-V[edge[1]]) for edge in G.edges]
    distance = sum(l)
    
    if t: 
        end3 = time.time()
        print(end3-start3)
        print('Done')

    return distance


def dir_dissimilarity(i, j, theta, G, large_number):
    '''
    Algorithm 1 from "A Family of Dissimilarity Measures between Nodes Generalizing both the Shortest-Path and the Commute-time Distances" by Luh Yen, Marco Saerens et al.
    '''
    if nx.has_path(G, i, j)==False:
        return np.inf
    
    else:
        # Find C and P_ref
        n = G.number_of_nodes()
        C = np.matrix(nx.adjacency_matrix(G).todense())
        P_ref = C/np.tile(np.sum(C, 1), (1, C.shape[1]))
        C[C==0] = large_number

        # Step 1
        C_tilde = C.copy()
        C_tilde[j] = np.full((1, n), large_number)

        # Step 2
        W_tilde = np.multiply(P_ref, np.exp(-theta*C_tilde))

        # Steps 3-5
        #assert np.max(np.absolute(np.linalg.eigvals(W_tilde)))<1, 'Spectral radius is greater than one :('

        # Step 6
        L = np.eye(n, n)-W_tilde
        e_i = lt.standard_unit_vector(n, i)
        e_j = lt.standard_unit_vector(n, j)
        '''if lt.is_pd(L):
            print('Matrix is PD, solving iteratively')
            z_i = lt.gauss_seidel(L, e_i, tolerance=1e-10, max_iterations=500).reshape((n, 1))
            z_j = lt.gauss_seidel(L, e_j, tolerance=1e-10, max_iterations=500).reshape((n, 1))
        else:
            print('Matrix is not PD, solving directly')
            z_i = np.linalg.solve(L, e_i).reshape((n, 1))
            z_j = np.linalg.solve(L, e_j).reshape((n, 1))'''

        z_i = np.linalg.solve(np.transpose(L), e_i).reshape((n, 1))
        z_j = np.linalg.solve(L, e_j).reshape((n, 1))

        # Step 7
        A = np.multiply(C_tilde, W_tilde)
        diss = np.dot(np.transpose(z_i), np.dot(A, z_j))/z_i[j]
        return diss[0, 0]


def diss_distance(i, j, theta, G, large_number):
    if nx.has_path(G, i, j)==False:
        return np.inf
    else:
        return (dir_dissimilarity(i, j, theta, G, large_number)+dir_dissimilarity(j, i, theta, G, large_number))/2