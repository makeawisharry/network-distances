import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random as rnd
import itertools

def avg_deg(G):
    '''
    find average degree of graph G'''
    degrees = sorted((d for n, d in G.degree()), reverse=True)
    return np.mean(degrees)

def fitting(G, plot=False, cumulative=True):
    '''
    find exponent of degree distribution
    plot if desired'''
    n = G.number_of_nodes()
    degrees = sorted((d for n, d in G.degree()), reverse=True)
    x, y = np.unique(degrees, return_counts=True)
    if cumulative:
        y = [sum(y[:i]) for i in range(len(y))]
    else:
        y = [y/n for y in y]
    
    #define model function 
    def exp_func(x,a,b):
        return b*(x**a)

    #actual fitting yields optimized parameters a and b
    paras = np.polyfit(np.log(x), np.log(y), 1)

    if plot:
        x_fit = x
        y_fit = [exp_func(z,paras[0],np.exp(paras[1])) for z in x_fit]
        
        plt.figure(figsize=(7,6))
        plt.scatter(x, y)
        plt.plot(x_fit,y_fit,'r')
        plt.title("Degree distribution",fontsize=18)
        plt.ylabel("$p_k$",fontsize=18)
        plt.xlabel("$k$",fontsize=18)

        ax = plt.gca()

        ax.set_yscale('log')
        ax.set_xscale('log')
        plt.tick_params(axis='both', which='major', colors='k',labelsize=16)
    
    return paras

def generate_network(N, e_k, omega=2.5):
    '''(code from https://github.com/uluturki/hpm)
    create graph with N nodes, e_k avg degree, omega degree exponent'''
    G = nx.empty_graph(N)

    # According to equation 4.29,  omega = (1 + 1/a)
    alpha = 1 / (omega - 1)

    # @TODO Can also solve the series sum analytically, approximation by definite integrals for example.
    # https://math.stackexchange.com/questions/1576502/calculate-finite-p-series

    # Series = SUM 1/n^(alpha) FOR n FROM 1 to N
    harmonic_series = np.power((1 / (np.arange(N) + 1)), alpha)
    c = (e_k * N) / np.sum(harmonic_series)

    # No need to recalculate the list of n, already have the harmonic series, just multiply by c, Eq. 4.28
    n_list = harmonic_series * c

    # n_expected is calculated to be e_k so it is unnecessary to recalculate it, to a very high precision.
    # unless you want to verify that e_k is correct.
    n_expected = e_k

    # Assign each node a hidden parameter n_i, they are generated above
    edges = itertools.combinations(range(N), 2)
    denominator = (n_expected * N)
    # @TODO might be possible to vectorize if we don't mind about the ram, keep in mind.
    for e in edges:
        # compute probability of connecting two nodes based on equation under
        # Figure 4.18
        p = ((n_list[e[0]] * n_list[e[1]]) / denominator)
        if rnd.random() < p:
            G.add_edge(*e)

    return G

def get_lcc(G):
    return G.subgraph(max(nx.connected_components(G), key=len))