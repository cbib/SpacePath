import sys
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from scipy.sparse import csr_matrix
import networkx as nx
import time
from tqdm import tqdm

def time_formatting(timeseconds):  # second copy, first in experiments_pws_scores
    hours = int(timeseconds / 3600)
    remaining_seconds = timeseconds % 3600
    minutes = int(remaining_seconds / 60)
    seconds = np.around(remaining_seconds % 60, 1)
    return hours, minutes, seconds


def get_cosine_simil(b, y):
    # print("b..." , b.reshape(1, -1))
    # print("y...", y.reshape(1, -1))
    tmp = cosine_similarity( y.reshape(1, -1), b.reshape(1, -1))
    return tmp

def print_dots(iteration, interval=5):
    if iteration > 0 and iteration % interval == 0:
        sys.stdout.write('../bank4zonist')
        sys.stdout.flush()

def compute_upright_triangle(X_arr):
    print("\nComputing upright triangle (this is slow but memory low-cost)")
    o_dict = {}
    eval_sample = 0
    while eval_sample < (X_arr.shape[0] -1):
        remaining = eval_sample + 1
        # Y = np.reshape(X_arr[eval_sample, :],
        #                (X_arr[eval_sample, :].shape[0], 1 ))
        Y = X_arr[eval_sample, :]
        X = X_arr[remaining:, :]
        tmp = np.apply_along_axis(get_cosine_simil, 1, arr=X, y=Y)
        o_dict[eval_sample] = np.array(tmp).flatten()
        print_dots(eval_sample)
        eval_sample += 1

    return o_dict


def compute_graph(o_dict, cutoff_weight):
    print("\nComputing graph")
    undirected_graph = nx.Graph()
    for u in o_dict.keys():
        for count in range(len(o_dict[u])):
            v = count + u + 1  # adjust v index
            if o_dict[u][count] >= cutoff_weight:
                # if u > v:
                #     undirected_graph.add_edge(v, u, weight=o_dict[u][count])
                # else:
                #     undirected_graph.add_edge(u, v, weight=o_dict[u][count])
                undirected_graph.add_edge(u, v, weight=o_dict[u][count])
    return undirected_graph


# start

cutoff_cos_neighbors = 0.7

X_arr = np.random.rand(80000, 50)  # 30000
X_arr[X_arr < 0.2 ] = 0

# X_arr = np.array([[4,6,8, 0], [1, 0.4, 1.2, 1.5], [-1, -0.4, -1.2, -1.5],
#                   [7, 5, 7, 0], [6, 7, 4, 6]])
# print(X_arr)

print("input dimensions: ", X_arr.shape)
t0 = time.time()
print("computing cosine similarities with my method")
o_dict = compute_upright_triangle(X_arr)
Goool = compute_graph(o_dict, cutoff_cos_neighbors)
print("elapsed time: ", time_formatting(time.time() - t0))
#else:
print("computing cosine similarities with numpy arrays")
t0 = time.time()
#print("this will take long, using sparse matrix")
#X_arr = csr_matrix(X_arr)  # time did not improved with sparse array
computed_matrix = cosine_similarity(X_arr)
print(computed_matrix.shape)
np.fill_diagonal(computed_matrix, 0)  # remove self loops
computed_matrix[computed_matrix < cutoff_cos_neighbors] = 0
G = nx.from_scipy_sparse_array(csr_matrix(computed_matrix))
del computed_matrix
print("elapsed time: ", time_formatting(time.time() - t0))

# for edge in G.edges():
#     print(edge[0], edge[1], G.get_edge_data(edge[0], edge[1]))   #  {'weight': 0.7363411622377004}
# print("!!!!!!")
# for edge in Goool.edges():
#     print(edge[0], edge[1], Goool.get_edge_data(edge[0], edge[1]))   #  {'weight': 0.7363411622377004}



print("\nComparing results: are final Graphs equal?")
import networkx.algorithms.isomorphism as iso
em = iso.numerical_edge_match("weight", 1)

print(nx.is_isomorphic(G, Goool, edge_match=em))





#
# similarity_matrix = computed_matrix.copy()
# similarity_matrix[similarity_matrix < cutoff_cos_neighbors] = 0
#
# print(similarity_matrix.shape)



# opposites_matrix = computed_matrix.copy()
# opposites_matrix[opposites_matrix > -0.7] = 0
# print(opposites_matrix.shape)
