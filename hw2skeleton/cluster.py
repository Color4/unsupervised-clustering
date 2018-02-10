from .utils import Atom, Residue, ActiveSite
from .io import read_active_sites, read_active_site, write_clustering, write_mult_clusterings

# my imports 
from pyxdameraulevenshtein import damerau_levenshtein_distance as distance
import random 
import numpy as np
import scipy 
import itertools





def compute_similarity(a, b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """

    return distance(a, b)

def avg(l):
    return sum(l)/len(l)

def avg_string(strings):
    """
    Returns the "average" string in the list
    AKA the string with the least average distance to all the other strings in the list
    :param strings: List of strings
    :return: The average string as defined above
    """

    answer = min([(y, avg([compute_similarity(y, x) for x in strings])) for y in strings], key=lambda x: x[1])[0]
    return answer

def convert_to_aa_str(sequence):

    aa = { 'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'ASX':'B', 'LEU':'L',
        'CYS':'C', 'GLU':'E','GLN':'Q', 'GLX':'Z', 'GLY':'G', 'HIS':'H', 
        'ILE':'I', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 
        'THR':'T', 'TYR':'Y', 'VAL':'V', 'TRP':'W'} 

    aa_seq = ""

    for item in sequence: 
        aa_seq+= aa[item.type]

    return aa_seq


def get_aa_seq(active_sites): 

    seq_dict = {}
    for site in active_sites: 
        seq_dict[convert_to_aa_str(site.residues)] = site

    return [convert_to_aa_str(site.residues) for site in active_sites], seq_dict


def cluster_by_partitioning(active_sites, k):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """

    seq, seq_dict = get_aa_seq(active_sites)

    cents = random.sample(seq, k)

    clusters = [[] for i in range(k)]
    change = True

    while change:

        change = False
        for s in seq:

            distances = [(idx, compute_similarity(s, cent)) for idx, cent in enumerate(cents)]
            min_idx, _ = min(distances, key=lambda t:t[1])
            cluster = clusters[min_idx]
            if s not in cluster:
                cluster.append(s)
                change = True
                for idx, other_cluster in enumerate(clusters):
                    if idx != min_idx and s in other_cluster:
                        other_cluster.remove(s)

        cents = [avg_string(cluster) for cluster in clusters]
    
    return convert_to_object_2d(clusters, seq_dict)


def convert_to_object_2d(clusters, seq_dict):
    clusters = [[seq_dict[item] for item in cluster] for cluster in clusters]
    return clusters

def convert_to_object_3d(clusters_list, seq_dict):
    clusters_list = [[[seq_dict[item] for item in cluster] for cluster in clusters] for clusters in clusters_list]
    return clusters_list



def cluster(clusters):

    size = len(clusters)
    sim_matrix = np.empty((size, size))

    for i in range(size):
        for j in range(size):
            sim_matrix[i][j] = single_linkage(clusters[i], clusters[j])
    i, j, distance = minimum_position(sim_matrix)
    return merge_clusters(i, j, clusters, size), distance

def single_linkage(cluster1, cluster2):

    tuples = list(itertools.product(cluster1, cluster2))
    distances = [distance(t[0], t[1]) for t in tuples]

    return min(distances)


def minimum_position(matrix):

    np.fill_diagonal(matrix, np.inf)
    min_loc = np.where(matrix == np.min(matrix))

    min_i = (min_loc[0])[0]
    min_j = (min_loc[1])[0]


    return min_i, min_j, matrix.min()


def merge_clusters(i, j, clusters, size):

    merged_cluster = []
    old_clusters = []
    new_clusters_all = []

    for cluster_count in range(size):
        if cluster_count == i or cluster_count == j:
            merged_cluster += clusters[cluster_count]
        else:
            old_clusters.append(clusters[cluster_count])
    return([merged_cluster] + old_clusters)
   

def cluster_hierarchically(active_sites, k):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                 

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    seq, seq_dict = get_aa_seq(active_sites)
    all_clusters = [[]]

    clusters = [[seq[i]] for i in range(len(seq))]
    distances = []

    max_d = 0
    max_d_loc = 0

    while len(clusters) > 1:
        clusters, distance = cluster(clusters)
        distances.append(distance)
        all_clusters += [clusters]
        if len(all_clusters)==k: 
            break 
       
    return(convert_to_object_3d(all_clusters, seq_dict))








 
