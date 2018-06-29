# my imports
from pyxdameraulevenshtein import damerau_levenshtein_distance as distance
import random 
import numpy as np
import itertools


def compute_similarity(a, b):
    """
    Compute the similarity between two given ActiveSite instances.

    Distance is calculated using the Damerau-Levenshtein distance. This algorithm quantifies the
    'distance' between two strings by counting the number of operations that are needed to transform one string
    into another. In this algorithm the operations allowed are inserting a letter, deleting a letter, subsituting a letter
    or transposing a letter. This is an improvement on the Levenshtein distance algorithm that preceeded it as
    transpositions are common errors, both in the written language, and may have a viable biological relevance.

    Implementation from-  https://github.com/gfairchild/pyxDamerauLevenshtein
    Algorithmic understanding from- https://www.mathworks.com/matlabcentral/cody/problems/2309-calculate-the-damerau-levenshtein-distance-between-two-strings

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    return distance(a, b)


def avg(l):
    """
    Computes the average of a list

    :param l: a list
    :return: average of list
    """
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

    """
    takes an amino acid abbreviation and converts it to its single letter amino acid code

    codon table taken from- http://130.88.97.239/bioactivity/aacodefrm.html
    :param sequence: a sequence of amino acid codes
    :return: sequence of one letter amino acids
    """

    aa = { 'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'ASX':'B', 'LEU':'L',
        'CYS':'C', 'GLU':'E','GLN':'Q', 'GLX':'Z', 'GLY':'G', 'HIS':'H', 
        'ILE':'I', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 
        'THR':'T', 'TYR':'Y', 'VAL':'V', 'TRP':'W'} 

    aa_seq = ""

    for item in sequence: 
        aa_seq += aa[item.type]

    return aa_seq

def get_aa_seq(active_sites):
    """

    :param active_sites:
    :return:
    """

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

    cents = generate_cents(seq, k)
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

    return convert_to_object_2d(clusters, seq_dict), clusters


def generate_cents(seq, k):

    cents = random.sample(seq, k)
    dups = True

    while dups:
        if not any(cents.count(x) > 1 for x in cents):
            dups = False
        else:
            cents = random.sample(seq, k)

    return cents


def convert_to_object_2d(clusters, seq_dict):
    clusters = [[seq_dict[item] for item in cluster] for cluster in clusters]
    return clusters


def convert_to_object_3d(clusters_list, seq_dict):

    all_clusters = []

    for clusters in clusters_list:
        new_cluster = []
        for cluster in clusters:
            new_cluster.append(seq_dict[cluster])
        all_clusters.append(new_cluster)

    return all_clusters


def single_linkage(cluster1, cluster2):

    avg_cluster1 = avg_string(cluster1)
    avg_cluster2 = avg_string(cluster2)

    return distance(avg_cluster1, avg_cluster2)


def cluster_hierarchically(active_sites, k):

    seq, seq_dict = get_aa_seq(active_sites)
    clusters = [[seq[i]] for i in range(len(seq))]

    while len(clusters) > k:

        i, j = find_most_similar(clusters)
        clusters = merge_clusters(clusters, i, j)

    return convert_to_object_3d(clusters, seq_dict), clusters


def find_most_similar(clusters):

    size = len(clusters)
    similarities = {}

    for i in range(size-1):
        for j in range(i+1, size-1):
            similarities[(i, j)] = single_linkage(clusters[i], clusters[j])

    return min(similarities, key=similarities.get)


def merge_clusters(clusters, i, j):

    new_clusters = []

    for item in clusters:
        if item not in [clusters[i], clusters[j]]:
            new_clusters.append(item)

    new_clusters.append(clusters[i] + clusters[j])

    return new_clusters
