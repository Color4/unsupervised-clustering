from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import numpy as np

def calculate_property(cluster, property):

    property_list = []

    for sequence in cluster:
        fn = getattr(ProteinAnalysis(sequence), property)
        property_list.append(fn())

    return property_list


def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)


def plot_clusters(sequence_list):
    color_list = get_cmap(len(sequence_list)+1)
    for cluster, i in zip(sequence_list, range(len(sequence_list))):

        isoelectric, aromaticity = calculate_property(cluster, "isoelectric_point"), calculate_property(cluster, "molecular_weight")
        plt.scatter(aromaticity, isoelectric, c=color_list(i), label="Cluster " + str(i))

    plt.legend()
    plt.ylabel("Isoelectric point")
    plt.xlabel("Molecular weight")
    plt.show()


def plot_rankings(sequence_list):

    isoelectric = []
    molecular_weight = []

    for cluster, i in zip(sequence_list, range(len(sequence_list))):
        isoelectric.append(calculate_property(cluster, "isoelectric_point"))
        molecular_weight.append(calculate_property(cluster, "molecular_weight"))

    labels = ["Cluster " + str(i+1) for i in range(len(sequence_list))]

    plt.subplot(1, 2, 1)
    plt.boxplot(isoelectric, labels=labels, patch_artist=True)
    plt.ylabel("Isolectric point")

    plt.subplot(1, 2, 2)
    plt.boxplot(molecular_weight, labels=labels, patch_artist=True)
    plt.ylabel("Molecular Weight")

    plt.show()


def jacard_index(cluster1, cluster2):

    union = list(set(cluster1).union(cluster2))
    intersection = list(set(cluster1) & set(cluster2))

    return len(intersection)/len(union)


def assess_similarity(cluster_list1, cluster_list2):
    cluster_size = len(cluster_list1)

    similarity = np.zeros((cluster_size, cluster_size))

    for i in range(cluster_size):
        for j in range(cluster_size):
            similarity[i, j] = jacard_index(cluster_list1[i], cluster_list2[j])

    return similarity



