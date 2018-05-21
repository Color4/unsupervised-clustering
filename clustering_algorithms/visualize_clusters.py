from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt


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

        isoelectric, aromaticity = calculate_property(cluster, "isoelectric_point"), calculate_property(cluster, "aromaticity")
        plt.scatter(aromaticity, isoelectric, c=color_list(i), label="Cluster " + str(i))

    plt.legend()
    plt.ylabel("Isoelectric point")
    plt.xlabel("Molecular weight")
    plt.show()