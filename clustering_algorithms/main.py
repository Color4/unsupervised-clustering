from clustering_algorithms.io import read_active_sites, write_clustering
from clustering_algorithms.cluster import cluster_by_partitioning, cluster_hierarchically
from clustering_algorithms.assess_clusters import plot_clusters, plot_rankings, jacard_index
import sys

active_sites = read_active_sites(sys.argv[2])

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering, sequences = cluster_by_partitioning(active_sites, 3)
    plot_clusters(sequences)
    write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clustering, sequences = cluster_hierarchically(active_sites, 3)
    plot_rankings(sequences)
    write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-T':

    _, partitioning = cluster_by_partitioning(active_sites, 3)
    _, hierarchical = cluster_hierarchically(active_sites, 3)




