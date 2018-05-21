from clustering_algorithms.io import read_active_sites, write_clustering, write_mult_clusterings
from clustering_algorithms.cluster import cluster_by_partitioning, cluster_hierarchically
from clustering_algorithms.visualize_clusters import plot_clusters
import sys


active_sites = read_active_sites(sys.argv[2])

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering, sequences = cluster_by_partitioning(active_sites, 4)
    plot_clusters(sequences)
    write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clustering = cluster_hierarchically(active_sites, 5)
    write_mult_clusterings(sys.argv[3], clustering)




