import sys
from clustering_algorithms.io import read_active_sites, write_clustering, write_mult_clusterings
from clustering_algorithms.cluster import cluster_by_partitioning, cluster_hierarchically
import sys


active_sites = read_active_sites(sys.argv[2])

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster_by_partitioning(active_sites, 5)
    write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clusterings = cluster_hierarchically(active_sites, 5)
    print(clusterings)
    write_mult_clusterings(sys.argv[3], clusterings)




