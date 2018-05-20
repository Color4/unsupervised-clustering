import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically
import sys

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
    sys.exit(0)

active_sites = read_active_sites(sys.argv[2])

# if sys.argv[4] == "-T":
# 	active_sites = read_active_sites("test/data/")

# 	print("Clustering using Partitioning method")
# 	for i in range(10): 
# 		clustering = cluster_by_partitioning(active_sites, 3)

# 		correct = 0 
# 		for cluster in clustering: 
# 			if "moleculeA1" and "moleculeA" in cluster: 
# 				correct +=1 
# 			if "moleculeE1" and "moleculeE" in cluster: 
# 				correct +=1 
# 	print(correct/10)
	


# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster_by_partitioning(active_sites, 5)
    write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clusterings = cluster_hierarchically(active_sites, 5)
    write_mult_clusterings(sys.argv[3], clusterings)




