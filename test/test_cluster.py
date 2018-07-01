from clustering_algorithms import cluster
from clustering_algorithms import io
import os
import unittest


class TestSimilarityClustering(unittest.TestCase):

    def test_no_similarity(self):
        filename_a = os.path.join("data", "276.pdb")
        filename_b = os.path.join("data", "4629.pdb")

        activesite_a = io.read_active_site(filename_a)
        activesite_b = io.read_active_site(filename_b)

        seq_a = cluster.convert_to_aa_str(activesite_a.residues)
        seq_b = cluster.convert_to_aa_str(activesite_b.residues)

        # should be 8 because there are 8 edit operations between the strings
        self.assertEqual(cluster.compute_similarity(seq_a, seq_b), 8)

    def test_identical_similarity(self):
        filename_a = os.path.join("data", "276.pdb")
        filename_b = os.path.join("data", "276.pdb")

        activesite_a = io.read_active_site(filename_a)
        activesite_b = io.read_active_site(filename_b)

        seq_a = cluster.convert_to_aa_str(activesite_a.residues)
        seq_b = cluster.convert_to_aa_str(activesite_b.residues)

        # should be 0 because the strings are identical
        self.assertEqual(cluster.compute_similarity(seq_a, seq_b), 0)


class ClusterTestCase(unittest.TestCase):

    def assertClustersEqual(self, test_clusters, expected_clusters):
        """
        Asserts that the test clusters are equal to the expected clusters, with possible variance in order both
        within the cluster and the order of the clusters together.
        :param test_clusters: List of clusters, each filled with activesite objects.
        :param expected_clusters: Expected results of clustering, list of list of strings.
        """
        clusters_b_set = frozenset([frozenset(x) for x in expected_clusters])
        clusters_a_set = frozenset([frozenset(x) for x in self.cluster_names(test_clusters)])
        try:
            self.assertSetEqual(clusters_a_set, clusters_b_set)
        except AssertionError:
            raise AssertionError("{} is not the same clustering as {}".format(test_clusters, expected_clusters))

    def cluster_names(self, clusters):
        """
        Returns the name of each cluster for testing.
        :param clusters: List of clusters, each filled with activesite objects.
        :return: List of clusters, each filled with activesite names for testing.
        """
        return [[site.name for site in cluster] for cluster in clusters]


class TestSimilarityClustering(ClusterTestCase):

    def test_partition_clustering_individual(self):
        # tractable subset
        pdb_ids = [276, 4629, 10701]

        active_sites = []
        for id in pdb_ids:
            filepath = os.path.join("data", "%i.pdb"%id)
            active_sites.append(io.read_active_site(filepath))

        clusters = cluster.cluster_by_partitioning(active_sites, len(pdb_ids))[0]
        # Should each be in their own cluster.
        self.assertClustersEqual(clusters, [["276"], ["10701"], ["4629"]])

    def test_partition_clustering_one(self):
        # tractable subset
        pdb_ids = [276, 4629, 10701]

        active_sites = []
        for id in pdb_ids:
            filepath = os.path.join("data", "%i.pdb" % id)
            active_sites.append(io.read_active_site(filepath))

        clusters = cluster.cluster_by_partitioning(active_sites, 1)[0]
        # Should all be in one cluster
        self.assertClustersEqual(clusters, [["276", "10701", "4629"]])

    def test_partition_clustering_two(self):
        # tractable subset
        pdb_ids = [37438, 55996, 81859]
        # Try it four times, since test is flaky (because of inherent randomness)
        # Normally would use flaky.
        for i in range(4):
            try:
                active_sites = []
                for id in pdb_ids:
                    filepath = os.path.join("data", "%i.pdb" % id)
                    active_sites.append(io.read_active_site(filepath))

                clusters = cluster.cluster_by_partitioning(active_sites, 2)[0]

                # 37438 should be clustered with 81859
                self.assertClustersEqual(clusters, [["37438", "55996"], ["81859"]])
                break
            except AssertionError:
                pass
        else:
            raise AssertionError("Partitioning failed on all 4 attempts")

class TestHierarchicalClustering(ClusterTestCase):


    def test_hierarchical_clustering_individual(self):
        # tractable subset
        pdb_ids = [276, 4629, 10701]

        active_sites = []
        for id in pdb_ids:
            filepath = os.path.join("data", "%i.pdb"%id)
            active_sites.append(io.read_active_site(filepath))

        clusters = cluster.cluster_hierarchically(active_sites, len(pdb_ids))[0]

        # Should each be in their own cluster.
        self.assertClustersEqual(clusters, [["276"], ["10701"], ["4629"]])

    def test_hierarchical_clustering_one(self):
        # tractable subset
        pdb_ids = [276, 4629, 10701]

        active_sites = []
        for id in pdb_ids:
            filepath = os.path.join("data", "%i.pdb" % id)
            active_sites.append(io.read_active_site(filepath))

        clusters = cluster.cluster_hierarchically(active_sites, 1)[0]
        print(clusters)
        # Should all be in one cluster
        self.assertClustersEqual(clusters, [["276", "10701", "4629"]])

    def test_hierarchical_clustering_two(self):

        # tractable subset
        pdb_ids = [37438, 55996, 81859]

        # Try it four times, since test is flaky (because of inherent randomness)
        active_sites = []
        for id in pdb_ids:
            filepath = os.path.join("data", "%i.pdb" % id)
            active_sites.append(io.read_active_site(filepath))

        clusters = cluster.cluster_hierarchically(active_sites, 2)[0]

        # 37438 should be clustered with 81859
        self.assertClustersEqual(clusters, [["37438", "55996"], ["81859"]])
