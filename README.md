# HW2 Skeleton

[![Build Status](https://travis-ci.org/christacaggiano/Unsupervised-clustering.svg?branch=master)](https://travis-ci.org/christacaggiano/Unsupervised-clustering)

# Use

`python clustering_algorithms/main.py [-P| -H | -T] <pdb directory> <output file>"`

* -P: k-means clustering
* -H: agglomerative hierarchical clustering
* -T: compare clustering methods


# Hierarchical and Partitioning Clustering  

## Similarity Metric
I chose to use a Damerau-Levenshtein edit distance<sup>[1](#PythonDL)</sup> algorithm to calculate the distance between two amino acid sequences. This algorithm quantifies the 'distance' between two strings by counting the number of operations that are needed to transform one string into another. In this algorithm the operations allowed are inserting a letter, deleting a letter, substituting a letter or transposing a letter. This is an improvement on the Levenshtein distance algorithm that proceeded it as transpositions are common errors, both in the written language, and may have a viable biological relevance.<sup>[2](#Explanation)</sup>


I decided that sequence similarity was the best indication of similarity between proteins as this naturally encompasses
some distinguishing properties of proteins- such as charge or hydrophobicity. Sequence similarity also is useful for
measuring how evolutionarily close two proteins are. Since we often think of proteins that are closer together in
evolutionary time as more similar, this would be a useful feature of my metric. For example, two proteins may have
similar levels of hydrophobicity in their active sites, but have some from extremely divergent evolutionary histories,
which may not be useful in assessing the functional similarities between proteins.

## Partitioning Algorithm
For a partitioning algorithm, I chose a k-means clustering algorithm that calculated centroids by taking the “average”
of all the sequences in a given cluster.

K-means clustering involves choosing *k* initial clusters. *K* is defined by the user, and generally has a broader meaning in your data. For example, if a biologist was clustering RNA-seq expression data from 3 cell types- like red blood cells, white blood cells, and platelets- the biologist may specify 3 initial clusters. The number of clusters defines the number of starting points for the algorithm, called centroids. These starting points may be chosen randomly, as is implemented this code, or may try to intelligently pick centroids depending on variables such as the mean or spread of the data.<sup>[3](#Fancyclustering)</sup> All data points are placed into a cluster by a similarity metric, in my case the Damerau-Levenshtein edit distance described above.

Once all data points are in a cluster, a new centroid is chosen. In my implementation, the "average" string was taken as the new centroid. Average was defined as the string that had the smallest distance to all other strings in the cluster.
Since I was working with strings
and not integers, I felt this was the fairest way to calculate an average. Averaging in this way, however, could be a
drawback when the sequences in a cluster are very distant. Additionally, this constrains the centroid to be included in
my set of sequences.

Then, the distance of all data points to the new centroids is calculated. If a data point is closer to the centroid of another cluster, it is deleted from its cluster, and swapped to the new cluster. This is repeated until no swaps are made.<sup>[4](#kmeansexplanation)</sup>


## Hierarchical Clustering Algorithm

Hierarchical clustering takes unordered data points, and given a similarity metric, puts them into a ranked order. Either a user can specify a number of clusters, as in partitioning clustering, or allow the algorithm to proceed until every data item is ordered.


My hierarchical clustering was a simple agglomerative clustering with a single linkage. This began by assigning each data point to its own cluster. An upper diagonal *n x n* matrix is computed, where the elements of the matrix are the similarity between each cluster and all others in the list.

The distance between two clusters was determined by single linkage, or the minimum distance between two clusters. Single linkage was chosen
because it seemed to make the most sense for the strings I was working with.<sup>[5](#linkage)</sup> Average or centroid linkage was difficult
to compute for a set of strings, and complete linkage was very biased toward outlying strings. The most similar clusters were then merged. This continues until we have *k* clusters, or until all items are merged.<sup>[6](#linkage)</sup>


I found that this was a
straightforward algorithm, without needing a specified number of clusters or being dependent on initial starting conditions,
as with k-means clustering. However, I felt that this algorithm was less biologically relevant than k-means as at the end,
I had a difficult time interpreting what the hierarchy of sequences actually means.


## Quality Metric
The most intuitive way I thought of to calculate a quality metric was by generating sequences that were extremely similar
– off by one or two amino acids- and performing multiple simulations to observe how often they clustered together. I found
that this metric was very sensitive in k-means (Image 1) where the initial starting conditions and number of clusters changed
the number of successes. In hierarchical clustering, because it was not dependent on random initial conditions the results
were more consistent (Image 2). Thus, this metric was less useful for hierarchical clustering. This metric could have been improved
by having a more robust test dataset, many test datasets, by using cross-validation, or a combination of all three.

![alt text][image1]
**Image 1:** Number of times (out of 10 simulations) that two sets of nearly identical sequences
were clustered together in my partitioning algorithm

![alt text][image2]
**Image 2:** The ranking of similar molecules by hierarchical algorithm.


# Comparisons of Methods

I computed the Jacard similarity to compare *k=3* clusters yielded by hierarchical and partitioning clustering. The Jacard similarity is the intersection of all the items in one cluster divided by its union.<sup>[7](#jacard)</sup>



|  |H1         | H2| H3       |
|--|----------:|--:|---------:|
|**P1** | 0         | 0 | 0.467742 |
| **P2** |0         | 1 | 0        |
|**P3**| 0.0294118 | 0 | 0.52381  |

**Table 1**: Jacard similarity of clusters produced by partitioning and hierarchical clustering methods.

Overall, the clusters were not very similar (Table 1). The exception was cluster #2 from each method was nearly identical, suggesting there is some shared importance. Since my partitioning algorithm was initialized randomly, however, so this similarity is not stable (Table 2).

|  |H1         | H2| H3       |
|--|----------:|--:|---------:|
|**P1**| 0        | 0.0384615 | 0.78125  |
|**P2** | 0.333333 | 0         | 0.031746 |
|**P3** | 0        | 0         | 0.16129  |


**Table 2**: Another iteration of partitioning algorithm.

##  Biological Relevance

I tested several different cluster sizes and examined the properties of the proteins in each cluster. Using the BioPython Protein Analysis method
I calculated the hydrophobicity and isoelectric point.<sup>[8](#biopython)</sup>

For k-means clustering, I plotted the hydrophobicity against the isolectric point for clusters of size 4, 3, and 2.

![alt text][image3]
**Image 3:** isolectric point vs average hydrophobicity K=4

![alt text][image4]
**Image 4:** isolectric point vs average hydrophobicity K=3

![alt text][image5]
**Image 5:** isolectric point vs average hydrophobicity K=2


I found these clusters to have little biological meaning. This is no clear grouping with the clusters determined by the algorithm. My algorithm's inability to predict clusters based on these properties, however, is not surprising. My similarity metric only penalizes differences between amino acids. It doesn't take in account the biological properties that could be shared between amino acids. For example, two hydrophobic amino acids, like glycine or alanine, would be penalized to the same degree as amino acids with different hydrophobicity, like glycine and glutamine. This would be an important thing to consider in future iterations of this clustering algorithm, however, since charge and hydrophobicity are both important elements of the protein active site's functionality.

However, I found that my clusters were well predicted by molecular weight.

![alt text][image7]
**Image 6:** molecular weight vs isoelectric point K=4

![alt text][image6]
**Image 7:** molecular weight vs isoelectric point K=4

This may be a result of some clusters being heavily composed of very light amino acids, such as glycine in the red cluster, whereas heavier amino acids were rare amino acids- such as tryptophan and thus, were immediately more dissimilar.

I repeated the same assessment for biological relevance with hiearchical clustering. I again tested the distribution of isolectric points and molecular weight for my clusters.

![alt text][image8]
**Image 8:** Distribution of isolectric point and molecular weight for all proteins.

![alt text][image9]
**Image 9:** Distribution of isolectric point and molecular weight for K=3

![alt text][image10]
**Image 10:** Distribution of isolectric point and molecular weight for K=2

In general, it seems that isoelectric point was approximately the same for all clusters, whereas molecular weight varied. This could be quantified further.

Another insight gleamed from this assessment is that my algorithm tended to isolate proteins into a solitary cluster. This could be an opportunity for future versions of this algorithm to learn what makes that protein so distinct, and use that information to better separate clusters.

I found kmeans clustering to be fairly biologically relevant. It tended to group active sites with similar motifs
together and grouped the sequences of similar lengths. Both of these things may affect protein functionality.
Hierarchical clustering was more consistent than k-means, as pairs of simulated sequences that were expected to be close
were always adjacent in the hierarchy (Image 2). However, I found this difficult to interpret with biological meaning.
What does it mean for molecule D and molecule A to be nearby in a hierarchy? In this case, I think exploring other
similarity metrics that more easily lend to numerical ranking, such as hydrophobicity or charge, would be better suited
for hierarchical clustering.



 # Bibliography
 <a name="PythonDL">1</a>: https://github.com/gfairchild/pyxDamerauLevenshtein

 <a name="Explanation">2</a>: https://www.mathworks.com/matlabcentral/cody/problems/2309-calculate-the-damerau-levenshtein-distance-between-two-strings

 <a name="Fancyclustering">3</a>: http://people.csail.mit.edu/tieu/notebook/kmeans/15_p600-hamerly.pdf

 <a name="kmeansexplanation">4</a>: http://bigdata-madesimple.com/possibly-the-simplest-way-to-explain-k-means-algorithm/

 <a name="linkage">5</a>: http://www.saedsayad.com/clustering_hierarchical.htm

 <a name="similarity">6</a>: http://www.analytictech.com/networks/hiclus.htm

 <a name="jacard">7</a>: http://www.statisticshowto.com/jaccard-index/

 <a name="biopython">8</a>: http://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParam.ProteinAnalysis-class.html



[image1]: https://raw.githubusercontent.com/christacaggiano/images/master/cluster_quality.png "Image1"

[image2]: https://raw.githubusercontent.com/christacaggiano/images/master/cluster_quality2.png "Image2"

[image3]: https://raw.githubusercontent.com/christacaggiano/images/master/Figure_1.png "Image3"

[image4]: https://raw.githubusercontent.com/christacaggiano/images/master/Figure_2.png "Image4"

[image5]: https://raw.githubusercontent.com/christacaggiano/images/master/Figure_3.png "Image5"

[image6]:
https://raw.githubusercontent.com/christacaggiano/images/master/fig5.png "Image 6"

[image7]:
https://raw.githubusercontent.com/christacaggiano/images/master/fig6.png "Image 7"

[image8]:
https://raw.githubusercontent.com/christacaggiano/images/master/hierarch1.png "Image 8"

[image9]:
https://raw.githubusercontent.com/christacaggiano/images/master/hierarch2.png "Image 9"

[image10]:
https://raw.githubusercontent.com/christacaggiano/images/master/hierarch3.png "Image 8"
