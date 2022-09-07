# network-distances

In this repository I explore one distance and one dissimilarity measure on two different complex networks. Using the shortest path length between two nodes as a distance measure is too sensitive to noise, whereas using a pure random walk distance measure has high variance. The dissimilarity measure described in [1] nicely generalises both of these.

## In this repository

- `distances.py`: in this file, I implemented the following two dissimilarity measures: 
  - `expected_rw_distance`: the distance function from this paper https://hal.inria.fr/hal-01944246/document (equations (10) and (17)).
  - `diss_distance`: the dissimilarity measure from [1].
- `thetas.ipynb`: `diss_distance` depends on a parameter $\theta$. In this file, I explored the influence of $\theta$ on the results.
- `toy-networks_distances-check.ipynb`: on toy networks, I randomly modify some edges and check the correlation between all the pairwise dissimilarities (only using `diss_distance`, as `expected_rw_distance` takes too long to compute) before and after modifications. This is to check whether the measure is robust under noise.
- the notebooks `ds1.ipynb`, `ds2.ipynb`, `ds2-run2.ipynb`, `ds2-run3.ipynb`, `ds2-different_thetas.ipynb` and `ds2-different_thetas-run2.ipynb` include experiments on applying the dissimilarity measures on two different datasets. For details, see the annotations in the notebooks.
- `network_tools.py` and `linalg_tools.py` are collections of short functions that I needed often.

## Conclusion and suggestions

`diss_distance` seems to be a promising dissimilarity measure to use on these kinds of networks. To further investigate this, one could do the same checks done in `ds1.ipynb` using different metadata. Furthermore, the parameter $\theta$ could be tuned to improve results.  
The implementation of `expected_rw_distance` is not optimal and needs to be improved, as it takes very long to compute (see also `toy-networks_distances-check.ipynb`). 

## References

[1] "A Family of Dissimilarity Measures between Nodes Generalizing both the Shortest-Path and the Commute-time Distances" by Luh Yen, Marco Saerens et al.
