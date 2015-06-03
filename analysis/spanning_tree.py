import copy
import numpy as np
import pyregion
from astropy import units as u
from astropy import coordinates
from astropy.utils.console import ProgressBar
from scipy.sparse.csgraph import minimum_spanning_tree
import paths

regions = pyregion.open(paths.rpath('pointsource_centroids.reg'))

coords = coordinates.SkyCoord([x.coord_list for x in regions], unit=('deg','deg'))

separation_matrix = np.array([coords.separation(y) for y in ProgressBar(coords)])

ordered_separation_matrix = np.argsort(separation_matrix.ravel()).reshape(separation_matrix.shape)
ordered_separation_matrix[separation_matrix==0] = 0


ordered_edges = np.argsort(separation_matrix.ravel())
lefts, rights = np.unravel_index(ordered_edges, separation_matrix.shape)

separation_matrix[separation_matrix==0] = np.inf

#Prim's algorithm:
vertices = [0]
all_vertices = list(range(len(coords)))
all_vertices.remove(0)
edges = []
while len(vertices) < len(coords):
    smallest_edge = (-1,-1)
    smallest_length = np.inf
    for vert in vertices:
        separation = copy.copy(separation_matrix[vert, :])
        for v in vertices:
            separation[v] = np.inf 
        newvert = np.argmin(separation)
        assert newvert not in vertices
        assert newvert in all_vertices
        edge = (vert, newvert)
        length = separation[newvert]
        if length < smallest_length:
            smallest_length = length
            smallest_edge = edge
    edges.append(set(smallest_edge))
    vertices.append(smallest_edge[1])
    all_vertices.remove(smallest_edge[1])


mst = minimum_spanning_tree(separation_matrix)
connections = np.where(mst.toarray())
lines = [(coords[ii], coords[jj]) for ii,jj in zip(*connections)]
plotlines = [((a.ra.deg, b.ra.deg), (a.dec.deg, b.dec.deg)) for a,b in lines]
xx,yy = np.array(zip(*plotlines))

lines2 = [(coords[ll], coords[rr]) for ll,rr in edges]
plotlines2 = [((a.ra.deg, b.ra.deg), (a.dec.deg, b.dec.deg)) for a,b in lines2]
xx2,yy2 = np.array(zip(*plotlines2))

import pylab as pl
pl.clf()
pl.plot(coords.ra.deg, coords.dec.deg, '*')
pl.plot(xx.T, yy.T, color='b', alpha=0.5, linewidth=2)
pl.plot(xx2.T, yy2.T, color='r', alpha=0.5, linewidth=2)

inds = np.where(np.tril(separation_matrix, -1))
mps = separation_matrix[inds].mean()
mbl = mst[mst.nonzero()].mean()
print("Mean branch length: {0}".format(mbl))
print("Mean point separation: {0}".format(mps))
print("Q parameter: {0}".format(mbl/mps))
