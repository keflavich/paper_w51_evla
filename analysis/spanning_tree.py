import copy
import numpy as np
import pyregion
from astropy import units as u
from astropy import coordinates
from astropy.utils.console import ProgressBar
from scipy.sparse.csgraph import minimum_spanning_tree
import paths
from common_constants import distance

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


def examinalyze(smat):
    mst = minimum_spanning_tree(smat)
    connections = np.where(mst.toarray())
    lines = [(coords[ii], coords[jj]) for ii,jj in zip(*connections)]
    plotlines = [((a.ra.deg, b.ra.deg), (a.dec.deg, b.dec.deg)) for a,b in lines]
    xx,yy = np.array(list(zip(*plotlines)))

    lines2 = [(coords[ll], coords[rr]) for ll,rr in edges]
    plotlines2 = [((a.ra.deg, b.ra.deg), (a.dec.deg, b.dec.deg)) for a,b in lines2]
    xx2,yy2 = np.array(list(zip(*plotlines2)))

    import pylab as pl
    pl.clf()
    pl.plot(coords.ra.deg, coords.dec.deg, '*')
    pl.plot(xx.T, yy.T, color='b', alpha=0.5, linewidth=2)
    pl.plot(xx2.T, yy2.T, color='r', alpha=0.5, linewidth=2)

    inds = np.where(np.tril(smat, -1))
    mps = smat[inds].mean()
    mbl = mst[mst.nonzero()].mean()
    print("Mean branch length: {0} arcsec {1}".format(mbl*3600,
                                                         (mbl*u.deg*distance).to(u.pc,
                                                                                 u.dimensionless_angles())))
    print("Mean point separation: {0} arcsec {1}".format(mps*3600,
                                                            (mbl*u.deg*distance).to(u.pc,
                                                                                    u.dimensionless_angles())))
    print("Q parameter: {0}".format(mbl/mps))

examinalyze(separation_matrix)


w51e2 = coordinates.SkyCoord('19 23 43.90 +14 30 34.8', unit=('hour', 'deg'), frame='fk5')
coords_e1e2 = coords[coords.separation(w51e2) < 15*u.arcsec]

separation_matrix_e1e2 = np.array([coords_e1e2.separation(y)
                                   for y in ProgressBar(coords_e1e2)
                                  ])
print("W51e1e2 cluster: ")
examinalyze(separation_matrix_e1e2)

w51e1 = coordinates.SkyCoord('19 23 43.77 +14 30 25.9', unit=('hour', 'deg'), frame='fk5')
coords_e1 = coords[coords.separation(w51e1) < 5*u.arcsec]

separation_matrix_e1 = np.array([coords_e1.separation(y)
                                   for y in ProgressBar(coords_e1)
                                  ])
print("W51e1 cluster: ")
examinalyze(separation_matrix_e1)
