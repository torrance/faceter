#! /usr/bin/env python
from __future__ import print_function, division

import argparse
import pickle
import sys

from astropy.coordinates import SkyCoord
import astropy.units as units
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np


def main(args):
    image = fits.open(args.image)[0]
    wcs = WCS(header=image)

    idxs = np.indices(image.data.shape)
    # WCS expects indices in header order (ie. reversed)
    idxs = np.flip(idxs, axis=0)
    # Rearrange indices to be of form N, D
    idxs = np.transpose(idxs, axes=[1, 2, 3, 4, 0])
    idxs = np.reshape(idxs, [-1, 4])

    # Calculate world coordinates of grid
    print("Calculating world coordinates of image grid... ", end="", file=sys.stderr)
    grid_coords = wcs.wcs_pix2world(idxs, 0)
    grid_coords = SkyCoord(grid_coords[:, 0], grid_coords[:, 1], unit=(units.degree, units.degree))
    print("Done", file=sys.stderr)

    # Use existing facet centers
    with open('facet_centers.pkl') as f:
        facet_centers = pickle.load(f)

    # Calculate distance of each pixle to each of facet_centers
    print("Calculating pixel distances to facet centers... ", end="", file=sys.stderr)
    dists = []
    for facet_center in facet_centers:
        dist = facet_center.separation(grid_coords)
        dists.append(dist)

    dists = np.array(dists)
    print("Done", file=sys.stderr)

    closest = np.argmin(dists, axis=0)
    min_dist = np.min(dists, axis=0)
    closest[min_dist > args.max] = -1
    closest = closest.reshape(image.data.shape)

    # Create clean mask for facet
    image.data[:] = 0
    image.data[closest == (args.facetid - 1)] = 1
    image.writeto('facet-%d-mask.fits' % args.facetid, overwrite=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--max', default=4.0, type=float, help='Maximum facet radius (degrees)')
    parser.add_argument('--facetid', type=int, required=True)
    parser.add_argument('--image', required=True)
    args = parser.parse_args()
    main(args)
