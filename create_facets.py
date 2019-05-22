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
    model = fits.open(args.model)[0]
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

    try:
        with open('facet_centers.pkl') as f:
            facet_centers = pickle.load(f)
        dists = np.load('dists.npy')
    except IOError as e:
        facet_centers = []
        dists = []

        while np.nanmax(image.data) > args.threshold:
            brightest = np.unravel_index(np.nanargmax(image.data), image.data.shape)
            coord = wcs.wcs_pix2world([np.flip(brightest)], 0)
            coord = SkyCoord(coord[0, 0], coord[0, 1], unit=(units.degree, units.degree))
            facet_centers.append(coord)
            print("\rFacet count: %d... " % len(facet_centers), end="", file=sys.stderr)

            dist = coord.separation(grid_coords)
            dists.append(dist)
            nearby_idx = (dist < args.min * units.degree).reshape(image.data.shape)
            image.data[nearby_idx] = np.nan

        with open('facet_centers.pkl', 'w') as f:
            pickle.dump(facet_centers, f)

        dists = np.array(dists)
        np.save('dists.npy', dists)
        print("Done", file=sys.stderr)

    dists = np.array(dists)

    closest = np.argmin(dists, axis=0)
    min_dist = np.min(dists, axis=0)
    closest[min_dist > args.max] = -1
    closest = closest.reshape(image.data.shape)

   #  Create model images of each facet (for prediction)
    for i, facet_center in enumerate(facet_centers):
        facet_model = model.copy()
        facet_model.data[closest != i] = 0

        # Add metadata for later retrieval
        # Facet center
        facet_model.header['FCTCEN'] = facet_center.to_string('hmsdms')
        # Maximum angular distance
        facet_model.header['FCTMAX'] = min_dist[closest.flatten() == i].max()

        facet_model.writeto('facet-%d-%s-model.fits' % (i+1, args.channel), overwrite=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--threshold', type=float, default=3.0)
    parser.add_argument('--min', default=0.25, type=float, help='Minium facet radius (degrees)')
    parser.add_argument('--max', default=4.0, type=float, help='Maximum facet radis (degrees)')
    parser.add_argument('--image', required=True)
    parser.add_argument('--model', required=True)
    parser.add_argument('--channel', required=True)
    args = parser.parse_args()
    main(args)
