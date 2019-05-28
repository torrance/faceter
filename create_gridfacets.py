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


def to_cartesian(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z


def to_spherical(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    return r, theta, phi


def rotate_about_x(x, y, z, angle):
    matrix = np.array([
        [1, 0, 0],
        [0, np.cos(angle), -np.sin(angle)],
        [0, np.sin(angle), np.cos(angle)],
    ])
    return np.matmul(matrix, [x, y, z])


def rotate_about_y(x, y, z, angle):
    matrix = np.array([
        [np.cos(angle), 0, np.sin(angle)],
        [0, 1, 0],
        [-np.sin(angle), 0, np.cos(angle)],
    ])
    return np.matmul(matrix, [x, y, z])


def rotate_about_z(x, y, z, angle):
    matrix = np.array([
        [np.cos(angle), -np.sin(angle), 0],
        [np.sin(angle), np.cos(angle), 0],
        [0, 0, 1],
    ])
    return np.matmul(matrix, [x, y, z])


def centroid(xs, ys):
    return (np.median(xs), np.median(ys))


def main(args):
    center = SkyCoord(args.center, unit=(units.hourangle, units.degree))

    model = fits.open(args.model)[0]
    wcs = WCS(header=model)

    idxs = np.indices(model.data.shape)
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

    facet_centers = []
    radius = np.radians(2 * args.radius)

    # See a point on uinit sphere radius away from origin
    r0, theta0, phi0 = 1, np.pi / 2  - radius, 0
    x0, y0, z0 = to_cartesian(r0, theta0, phi0)
    print("x0 y0 z0", x0, y0, z0)

    rotations = np.radians(range(0, 360, 45))
    for angle in rotations:
        x, y, z = rotate_about_x(x0, y0, z0, angle)
        print(x, y, z)
        facet_centers.append([x, y, z])

    facet_centers.append([1, 0, 0])
    facet_centers = np.array(facet_centers)
    print(facet_centers)

    # Now we've populated our facet centers, centered at zero, in Cartesian coords
    # Rotate to center coordinate
    theta = center.dec.rad
    phi = center.ra.rad
    print("Rotating by: ", np.degrees(theta))

    x, y, z = facet_centers.T
    print("Start: ", to_spherical(x, y, z))
    x, y, z = rotate_about_y(x, y, z, -theta)
    print("After theta rotation: ", to_spherical(x, y, z))
    x, y, z = rotate_about_z(x, y, z, phi)
    print("After phi rotation: ", to_spherical(x, y, z))

    r, theta, phi = to_spherical(x, y, z)

    facet_centers = SkyCoord(phi, np.pi / 2 - theta, unit=(units.radian, units.radian))
    with open('facet_centers.pkl', 'w') as f:
        pickle.dump(facet_centers, f)

    with open('facet_gridcenters.reg', 'w') as f:
        print("global color=red dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1", file=f)
        print("icrs", file=f)

        for i, center in enumerate(facet_centers):
            print("point %fd %fd # point=circle text={%d}" % (center.ra.deg, center.dec.deg, (i + 1)), file=f)

    print("Calculating dists array...")
    dists = []
    for center in facet_centers:
        dist = center.separation(grid_coords)
        dists.append(dist)
    dists = np.array(dists)

    closest = np.argmin(dists, axis=0)
    min_dist = np.min(dists, axis=0)
    closest[min_dist > args.max] = -1
    closest = closest.reshape(model.data.shape)

    # Create model images of each facet (for prediction)
    centroids = []
    for i, facet_center in enumerate(facet_centers):
        print("Calculating facet %d" % (i+1))
        facet_model = model.copy()
        facet_model.data[closest != i] = 0

        print(facet_model.data.sum(), len(facet_model.data.flatten()))

        # Calculate centroid
        idx = np.argwhere(closest == i)
        xp, yp = np.median(idx.T[2]), np.median(idx.T[3])
        ra, dec, freq, stokes = wcs.wcs_pix2world([[yp, xp, 0, 0]], 0)[0]
        print(ra, dec, freq, stokes)
        centroid = SkyCoord(ra, dec, unit=(units.degree, units.degree))
        centroids.append(centroid)

        # Calculate maximum distance from facet centroid
        max_dist = np.array(centroid.separation(grid_coords[closest.flatten() == i])).max()
        print("Max dist: ", max_dist)

        # Add metadata for later retrieval
        # Facet center
        facet_model.header['FCTCEN'] = centroid.to_string('hmsdms')
        # Maximum angular distance
        facet_model.header['FCTMAX'] = max_dist

        facet_model.writeto('facet-%d-%s-model.fits' % (i+1, args.channel), overwrite=True)

    with open('facet_centroids.reg', 'w') as f:
        print("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1", file=f)
        print("icrs", file=f)

        for i, center in enumerate(centroids):
            print("point %fd %fd # point=circle text={%d}" % (center.ra.deg, center.dec.deg, (i + 1)), file=f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', required=True)
    parser.add_argument('--channel', required=True)
    parser.add_argument('--max', type=float, default=999)
    parser.add_argument('--center', type=str, required=True)
    parser.add_argument('--radius', type=float, default=5)
    args = parser.parse_args()
    main(args)
