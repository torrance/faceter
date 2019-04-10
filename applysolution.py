#! /usr/bin/env python
from __future__ import print_function, division

import argparse
import sys

from casacore.tables import table, maketabdesc, makearrcoldesc, taql
import numpy as np
import numpy.linalg as linalg

import mwapy.aocal as aocal


def main(args):
    mset = table(args.mset, readonly=False)
    solution = aocal.fromfile(args.solution)
    solution = np.reshape(solution, (solution.shape[0], solution.shape[1], solution.shape[2], 2, 2))
    src = args.src.upper()
    dest = args.dest.upper()
    nrows = args.nrows


    if src not in mset.colnames():
        print("%s is not a valid column name" % src, file=sys.stderr)
        exit(1)

    if dest not in mset.colnames():
        col_dmi = mset.getdminfo('DATA')
        col_dmi['NAME'] = dest
        shape = mset.getcell('DATA', 0).shape
        mset.addcols(
            maketabdesc(makearrcoldesc(dest, 0.0+0.0j, valuetype='complex', shape=shape)),
            col_dmi
        )
        mset.flush()

    if args.reverse:
        solution = linalg.inv(solution)

    nchans, _ = mset.getcell(src, 0).shape
    width = nchans // solution.shape[2]
    idxs = np.array(range(0, nchans), dtype=int) // width
    lhs = solution[:, :, idxs]
    rhs = np.transpose(np.conj(solution), [0, 1, 2, 4, 3])[:, :, idxs]

    assert(solution.shape[0] == 1)  # Timesteps not handled yet
    assert(nchans % solution.shape[2] == 0)  # Width must be integer

    for n in range(0, len(mset) // nrows + 1):
        offset = n * nrows
        submset = taql("select * from $mset limit $nrows offset $offset")
        print("Processing rows %d - %d" % (offset, offset + len(submset)), file=sys.stderr)

        print("- Reading data from disk... ", end="", file=sys.stderr)
        sys.stderr.flush()
        ants1 = submset.getcol('ANTENNA1')
        ants2 = submset.getcol('ANTENNA2')
        data = submset.getcol(src)
        print("Done", file=sys.stderr)

        print("- Applying corrections... ", end="", file=sys.stderr)
        sys.stderr.flush()
        data = np.reshape(data, (data.shape[0], data.shape[1], 2, 2))
        corrected = np.matmul(lhs[0, ants1], np.matmul(data, rhs[0, ants2]))
        corrected = np.reshape(corrected, (data.shape[0], data.shape[1], 4))
        print("Done", file=sys.stderr)

        print("- Writing data back to disk... ", end="", file=sys.stderr)
        sys.stderr.flush()
        submset.putcol(dest, corrected)
        print("Done", file=sys.stderr)

    mset.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--src', default='DATA')
    parser.add_argument('--dest', default='CORRECTED_DATA')
    parser.add_argument('--reverse', action='store_true')
    parser.add_argument('--nrows', type=int, default=int(50e3))
    parser.add_argument('mset')
    parser.add_argument('solution')
    args = parser.parse_args()
    main(args)
