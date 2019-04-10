#! /usr/bin/env python

import argparse

from casacore.tables import table, maketabdesc, makearrcoldesc


def main(eqn):
    # Check sanity of equation
    print(eqn)
    assert(len(eqn) == 3 or len(eqn) == 5)
    dest = eqn[0]
    assert(eqn[1] == '=')
    col1 = eqn[2]
    if len(eqn) == 5:
        assert(eqn[3] in ['+', '-'])
        op = eqn[3]
        col2 = eqn[4]
    else:
        col2 = None

    dest_mset = table(dest.split('::')[0], readonly=False)
    dest_col = dest.split('::')[1].upper()
    col1_mset = table(col1.split('::')[0])
    col1_col = col1.split('::')[1].upper()
    if col2:
        col2_mset = table(col2.split('::')[0])
        col2_col = col2.split('::')[1].upper()

    if col1_col not in col1_mset.colnames():
        print("column not recognised: ", col1)
        exit(1)
    if col2 and col2_col not in col2_mset.colnames():
        print("column not recognised: ", col2)
        exit(1)

    if dest_col not in dest_mset.colnames():
        col_dmi = dest_mset.getdminfo('DATA')
        col_dmi['NAME'] = dest_col
        shape = dest_mset.getcell('DATA', 0).shape
        dest_mset.addcols(
            maketabdesc(makearrcoldesc(dest_col, 0.0+0.0j, valuetype='complex', shape=shape)),
            col_dmi
        )
        dest_mset.flush()

    # Copy operation
    if len(eqn) == 3:
        data1 = col1_mset.getcol(col1_col)
        dest_mset.putcol(dest_col, data1)

    # Arithmetic operation
    if len(eqn) == 5:
        data1 = col1_mset.getcol(col1_col)
        data2 = col2_mset.getcol(col2_col)
        if op == '+':
            result = data1 + data2
        elif op == '-':
            result = data1 - data2

        dest_mset.putcol(dest_col, result)

    dest_mset.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("eqn", nargs='+')
    args = parser.parse_args()
    main(args.eqn)
