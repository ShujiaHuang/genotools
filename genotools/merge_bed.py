"""Merge all overlap regions

Author: Shujia Huang
Date: 2021-12-17 09:59:04
"""
import sys


def merge_region(position_region, delta=1):
    """Merge a batch of sorted region

    Parameters
    ----------
    ``position_region``: a list like, required
        A regions (2D) array, format like: [[start1,end1], [start2,end2], ...]

    ``delta``: Integer, optinal

    Example
    -------

    >>> merge_region([[1,1], [2,3], [4,6], [4,5], [8, 20], [9, 12]])
    ... [[1, 6], [8, 20]]

    """
    # sorted
    position_region.sort(key=lambda x: x[0])

    m_region = []
    pre_pos, start, end = None, None, None

    flag = False
    for s, e in position_region:
        if s > e:
            sys.stderr.write("[ERROR]Your region start > end. It is not allow when call Merge function\n")
            sys.exit(1)

        if pre_pos is None:
            # The light is on => Get the region!
            if flag:
                m_region.append([start, end])

            start, end = s, e
            flag = True

        else:
            if pre_pos > s:
                sys.stderr.write("[ERROR] The array hasn't been sorted.\n")
                sys.exit(1)

            if delta + end >= s:
                end = e if e > end else end

            else:
                m_region.append([start, end])
                start, end = s, e

        pre_pos = s

    if flag:
        m_region.append([start, end])

    return m_region


if __name__ == "__main__":
    in_fname = sys.argv[1]

    regions = {}
    with open(in_fname) as IN:

        for line in IN:
            """
            chr1	1577362	1877362
            chr1	1945970	2245970
            chr1	2055646	2355646
            """
            col = line.strip().split()
            if col[0] not in regions:
                regions[col[0]] = []

            regions[col[0]].append([int(col[1]), int(col[2])])

    for k, region in regions.items():
        for d in merge_region(region):
            print("%s\t%d\t%d" % (k, d[0], d[1]))




