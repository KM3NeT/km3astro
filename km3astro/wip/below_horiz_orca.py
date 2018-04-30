#!/usr/bin/env python
"""Test whether at some time the sun is below the horizon in ORCA."""


with pd.HDFStore('orca_sun_isup.h5') as h5:
    o_rise = h5['o_rise'].values
    o_set = h5['o_set'].values


def above_horiz(times):
    return [(t >= o_rise) & (t <= o_set) for t in times]


def is_below(times):
    ab = np.array(above_horiz(times))
    risen = np.logical_not(np.sum(ab, axis=1))
    return risen
