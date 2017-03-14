"""Time utilities: time sample from interval, date conversion etc.
"""
from datetime import datetime

import numpy as np


def second_from_interval(start, stop, n=1):
    """Sample random times from an interval (in seconds)."""
    sec = np.timedelta64(1, 's')
    n_seconds = (stop - start) / sec
    samples = np.random.randint(low=0, high=n_seconds, size=n) * sec
    return start + samples


def equidistant_from_interval(start, stop, step=np.timedelta64(2, 'm')):
    """Draw equidistant samples (fixed stepsize) from interval."""
    start = np.datetime64(start)
    stop = np.datetime64(stop)
    duration = stop - start
    n_steps = np.ceil(duration / step)
    samples = np.arange(n_steps) * step
    return start + samples


def random_date(year=2015, **randargs):
    """Create random dates in the given year.

    Parameters
    ----------
    year: int, default: 2015
    n: int, default: 1
    """
    start = np.datetime64('{}-01-01'.format(year))
    stop = np.datetime64('{}-01-01'.format(year + 1))
    return second_from_interval(start, stop, **randargs)


def np_to_datetime(intime):
    """Convert numpy/pandas datetime64 to list[datetime]."""
    nptime = np.atleast_1d(intime)
    np_corr = (nptime - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
    return [datetime.utcfromtimestamp(t) for t in np_corr]
