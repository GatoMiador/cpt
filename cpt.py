#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python wrapper for the C++ implementation of the Conservative Power Theory
(CPT) in cpt.h.

https://ieeexplore.ieee.org/document/5638628

Created on Sun Sep  17 14:18:26 2023

@author: marcos

Example
-------
>>> import numpy as np, cpt
>>> c = cpt.CPT(channels=3, sampling_rate=12000, freq=60)
>>> r = c.feed_block(u, i)                  # u, i: arrays of shape (N, 3)
>>> r['totals']['P'][-1]                    # overall active power, last sample
>>> r = c.feed_block(u, i, fields=['t.P', 't.Q'], step=200)   # once per cycle
"""

import numpy as np

import _cpt

#: Per phase quantities, in output order.
PHASE_FIELDS = ('u', 'i', 'p', 'P', 'w', 'W', 'U', 'ia', 'ir', 'iv', 'iba', 'ibr', 'iua', 'iur')

#: Overall quantities, in output order (selected as 't.<name>' in 'fields').
TOTAL_FIELDS = ('P', 'Q', 'N', 'V', 'A', 'rf', 'uf', 'df', 'pf')

#: Quantities defined only for polyphase systems.
POLYPHASE_FIELDS = ('iba', 'ibr', 'iua', 'iur', 't.N', 't.uf')

_DTYPES = {'double': np.float64, 'float': np.float32}


def configurations():
    """Returns the available configurations as (phases, sampling_rate, freq, dtype) tuples."""
    return list(_cpt.configurations)


class CPT:
    """Computes the CPT quantities for a stream of voltage and current samples.

    Parameters
    ----------
    channels : int
        Number of phases (1 or 3).
    sampling_rate : int
        Sampling rate in Hz; must be one of the compiled configurations.
    freq : int
        Mains frequency in Hz (50 or 60).
    dtype : str
        'double' or 'float', the number type used by the C++ code.

    The averages use a moving window of one mains cycle
    (sampling_rate / freq samples), so the results settle after two cycles.
    """

    def __init__(self, channels=3, sampling_rate=61440, freq=60, dtype='double'):
        if dtype not in _DTYPES:
            raise ValueError(f"dtype must be one of {list(_DTYPES)}, not {dtype!r}")

        name = f'CPT_{channels}_{sampling_rate}_{freq}_{dtype[0]}'
        cls = getattr(_cpt, name, None)
        if cls is None:
            raise ValueError(
                f'configuration not compiled: channels={channels}, sampling_rate={sampling_rate}, '
                f'freq={freq}, dtype={dtype!r}; available: {configurations()}')

        self.channels = channels
        self.sampling_rate = sampling_rate
        self.freq = freq
        self.dtype = _DTYPES[dtype]
        self._core = cls()

    @property
    def samples(self):
        """Number of samples fed since construction or the last reset()."""
        return self._core.samples

    def available(self):
        """Names of the quantities this configuration computes ('t.' prefix for the overall ones)."""
        return self._core.available()

    def reset(self):
        """Clears the internal state, as if the object had just been created."""
        self._core.reset()

    def _as_input(self, x):
        a = np.asarray(x, dtype=self.dtype)
        if a.ndim == 1 and self.channels == 1:
            a = a.reshape(-1, 1)
        if a.ndim != 2 or a.shape[1] != self.channels:
            raise ValueError(f'expected shape (N, {self.channels}), got {np.shape(x)}')
        return np.ascontiguousarray(a)

    def feed_block(self, u, i, fields=None, step=1):
        """Feeds a block of samples and returns the results as NumPy arrays.

        Parameters
        ----------
        u, i : array_like
            Voltages and currents, shape (N, channels); shape (N,) is accepted
            for single phase.
        fields : list of str, optional
            Quantities to return, from PHASE_FIELDS and 't.' + TOTAL_FIELDS.
            Default: all quantities the configuration computes.
        step : int
            Keeps one result every 'step' samples, counted from the first
            sample fed (e.g. sampling_rate // freq for one result per cycle).
            The state is updated with every sample either way.

        Returns
        -------
        dict
            'n': sample number of each row; one array of shape (rows, channels)
            per phase quantity; 'totals': dict with one array of shape (rows,)
            per overall quantity.
        """
        return self._core.feed_block(self._as_input(u), self._as_input(i), fields, int(step))

    def feed(self, u, i):
        """Feeds one sample per phase and returns the results as lists and floats.

        Kept for compatibility with the former pure Python implementation;
        feed_block() is much faster for more than a few samples. Quantities
        not defined for single phase systems are returned as zeros.
        """
        u = np.asarray(u, dtype=self.dtype).reshape(1, self.channels)
        i = np.asarray(i, dtype=self.dtype).reshape(1, self.channels)
        r = self._core.feed_block(u, i)

        res = {k: r[k][0].tolist() if k in r else [0.0] * self.channels for k in PHASE_FIELDS}
        res['totals'] = {k: float(r['totals'][k][0]) if k in r['totals'] else 0.0 for k in TOTAL_FIELDS}
        return res


if __name__ == '__main__':
    import math

    cpt = CPT()
    n = np.arange(4096)
    angle = 0  # -math.pi/2
    phases = np.array([0, -2 * math.pi / 3, 2 * math.pi / 3])
    theta = 2 * math.pi * n[:, None] / 1024 + phases
    u = 127 * math.sqrt(2) * np.sin(theta)
    i = 1 * math.sqrt(2) * np.sin(theta + angle)

    res = cpt.feed_block(u, i)
    print({k: float(v[-1]) for k, v in res['totals'].items()})
