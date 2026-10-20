# -*- coding: utf-8 -*-
"""Tests of the Python wrapper (cpt.py + _cpt) itself: plumbing, API and
equivalence with the former pure Python implementation.

These tests do not check the CPT theory; they check that the wrapper exposes
exactly what cpt.h computes.

Run with: pip install .[test] && pytest
"""

import math

import numpy as np
import pytest

import cpt
import reference_cpt

FS, F = 12000, 60
CYCLE = FS // F


def three_phase(n, fs=FS, f=F):
    """Unbalanced and distorted three phase voltages and currents, with DC on the voltage."""
    t = np.arange(n) / fs
    a = 2 * np.pi * f * t[:, None] + np.array([0, -2 * np.pi / 3, 2 * np.pi / 3])
    u = 127 * np.sqrt(2) * np.sin(a) * np.array([1.0, 0.95, 1.05]) + 8 * np.sin(5 * a) + 0.5
    i = np.sqrt(2) * (np.array([10, 6, 3]) * np.sin(a - np.array([0.5, 0.2, -0.3])) + 2 * np.sin(3 * a + 0.4))
    return u, i


def single_phase(n, fs=FS, f=F, angle=-math.pi / 6):
    """127 V, 10 A sinusoids, current lagging by 'angle'."""
    t = np.arange(n) / fs
    u = 127 * np.sqrt(2) * np.sin(2 * np.pi * f * t)
    i = 10 * np.sqrt(2) * np.sin(2 * np.pi * f * t + angle)
    return u, i


def assert_results_equal(a, b):
    assert a.keys() == b.keys()
    for k in a:
        if k == 'totals':
            assert_results_equal(a[k], b[k])
        else:
            np.testing.assert_array_equal(a[k], b[k], err_msg=k)


def test_matches_former_python_implementation():
    n = 6 * CYCLE
    u, i = three_phase(n)

    ref = reference_cpt.CPT(channels=3, sampling_rate=FS, freq=F)
    expected = [ref.feed(list(u[k]), list(i[k])) for k in range(n)]
    got = cpt.CPT(3, FS, F).feed_block(u, i)

    for k in cpt.PHASE_FIELDS:
        np.testing.assert_allclose(got[k], [r[k] for r in expected], rtol=1e-12, atol=1e-9, err_msg=k)
    for k in cpt.TOTAL_FIELDS:
        np.testing.assert_allclose(got['totals'][k], [r['totals'][k] for r in expected],
                                   rtol=1e-12, atol=1e-9, err_msg='t.' + k)


def test_chunks_equal_single_block():
    u, i = three_phase(5 * CYCLE)
    whole = cpt.CPT(3, FS, F).feed_block(u, i)

    c = cpt.CPT(3, FS, F)
    cuts = [0, 1, 77, 200, 201, 650, len(u)]
    parts = [c.feed_block(u[a:b], i[a:b]) for a, b in zip(cuts, cuts[1:])]

    np.testing.assert_array_equal(np.concatenate([p['n'] for p in parts]), whole['n'])
    for k in cpt.PHASE_FIELDS:
        np.testing.assert_array_equal(np.concatenate([p[k] for p in parts]), whole[k])
    for k in cpt.TOTAL_FIELDS:
        np.testing.assert_array_equal(np.concatenate([p['totals'][k] for p in parts]), whole['totals'][k])


def test_step_keeps_every_nth_sample_across_chunks():
    u, i = three_phase(5 * CYCLE)
    whole = cpt.CPT(3, FS, F).feed_block(u, i)

    c = cpt.CPT(3, FS, F)
    parts = [c.feed_block(u[a:b], i[a:b], fields=['P', 't.P'], step=CYCLE)
             for a, b in [(0, 150), (150, 590), (590, len(u))]]

    n = np.concatenate([p['n'] for p in parts])
    np.testing.assert_array_equal(n, np.arange(CYCLE - 1, len(u), CYCLE))
    np.testing.assert_array_equal(np.concatenate([p['P'] for p in parts]), whole['P'][n])
    np.testing.assert_array_equal(np.concatenate([p['totals']['P'] for p in parts]), whole['totals']['P'][n])


def test_fields_selects_outputs():
    u, i = three_phase(CYCLE)
    r = cpt.CPT(3, FS, F).feed_block(u, i, fields=['ia', 't.pf'])
    assert set(r) == {'n', 'ia', 'totals'}
    assert set(r['totals']) == {'pf'}

    r = cpt.CPT(3, FS, F).feed_block(u, i, fields=['ia'])
    assert set(r) == {'n', 'ia'}


def test_feed_equals_feed_block():
    u, i = three_phase(2 * CYCLE + 3)
    block = cpt.CPT(3, FS, F).feed_block(u, i)

    c = cpt.CPT(3, FS, F)
    for k in range(len(u)):
        r = c.feed(u[k], i[k])
    for k in cpt.PHASE_FIELDS:
        assert r[k] == block[k][-1].tolist()
    for k in cpt.TOTAL_FIELDS:
        assert r['totals'][k] == block['totals'][k][-1]
    assert c.samples == len(u)


def test_reset_restores_initial_state():
    u, i = three_phase(3 * CYCLE)
    c = cpt.CPT(3, FS, F)
    first = c.feed_block(u, i)
    c.feed_block(*three_phase(CYCLE))
    c.reset()
    assert c.samples == 0
    assert_results_equal(c.feed_block(u, i), first)


def test_single_phase():
    c = cpt.CPT(1, FS, F)
    assert not set(cpt.POLYPHASE_FIELDS) & set(c.available())

    u, i = single_phase(4 * CYCLE)
    r = c.feed_block(u, i, step=CYCLE)          # 1-D input accepted
    assert r['P'].shape == (4, 1)
    assert r['totals']['P'][-1] == pytest.approx(1270 * math.cos(math.pi / 6), rel=1e-4)
    assert r['totals']['Q'][-1] == pytest.approx(1270 * math.sin(math.pi / 6), rel=1e-4)

    # Polyphase quantities are zero in the compatibility API, never garbage
    r = c.feed([u[0]], [i[0]])
    for k in ('iba', 'ibr', 'iua', 'iur'):
        assert r[k] == [0.0]
    assert r['totals']['N'] == 0.0
    assert r['totals']['uf'] == 0.0


def test_float_configuration():
    u, i = three_phase(4 * CYCLE)
    d = cpt.CPT(3, FS, F, dtype='double').feed_block(u, i)
    f = cpt.CPT(3, FS, F, dtype='float').feed_block(u, i)
    assert f['P'].dtype == np.float32
    for k in ('P', 'Q', 'A'):
        np.testing.assert_allclose(f['totals'][k][-1], d['totals'][k][-1], rtol=1e-4, err_msg=k)


@pytest.mark.parametrize('phases, fs, f, dtype', cpt.configurations())
def test_every_configuration(phases, fs, f, dtype):
    n = 3 * fs // f
    t = np.arange(n) / fs
    a = 2 * np.pi * f * t[:, None] - np.arange(phases) * 2 * np.pi / 3
    u = 100 * np.sqrt(2) * np.sin(a)
    i = 2 * np.sqrt(2) * np.sin(a)
    r = cpt.CPT(phases, fs, f, dtype).feed_block(u, i, fields=['t.P'], step=fs // f)
    assert r['totals']['P'][-1] == pytest.approx(200 * phases, rel=1e-4)


def test_errors():
    c = cpt.CPT(3, FS, F)
    u, i = three_phase(10)
    with pytest.raises(ValueError, match='shape'):
        c.feed_block(u[:, :2], i[:, :2])
    with pytest.raises(ValueError, match='same number'):
        c.feed_block(u, i[:5])
    with pytest.raises(ValueError, match='unknown field'):
        c.feed_block(u, i, fields=['X'])
    with pytest.raises(ValueError, match='step'):
        c.feed_block(u, i, step=0)
    with pytest.raises(ValueError, match='single phase'):
        cpt.CPT(1, FS, F).feed_block(u[:, 0], i[:, 0], fields=['t.N'])
    with pytest.raises(ValueError, match='not compiled'):
        cpt.CPT(3, 44100, 60)
    with pytest.raises(ValueError, match='dtype'):
        cpt.CPT(3, FS, F, dtype='int')
