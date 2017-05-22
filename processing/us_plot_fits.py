#!/usr/bin/env python
# Copyright 2017 Tom Eulenfeld, MIT license
"""
Reprocess some events to plot envelope fits.
"""

from copy import deepcopy
import json
import numpy as np
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import pickle

from us_qopen import CONF, EVENTRESULTSJSON, calc_stuff
from qopen.core import get_pair, Gsmooth
from qopen.rt import G as G_func


def set_gridlabels(ax, i, n, N, xlabel='frequency (Hz)', ylabel=None):
    if i % n != 0 and ylabel:
        plt.setp(ax.get_yticklabels(), visible=False)
    elif i // n == (n - 1) // 2 and ylabel:
#        ax.set_ylabel(ylabel)
        ax.annotate(ylabel, (0, 1), (-35, 0), 'axes fraction', 'offset points',
                    rotation='vertical', clip_on=False, va='center', ha='right')
    if i < N - n and xlabel:
        plt.setp(ax.get_xticklabels(), visible=False)
    elif i % n == (n - 1) // 2 and i >= N - n - 1 and xlabel:
        ax.set_xlabel(xlabel)
#        ax.annotate(xlabel, (1, 0), (5, -20), 'axes fraction', 'offset points',
#                    clip_on=False, va='top', ha='center')
    if i % n == 0:
        for tick in ax.get_yaxis().get_major_ticks():
            tick.set_pad(0.)
            tick.label1 = tick._get_text1()
        for i, l in enumerate(ax.yaxis.get_ticklabels()):
            if i % 2 == 0:
                l.set_visible(False)


def _get_times(tr):
    t0 = tr.stats.starttime - tr.stats.origintime
    return np.arange(len(tr)) * tr.stats.delta + t0


def plot_fits(energies, g0, b, W, R, v0, info, smooth=None,
              smooth_window='bartlett', title=None, mag=None):
    fs = 150 / 25.4
    plt.figure(figsize=(fs, 0.6*fs))
    tcoda, tbulk, Ecoda, Ebulk, Gcoda, Gbulk = info
    N = len(energies)
    nx, ny = 3, 2
    gs = gridspec.GridSpec(ny, nx, wspace=0.05, hspace=0.05)
    share = None
    if b is None:
        b = 0
    tmaxs = []
    ymaxs = []
    ymins = []
    c1 = 'mediumblue'
    c2 = 'darkred'
    c1l = '#8181CD'
    c2l = '#8B6969'
    for i, energy in enumerate(energies):
        evid, station = get_pair(energy)
        ax = plt.subplot(gs[i // nx, i % nx], sharex=share, sharey=share)
        plot = ax.semilogy

        def get_Emod(G, t):
            return R[station] * W[evid] * G * np.exp(-b * t)
        st = energy.stats
        r = st.distance
        t = _get_times(energy) + r / v0 - (st.sonset - st.origintime)

        if smooth:
            plot(t, energy.data_unsmoothed, color='0.7')
        plot(t, energy.data, color=c1l)
        G_ = Gsmooth(G_func, r, t, v0, g0, smooth=smooth,
                     smooth_window=smooth_window)
        Emod = get_Emod(G_, t)
        index = np.argwhere(Emod < 1e-30)[-1]
        Emod[index] = 1e-30

        plot(t, Emod, color=c2l)

        plot(tcoda[i], Ecoda[i], color=c1)
        Emodcoda = get_Emod(Gcoda[i], tcoda[i])
        plot(tcoda[i], Emodcoda, color=c2)

        if tbulk and len(tbulk) > 0:
            plot(tbulk[i], Ebulk[i], 'o', color=c1, mec=c1, ms=4)
            Emodbulk = get_Emod(Gbulk[i], tbulk[i])
            plot(tbulk[i], Emodbulk, 'o', ms=3,
                 color=c2, mec=c2)

        l = '%s\nr=%dkm' % (station, r / 1000)
        ax.annotate(l, (1, 1), (-5, -5), 'axes fraction',
                    'offset points', ha='right', va='top', size='x-small')
        set_gridlabels(ax, i, nx, N, xlabel='time (s)',
                       ylabel=r'E (Jm$^{-3}$Hz$^{-1}$)')
        kw = dict(color='darkgreen', alpha=0.5, lw=0, zorder=10000)
        ax.axvspan(tcoda[i][0]-10, tcoda[i][0]-0.8, 0.05, 0.08, **kw)
        ax.axvspan(tcoda[i][0]+0.8, tcoda[i][-1], 0.05, 0.08, **kw)
        tmaxs.append(t[-1])
        ymaxs.append(max(np.max(Emod), np.max(energy.data)))
        ymins.append(min(np.min(Emodcoda), np.max(Ecoda[i])))
        if share is None:
            share = ax
    f = 6
    lmag = (r'  $M_{\mathrm{R}}%.1f$' % mag) if mag else ''
    l = ('id %s' + lmag + '\n'
         r'$Q^{-1}_{\mathrm{sc}} = %.1f\times 10^{%d}$' + '\n'
         r'$Q^{-1}_{\mathrm{i}} = %.1f\times 10^{%d}$')
    l = l % (evid, g0 * v0 / (2 * np.pi * f) * 1e3, -3,
             b / (2 * np.pi * f) * 1e3, -3)
    ax.annotate(title, (0.715, 0.43), xycoords='figure fraction',
                ha='left', va='top', annotation_clip=False, size='large')
    ax.annotate(l, (0.715, 0.34), xycoords='figure fraction',
                ha='left', va='top', annotation_clip=False)
    ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.set_yticks(10. ** np.arange(-18, -7, 2))
    ax.set_xlim((-5, 120))
    ax.set_ylim((1e-15 / 1.5, 1e-7 * 1.5))


if __name__ == '__main__':
    with open(EVENTRESULTSJSON) as f:
        events = json.load(f)
    iter_ = [
        ('CI.BFS', 34.2388, -117.6585, 0.3, 3),
        ('TA.S53A', 37.6815, -82.1264, 0.08, 2.5)
        # ('BK.WENL', 37.6221, -121.7570, 1, 3),
        # ('TA.I16A', 43.8756, -111.4868, 1.2, 3),
        # ('TA.R54A', 38.1909, -80.9904, 1, 3),
        # ('TA.T48A', 37.1094, -86.3943, 0.6, 2.5)
        ]
    conf = deepcopy(CONF)
    conf['freqs'] = {"width": 1, "cfreqs": [3.0, 6.0, 12.0]}
    conf['parallel'] = False
    conf['logfile'] = '../tmp/us_qopen_plot_fits.log'
    conf['plot_optimization'] = True
    conf['plot_fits'] = True
    conf['plot_eventresult'] = True
    # start Qopen again for some events and plot envelope fits and optimizations
    for sta, lat, lon, dif, mag in iter_:
        arr = np.array
        ind = ((((arr(events['lon']) - lon) * np.cos(lat / 180 * np.pi)) ** 2 +
                (arr(events['lat']) - lat) ** 2) ** 0.5 < dif) & (
                arr(events['MR']) > mag)
        ids = arr(events['id'])[ind]
        print('Selected %d events near station %s.' % (len(ids), sta))
        conf['plot_optimization_options']['fname'] = '../tmp/plot_fits/' + sta + '/%s_opt.png'
        conf['plot_fits_options']['fname'] = '../tmp/plot_fits/' + sta + '/%s_fits.png'
        conf['plot_eventresult_options']['fname'] = '../tmp/plot_fits/' + sta + '/%s_res.png'
        calc_stuff(conf, eventids=ids)

    # make a nice plot with envelope fits for two events
    conf = deepcopy(CONF)
    conf['freqs'] = {"width": 1, "cfreqs": [6.0]}
    conf['parallel'] = False
    conf['logfile'] = '../tmp/us_qopen_plot_fits.log'
    conf['dump_fitpkl'] = '../tmp/plot_fits/%s_fit.pkl'
    calc_stuff(conf, eventids=('2511938', '4808371'))
    annotate_kw = dict(xy=(0.02, 0.89), xycoords='figure fraction',
                       fontsize='large', annotation_clip=False)
    fname = '../tmp/plot_fits/2511938_04.00Hz-08.00Hz_fit.pkl'
    with open(fname, 'rb') as f:
        tup = pickle.load(f)
    plot_fits(*tup, title='western event', mag=3.2)
    plt.annotate('a)', **annotate_kw)
    plt.savefig('../tmp/fits1.pdf', bbox_inches='tight')
    fname = '../tmp/plot_fits/4808371_04.00Hz-08.00Hz_fit.pkl'
    with open(fname, 'rb') as f:
        tup = pickle.load(f)
    plot_fits(*tup, title='eastern event', mag=2.6)
    plt.annotate('b)', **annotate_kw)
    plt.savefig('../tmp/fits2.pdf', bbox_inches='tight')
    # join PDF figures with (run command inside processing folder)
    # pdfjam --papersize '{142mm,178mm}' --nup '1x2' ../tmp/fits?.pdf -o ../results/figs_paper/fits.pdf