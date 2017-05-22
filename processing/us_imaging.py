#!/usr/bin/env python
# Copyright 2017 Tom Eulenfeld, MIT license
"""Interpolate results, plot maps and other figures"""
from collections import defaultdict
import json

import cartopy.feature as cfeature
import matplotlib as mpl
from matplotlib.legend_handler import HandlerTuple
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.interpolate
from qopen.util import linear_fit

from util.geometry import points_outside_network
from util.projection import AEA, PC, add_scale, load_tif, write_tif


OUT = '../results/'
TMP = '../tmp/'
DATA_QVIEWER = '../visualization/Qviewer/data/'
TIFFOUT = OUT + 'tif/%s.tif'
PNGOUT = DATA_QVIEWER + '%s.png'
FIGOUT1 = OUT + 'figs_suppl/'
FIGOUT2 = OUT + 'figs_paper/'

EVENTRESULTSJSON = OUT + 'eventresults.json'
RESULTSJSON = OUT + 'results.json'
RESULTSJSON_UNCORRECTED = OUT + 'results_uncorrected.json'
CONDJSON = TMP + 'temp_nobscond.json'
CONF_FNAME = '../processing/conf_us.json'

OBS = ('nobs', 'Qsc', 'Qi', 'R', 'Qsc_error', 'Qi_error', 'R_error')
OBSSLOPE = ('Qsc_slope', 'Qi_slope')
OBSDERIV = ('B0', 'Qtot')
BOREHOLE_STATIONS = (
    'BK.BDM BK.CMB BK.CVS BK.GASB BK.HOPS BK.HUMO BK.KCC BK.MCCM BK.MNRC '
    'BK.MOD BK.ORV BK.PACP BK.PKD BK.SAO BK.WDC BK.WENL BK.YBH '
    'CI.MUR CI.SNCC II.PFO '
    'IU.ANMO IU.CCM IU.DWPF IU.HKT IU.RSSD IU.SSPA IU.TUC IU.WCI '
    'TA.Y22E US.ECSD US.GOGA').split()

USE_TIF = True
DX = 10000
XLIM = (-2500000, 2400000)
YLIM = (-1500000, 1700000)
GEOTRANSFORM = [XLIM[0], DX, 0, YLIM[1], 0, -DX]
GRID = tuple(np.meshgrid(np.arange(XLIM[0] + DX / 2, XLIM[1], DX),
                         np.arange(YLIM[1] - DX / 2, YLIM[0], -DX)))

ALPHA = 1 / 150e3
BUFFER = 50e3

TODROP = None
COND_INTERPOLATION = {'nobs': lambda x: x < 4}

CMAP = plt.get_cmap(lut=11)  # use matplotlib >= 2 for viridis colormap
DPI1 = 200
DPI2 = 100
ASPECT = 1.0 * (XLIM[1] - XLIM[0]) / (YLIM[1] - YLIM[0])
AXSIZE = (80 / 25.4, 80 / 25.4 / ASPECT)

with open(RESULTSJSON) as f:
    RES = json.load(f)
FREQ = np.array(RES['freqs'])
IGNORE_STATIONS = False


TICKS_Q = (-1.7, -1.9, -2.1, -2.3, -2.5, -2.6, -2.7, -2.8)
TICKS_QTOT = (-1.3, -1.5, -1.7, -1.9, -2.1, -2.2, -2.3, -2.4)

TICKS = {'Qi': [[t - 1, t - 0.5, t] for t in TICKS_Q],
         'Qsc': [[t - 1, t - 0.5, t] for t in TICKS_Q],
         'Qtot': [[t - 0.8, t - 0.5, t - 0.2] for t in TICKS_QTOT],
         'B0': [[0, 0.5, 1] for t in TICKS_Q],
         'Qi_slope': [-1.8, -1.0, -0.2],
         'Qsc_slope': [-1.8, -1.0, -0.2],
         'R': [[-2, 0, 2] for i in FREQ],
         'Qi_error': [[0, 0.5, 1.0] for i in FREQ],
         'Qsc_error': [[0, 0.5, 1.0] for i in FREQ],
         'R_error': [[0, 0.5, 1.0] for i in FREQ],
         #         'Qi_error2': [[0, 0.5, 1.0] for i in FREQ],
         #         'Qsc_error2': [[0, 0.5, 1.0] for i in FREQ],
         #         'R_error2': [[0, 0.5, 1.0] for i in FREQ],
         'nobs': [[0, 1, 2, 3] for i in FREQ],
         }

TITLE = {
    'nobs': r'$\log_{10}(n_{\mathrm{obs}})$',
    'Qi': r'$\log_{10}(Q_{\mathrm{i}}^{-1})$',
    'Qsc': r'$\log_{10}(Q_{\mathrm{sc}}^{-1})$',
    'Qtot': r'$\log_{10}(Q_{\mathrm{tot}}^{-1})$',
    'B0': r'$B_0$',
    'R': r'$\log_{10}R$',
    'Qi_error': r'$\mathrm{error}$ $\mathrm{of}$ $\log_{10}(Q_{\mathrm{i}}^{-1})$',
    'Qsc_error': r'$\mathrm{error}$ $\mathrm{of}$ $\log_{10}(Q_{\mathrm{sc}}^{-1})$',
    'R_error': r'error of $\log_{10}R$'}

LABEL = {'nobs': '$%s\,\mathrm{Hz}$',
         'Qi': '$%s\,\mathrm{Hz}$',
         'Qsc': r'$%s\,\mathrm{Hz}$',
         'Qtot': r'$%s\,\mathrm{Hz}$',
         'B0': r'$%s\,\mathrm{Hz}$',
         'Qi_error': r'$%s\,\mathrm{Hz}$',
         'Qsc_error': r'$%s\,\mathrm{Hz}$',
         'Qi_slope': r'$\mathrm{slope}$ $\gamma$',
         'Qsc_slope': r'$\mathrm{slope}$ $\gamma$',
         'R': r'$%s\,\mathrm{Hz}$',
         'R_error': r'$%s\,\mathrm{Hz}$'
         }
DESC = {
    'nobs': 'number of observations',
    'Qi': 'inverse of Q_i (quality factor of intrinsic attenuation)',
    'Qsc': 'inverse of Q_sc (quality factor of scattering attenuation)',
    'Qi_error': 'error of inverse of Q_i (quality factor of intrinsic attenuation)',
    'Qsc_error': 'error of inverse of Q_sc (quality factor of scattering attenuation)',
    'Qi_slope': 'slope of frequency dependency of -log_10 of Q_i (quality factor of intrinsic attenuation)',
    'Qsc_slope': 'slope of frequency dependency of -log_10 of Q_sc (quality factor of scattering attenuation)',
    'R': 'site amplification factor R',
    'R_error': 'error of site amplification factor R'}


def Qtot(Q1, Q2):
    return np.log10(10**Q1 + 10**Q2)


def B0(Q1, Q2):
    return 10**Q2 / (10**Q1 + 10**Q2)


def drop_stations(reset=False):
    global TODROP
    TODROP = {i: set() for i in range(len(FREQ))}
    TODROP[None] = set()
    if reset:
        return
    if IGNORE_STATIONS:
        for sta in IGNORE_STATIONS:
            TODROP[None].add(sta)
    num = defaultdict(int)
    for freqi in range(len(FREQ)):
        for obs in ('nobs', 'R', 'Qsc', 'Qi'):
            if obs not in COND_INTERPOLATION:
                continue
            for sta, val in RES[obs].items():
                if obs == 'nobs':
                    m = RES[obs][sta][freqi]
                else:
                    m = RES[obs][sta]['error'][freqi]
                if obs == 'R' and RES[obs][sta]['mean'][freqi] is None:
                    continue
                if ((m is None or COND_INTERPOLATION[obs](m)) and
                        sta not in TODROP[freqi]):
                    TODROP[freqi].add(sta)
                    num[sta] += 1
    for sta in num:
        if len(FREQ) - num[sta] < 4:
            TODROP[None].add(sta)


def collect_data(obs, freqi, return_droped=False):
    if obs in ('Qtot', 'B0'):
        r1 = collect_data('Qi', freqi, return_droped=False)
        r2 = collect_data('Qsc', freqi, return_droped=False)
        f = Qtot if obs == 'Qtot' else B0
        totdata = list(map(f, r1[2], r2[2]))
        return (r1[0], r1[1], totdata) + (None,) * 3 * return_droped
    lats, lons, data = [], [], []
    latsd, lonsd, datad = [], [], []
    if '_' in obs:
        obs, obs2 = obs.split('_')
    else:
        obs2 = 'mean'
    for key, val in RES[obs].items():
        if obs == 'R' and key in BOREHOLE_STATIONS:
            continue
        if obs == 'nobs':
            m = val[freqi]
            if m == 0:
                continue
            m = np.log10(m)
        else:
            m = val[obs2]
            if obs2 in ('mean', 'error'):
                m = m[freqi]
        if not (m is None or np.isnan(m)):
            lat, lon = RES['stations'][key]
            if key in TODROP[freqi] and obs != 'nobs':
                if not return_droped:
                    continue
                latsd.append(lat)
                lonsd.append(lon)
                datad.append(m)
            else:
                lats.append(lat)
                lons.append(lon)
                data.append(m)
    if return_droped:
        return lons, lats, data, lonsd, latsd, datad
    else:
        return lons, lats, data


def interp(lons, lats, data):
    xy = AEA.transform_points(PC, np.array(lons), np.array(lats))
    x = xy[:, 0]
    y = xy[:, 1]
    args = ((x, y), np.array(data), GRID)
    zi = scipy.interpolate.griddata(*args, method='linear')
    zi2 = scipy.interpolate.griddata(*args, method='nearest')
    zi[np.isnan(zi)] = zi2[np.isnan(zi)]
    outside = points_outside_network(x, y, *GRID, alpha=ALPHA, buf=BUFFER)
    zi[outside] = np.nan
    return zi


def plot1ax(obs, freqi=0, vmin=None, vmax=None, new_map=True,
            interpolate=True,
            cmap=CMAP, station_plot=False, png_plot=False, plotdrop=False):
    ax = plt.gca()
    if new_map:
        basemap(ax=ax)
    try:
        ticks = TICKS[obs]
    except KeyError:
        ticks, vmin, vmax = None, None, None
    if freqi is not None and ticks is not None:
        ticks = ticks[freqi]
    if ticks is not None:
        vmin, vmax = ticks[0], ticks[-1]

    if interpolate:
        if USE_TIF:
            if obs in ('Qtot', 'B0'):
                fname = get_tifname('Qi', freqi)
                zi1, extent = load_tif(fname)
                fname = get_tifname('Qsc', freqi)
                zi2, extent2 = load_tif(fname)
                f = Qtot if obs == 'Qtot' else B0
                zi = f(zi1, zi2)
            else:
                fname = get_tifname(obs, freqi)
                zi, extent = load_tif(fname)
            assert extent == XLIM + YLIM
        else:
            zi = interp(*collect_data(obs, freqi))
            extent = XLIM + YLIM
        sc = ax.imshow(zi, vmin=vmin, vmax=vmax, extent=extent, cmap=cmap,
                       rasterized=True, origin='upper')
    else:
        lons, lats, data, lonsd, latsd, datad = \
            collect_data(obs, freqi, return_droped=True)
        if plotdrop:
            ax.scatter(lonsd, latsd, s=8, c=datad, cmap=cmap, vmin=vmin,
                       vmax=vmax, linewidths=(0,), zorder=10, transform=PC)
            ax.scatter(lonsd, latsd, edgecolor='#CC6600', facecolor='none',
                       marker='o', s=8, linewidths=(0.4,), zorder=11,
                       transform=PC)
        sc = ax.scatter(lons, lats, s=8, c=data, cmap=cmap, vmin=vmin,
                        vmax=vmax, linewidths=(0,), zorder=10, transform=PC)
    # plot colorbar
    box = ax.get_position()
    dw, dh = box.x1 - box.x0, box.y1 - box.y0
    kw = {}
    if '_error' in obs:
        kw['extend'] = 'max'
    if obs == 'nobs':
        kw['format'] = '%d'
    if station_plot:
        rect = [box.x1 - 0.1 * dw, box.y0 + 0.05 * dh, 0.012 * dw, 0.3 * dh]
        kw['extend'] = 'max'
        mpl.rcParams['ytick.major.size'] = 3
        mpl.rcParams['axes.linewidth'] = 1
    else:
        rect = [box.x1 - 0.14 * dw, box.y0 + 0.07 * dh, 0.015 * dw, 0.3 * dh]
        mpl.rcParams['ytick.major.size'] = 2
        mpl.rcParams['axes.linewidth'] = 0.5
        mpl.rcParams['ytick.labelsize'] = 'x-small'
    cax = plt.axes(rect)
    kw['ticks'] = ticks
    kw['cax'] = cax
    plt.colorbar(sc, **kw)
    mpl.rcdefaults()
    # plot label
    if station_plot:
        label = TITLE[obs]
        ax.annotate(label, (0.99, 0.35), xycoords='axes fraction',
                    ha='right', va='bottom')
    else:
        label = LABEL[obs]
        if png_plot and obs in TITLE:
            label = TITLE[obs] + '\n' + r'$\mathrm{at}$ ' + label
            if 'error' in obs:
                obs2 = obs.strip('_error')
                label = (TITLE[obs2] + '\n' +
                         r'$\mathrm{error}$ $\mathrm{at}$ ' + LABEL[obs2])
        if png_plot and 'slope' in obs:
            label = (label + '\n' + r'$\mathrm{of}$ ' +
                     TITLE[obs.strip('_slope')])
        if '%' in label:
            label = label % FREQ[freqi]
        if png_plot == 2:  # fast hack
            return sc
        label = label.replace('.', '.\!') # hack for wrong spacing in mpl
        ax.annotate(label, (0.04, 0.02), xycoords='axes fraction',
                    size='small', va='bottom')
    return sc


def basemap(ax=None, label_scale=False):
    if ax is None:
        ax = plt.gcf().add_axes([0, 0, 1, 1], projection=AEA)
    ax.set_extent(XLIM + YLIM, crs=AEA)
    ax.background_patch.set_color('0.75')
    states = cfeature.NaturalEarthFeature(
        category='cultural', name='admin_1_states_provinces_lakes',
        scale='110m', facecolor='none', edgecolor='k')
    ax.add_feature(states, linewidth=0.5)
    add_scale(ax, 500, (-87, 27), label=label_scale, size='large')
    return ax


def plot_stations():
    with open(EVENTRESULTSJSON) as f:
        events = json.load(f)
    fig = plt.figure(figsize=(8, 5))
    ax = basemap(label_scale=True)
    lat, lon = zip(*RES['stations'].values())
    ax.plot(lon, lat, 'v', ms=4, mfc='none', marker='v', mec='k', mew=0.5,
            zorder=100, transform=PC)
    plot1ax('nobs', freqi=2, cmap=plt.get_cmap('hot_r', lut=11),
            station_plot=True, new_map=False)
    ax.plot(events['lon'], events['lat'], '.', ms=1, color='#377eb8',
            zorder=50, transform=PC)  # '#377eb8'
    ax.plot((-117.4748, -82.149), (34.1322, 37.7211), 'x', ms=5, color='white',
            transform=PC, zorder=101)
    fig.savefig(FIGOUT2 + 'map_stations.pdf', bbox_inches='tight', dpi=DPI2)


def plot_grid(obsgrid, label=False, title=False, fig=None, plotdrop=False):
    Ni, Nj = (len(obsgrid[0]), len(obsgrid))
    figsize = (Ni * AXSIZE[0] + 1, Nj * AXSIZE[1] + 1)
    w = AXSIZE[0] / figsize[0]
    h = AXSIZE[1] / figsize[1]
    wpad = 0.1 / figsize[0]
    hpad = 0.1 / figsize[1]
    if fig is None:
        fig = plt.figure(figsize=figsize)
    for i in range(Ni):
        for j in range(Nj):
            try:
                obs = obsgrid[j][i]
            except:
                continue
            if obs is None:
                continue
            rect = [wpad + i * (wpad + w),
                    hpad + (Nj - j - 1) * (hpad + h),
                    w,
                    h]
            fig.add_axes(rect, projection=AEA)
            if title:
                plt.annotate(TITLE[obs[0]], (1.05, 1.05),
                             xycoords='axes fraction',
                             va='bottom', ha='center')
                title = False
            plot1ax(obs[0], obs[1], interpolate=obs[2], png_plot=label,
                    plotdrop=plotdrop)


def image_suppl(obs):
    obsgrid = [((obs, 0, False), (obs, 0, True)),
               ((obs, 2, False), (obs, 2, True)),
               ((obs, 4, False), (obs, 4, True)),
               ((obs, 6, False), (obs, 6, True))]
    if obs in ('Qsc', 'Qi'):
        obs2 = obs + '_slope'
        obsgrid.append(((obs2, None, False), (obs2, None, True)))
    plot_grid(obsgrid, label=False, title=True)
    plt.savefig(FIGOUT1 + obs + '.pdf', dpi=DPI1, bbox_inches='tight')
    plt.close()

def image_paper():
    obsgrid = [(('Qsc', 0, True), ('Qi', 0, True)),
               (('Qsc', 2, True), ('Qi', 2, True)),
               (('Qsc', 4, True), ('Qi', 4, True)),
               (('Qsc', 6, True), ('Qi', 6, True))]
    plot_grid(obsgrid, label=True)
    plt.savefig(FIGOUT2 + 'Q.pdf', dpi=DPI1, bbox_inches='tight')

    obsgrid = [(('Qsc_slope', None, True), ('Qi_slope', None, True))]
    plot_grid(obsgrid, label=True)
    plt.savefig(FIGOUT2 + 'Q_slope.pdf', dpi=DPI1, bbox_inches='tight')

    obsgrid = [(('R', 0, True), ('R', 4, True)),
               (('R', 2, True), ('R', 6, True))]
    plot_grid(obsgrid, label=True)
    plt.savefig(FIGOUT2 + 'R.pdf', dpi=DPI1, bbox_inches='tight')

    obsgrid = [(('Qsc', 4, False), ('Qsc', 4, True))]
    plot_grid(obsgrid, label=True, plotdrop=True)
    plt.savefig(FIGOUT2 + 'Qsc_6Hz.pdf', dpi=DPI1, bbox_inches='tight')

    obsgrid = [(('Qtot', 0, True), ('Qtot', 6, True))]
    plot_grid(obsgrid, label=True)
    plt.savefig(FIGOUT2 + 'Qtot.pdf', dpi=DPI1, bbox_inches='tight')

    drop_stations(reset=True)
    obsgrid = [(None, ('R', 4, False))]
    plot_grid(obsgrid, label=True)
    global RES
    with open(RESULTSJSON_UNCORRECTED) as f:
        RES = json.load(f)
    obsgrid = [(('R', 4, False), None)]
    plot_grid(obsgrid, label=True, fig=plt.gcf())
    plt.savefig(FIGOUT2 + 'R_6Hz.pdf', dpi=DPI1, bbox_inches='tight')
    plt.close('all')


def get_tifname(obs, freqi=None):
    s = obs if freqi is None else '%s_%04.1fHz' % (obs, FREQ[freqi])
    return TIFFOUT % s


def calc_tifs():
    for obs in OBS + OBSSLOPE:
        if '_slope' in obs:
            zi = interp(*collect_data(obs, None))
            desc = DESC[obs]
            fname = get_tifname(obs)
            write_tif(fname, zi, GEOTRANSFORM, DESC[obs])
        else:
            for freqi in range(len(FREQ)):
                zi = interp(*collect_data(obs, freqi))
                desc = 'log_10 of %s at %sHz' % (DESC[obs], FREQ[freqi])
                fname = get_tifname(obs, freqi)
                write_tif(fname, zi, GEOTRANSFORM, desc)


def plot_usstates():
    fig = plt.figure(figsize=AXSIZE)
    ax = basemap()
    ax.outline_patch.set_visible(False)
    fig.savefig(PNGOUT % 'us', dpi=250)
    plt.close(fig)


def plot_png(obs, freqi, interpolate=True):
    fig = plt.figure(figsize=AXSIZE)
    ax = basemap()
    ax.outline_patch.set_visible(False)
    plot1ax(obs, freqi, interpolate=interpolate, png_plot=True)
    label = obs
    if freqi is not None:
        label = label + '_%04.1fHz' % FREQ[freqi]
    if interpolate:
        label = label + '_int'
    fig.savefig(PNGOUT % label, dpi=250)
    plt.close(fig)


def plot_all_pngs():
    plot_usstates()
    for obs in OBS + OBSSLOPE + OBSDERIV:
        for freqi in range(len(FREQ)):
            if '_slope' in obs and freqi != 0:
                continue
            elif '_slope' in obs:
                freqi = None
            for interpolate in (False, True):
                plot_png(obs, freqi, interpolate=interpolate)


def write_js():
    stations = list(RES['stations'].keys())
    normcoords = []
    for coord in RES['stations'].values():
        lat, lon = coord
        x, y = AEA.transform_point(lon, lat, PC)
        x = 1. * (x - XLIM[0]) / (XLIM[1] - XLIM[0])
        y = 1. * (y - YLIM[0]) / (YLIM[1] - YLIM[0])
        normcoords.append((x, y))
    with open(DATA_QVIEWER + '/stations.js', 'w') as f:
        f.write('var stations = ')
        json.dump([stations, normcoords], f)
        f.write(';')
    with open(DATA_QVIEWER + '/results.js', 'w') as f:
        f.write('var results = ')
        json.dump(RES, f)
        f.write(';')


def plot_mean():
    import seaborn as sns
    sns.set_style('ticks')
    d = [[] for _ in range(9)]
    for i in range(len(FREQ)):
        lons1, lats1, data1 = collect_data('Qi', i)
        lons2, lats2, data2 = collect_data('Qsc', i)
        data3 = Qtot(np.array(data1), np.array(data2))
        d[0].append(np.mean(data1))
        d[1].append(np.mean(data2))
        d[2].append(np.mean(data3))
        ind = np.array(lons1) < -103
        d[3].append(np.mean(np.array(data1)[ind]))
        d[4].append(np.mean(np.array(data1)[~ind]))
        d[5].append(np.mean(np.array(data2)[ind]))
        d[6].append(np.mean(np.array(data2)[~ind]))
        d[7].append(np.mean(np.array(data3)[ind]))
        d[8].append(np.mean(np.array(data3)[~ind]))
    fs = 100 / 25.4
    fig = plt.figure(figsize=(fs, 0.9 * fs))
    ax = fig.add_subplot(111)
    d = [10 ** np.array(dx) for dx in d]
    c1, c2 = '#e41a1c #275982'.split()
    kw = {'lw': 2, 'dash_capstyle': 'round'}
    ax.loglog(FREQ, d[2], color='#4E4E4E', label='total', **kw)
    kw['dashes'] = (3, 2)
    ax.loglog(FREQ, d[0], '--', color='#4E4E4E', label='intrinsic', **kw)
    ax.loglog(FREQ, d[3], '--', color='#CF0000', alpha=0.5, **kw)
    ax.loglog(FREQ, d[4], '--', color='#006AFF', alpha=0.5, **kw)
    kw['dashes'] = (2, 2, 0.1, 2)
    ax.loglog(FREQ, d[1], '-.', color='#4E4E4E', label='scattering', **kw)
    ax.loglog(FREQ, d[5], '-.', color='#CF0000', alpha=0.5, **kw)
    ax.loglog(FREQ, d[6], '-.', color='#006AFF', alpha=0.5, **kw)
    kw.pop('dashes')
    ax.loglog(FREQ, d[7], color='#CF0000', label='west', alpha=0.5, **kw)
    ax.loglog(FREQ, d[8], color='#006AFF', label='east', alpha=0.5, **kw)
    ax.legend(fontsize='small', labelspacing=0.3, handlelength=1.8,
              bbox_to_anchor=(1.1, 0.9))
    ax.set_xlim(1.5 / 1.2, 17 * 1.2)
    ax.set_ylim(4e-4, 2e-2)
    ax.set_xticks([1.5, 3, 6, 12])
    ax.set_xticklabels([1.5, 3., 6., 12.])
    ax.set_xticks([2.1, 4.2, 8.5, 17.], minor=True)
    ax.set_xlabel('frequency (Hz)')
    ax.set_ylabel('attenuation Q$^{-1}$')
    slope, intercept = linear_fit(np.log10(d[2]), np.log10(FREQ))
    l = r'$Q_{\mathrm{tot}}=%.0f\times f^{%.1f}$' % (10**(-intercept), -slope)
    ax.annotate(l, (1.7, 1.3e-2))
    sns.despine()
    fig.tight_layout()
    fig.savefig(FIGOUT2 + 'Qmean.pdf', bbox_inches='tight')
    plt.close(fig)


def _format_axes(ax1, ax2):
    import seaborn as sns
    ax1.set_xscale('log')
    ft = FREQ
    ax1.set_xlim(ft[0] / 1.2, ft[-1] * 1.2)
    ax1.set_xticks(ft[::2])
    ax1.set_xticks(ft[1::2], minor=True)
    ax1.set_xticklabels(ft[::2])
    ax1.set_xlabel('frequency (Hz)')
    ax2.set_xlabel('frequency (Hz)')
    ax1.set_ylabel(r'$\log_{10}Q^{-1}$')
    ax1.annotate('scattering', (0.05, 0.05), xycoords=('axes fraction'))
    ax2.annotate('intrinsic', (0.05, 0.05), xycoords=('axes fraction'))
    plt.setp(ax2.get_yticklabels(), visible=False)
    sns.despine(ax=ax1)
    sns.despine(ax=ax2)


def plot_comparison_CI_BFS(ax1, ax2):
    station = 'CI.BFS'
    Qsc = RES['Qsc'][station]['mean']
    Qscerr = RES['Qsc'][station]['error']
    Qi = RES['Qi'][station]['mean']
    Qierr = RES['Qi'][station]['error']
    # AA98 table 2
    Qfrome = lambda f, e: np.log10(3.5 * np.array(e) / (2 * np.pi * f))
    l1 = 'AA98 (CSP)'
    f1 = np.array([1.28, 2.02, 3.2, 5.1, 8.0, 12.56, 20])
    f1 = (f1[:-1] + f1[1:]) / 2
    e1sc = [0.049, 0.033, 0.013, 0.011, 0.013, 0.022]
    e1i = [0.049, 0.044, 0.031, 0.027, 0.028, 0.036]
    # AA98 table 1a
    l2 = 'AA98 (CJP)'
    f2 = np.array([8, 12.65, 20, 31.62])
    f2 = (f2[:-1] + f2[1:]) / 2
    e2sc = [0.0048, 0.0084, 0.013]
    e2i = [0.021, 0.031, 0.041]
    # LA94  table 2
    l3 = 'LA94 (CJP)'
    f3 = [7, 14, 28]
    Q3sc = -np.log10([60000, 28000, 85000])
    Q3i = -np.log10([610, 750, 1150])

    l = 'this study (%s)' % station
    ax1.errorbar(FREQ , Qsc, Qscerr, fmt='d', label=l, color='0.4', zorder=0)
    ax1.plot([], [])
    ax1.plot(f1, Qfrome(f1, e1sc), 'o', label=l1)
    ax1.plot(f2, Qfrome(f2, e2sc), '*', label=l2, zorder=3)
    ax1.plot(f3, Q3sc, 's', label=l3)
    ax2.errorbar(FREQ , Qi, Qierr, fmt='d', label=l, color='0.4', zorder=0)
    ax2.plot([], [])
    ax2.plot(f1, Qfrome(f1, e1i), 'o', label=l1)
    ax2.plot(f2, Qfrome(f2, e2i), '*', label=l2)
    ax2.plot(f3, Q3i, 's', label=l3)
    ax2.legend(loc='lower right', fontsize='x-small')


def _plot_own_California_Q(ax, key='Qsc'):
    Q = []
    for station, (lat, lon) in RES['stations'].items():
        if lat < 39 and lon < -115 and lat < -0.75 * lon - 51:
            try:
                Q.append(RES[key][station]['mean'])
            except KeyError:
                pass
    print('Number of used California stations: ', len(Q))
    Q = np.array(Q, dtype=np.float)
    p25 = np.nanpercentile(Q, 25, axis=0)
    p75 = np.nanpercentile(Q, 75, axis=0)
    l1, = ax.plot(FREQ, np.nanmedian(Q, axis=0), '0.2')
    l2, = ax.plot(FREQ, p25, color='0.8')
    ax.plot(FREQ, p75, color='0.8')
    l3, = ax.plot(FREQ, np.nanmin(Q, axis=0), '0.8', ls='--')
    ax.plot(FREQ, np.nanmax(Q, axis=0), '0.8', ls='--')
    return l3, l2, l1, l2, l3


# http://stackoverflow.com/a/40363560
class HandlerTupleVertical(HandlerTuple):
    def __init__(self, **kwargs):
        HandlerTuple.__init__(self, **kwargs)

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        numlines = len(orig_handle)
        handler_map = legend.get_legend_handler_map()
        height_y = (3 * height / numlines)
        leglines = []
        for i, handle in enumerate(orig_handle):
            handler = legend.get_legend_handler(handler_map, handle)
            legline = handler.create_artists(
                legend, handle, xdescent,
                height + (2 * i + 1 - numlines) * height_y,
                width, 2 * height, fontsize, trans)
            leglines.extend(legline)
        return leglines


def plot_comparison_California(ax1, ax2):
    # Mayeda et al. 1992 table 1
    f1 = [1.5, 3, 6, 9, 12, 15]
    Qsc_LV = [0.01916, 0.00903, 0.00179, 0.00081, 0.00062, 0.00064]
    Qi_LV = [0.01312, 0.00317, 0.00179, 0.00095, 0.00070, 0.00061]
    Qsc_CC = [0.01275, 0.00598, 0.00181, 0.00135, 0.00085, 0.00051]
    Qi_CC  = [0.00359, 0.00336, 0.00296, 0.00240, 0.00172, 0.00125]
    # Jin et al. 1994 table 1
    f2 = [0.75, 1.5, 3, 6, 12, 24]
    Qsc_SC = [
        [0.018, 0.006, 0.0011, 0.0007, None, None],  # GSC
        [0.011, 0.004, 0.0011, 0.0008, None, None],  # ISA
        [0.067, 0.008, 0.003, 0.0008, 0.0009, 0.0005],  # PAS
        [0.033, 0.007, 0.0024, 0.0011, None, None],  # PFO
        [0.021, 0.008, 0.002, 0.0006, 0.0003, 0.00008]]  # SVD
    Qi_SC = [
        [0.012, 0.010, 0.0038, 0.0019, None, None],  # GSC
        [0.013, 0.008, 0.0032, 0.0021, None, None],  # ISA
        [0.020, 0.008, 0.005, 0.0026, 0.0018, 0.0011],  # PAS
        [0.014, 0.008, 0.0049, 0.0029, None, None],  # PFO
        [0.010, 0.011, 0.0043, 0.0017, 0.0009, 0.0003]]  # SVD
    Qsc_SC = np.transpose(np.array(Qsc_SC, dtype=np.float))
    Qi_SC = np.transpose(np.array(Qi_SC, dtype=np.float))
    l1 = 'Long Valley'
    l2 = 'Central California'
    l3 = 'Southern California'

    _plot_own_California_Q(ax1, 'Qsc')
    ax1.plot(f1, np.log10(Qsc_LV), 'o', label=l1)
    # #AEB062 #C44E52
    ax1.plot(f1, np.log10(Qsc_CC), '*', mfc='#AEB062', label=l2, zorder=3)
    kw = dict( ls='', ms=5, marker='.', mfc='#C44E52', zorder=5)
    ax1.plot(f2, np.log10(Qsc_SC), **kw)
    h = _plot_own_California_Q(ax2, 'Qi')
    ax2.plot(f1, np.log10(Qi_LV), 'o', label=l1)
    ax2.plot(f1, np.log10(Qi_CC), '*', mfc='#AEB062', label=l2)
    ax2.plot(f2, np.log10(Qi_SC), **kw)
    ax2.plot([], [], label=l3, **kw)
    handles, labels = ax2.get_legend_handles_labels()
    handles.append(h)
    labels.append('this study (0, 25, 50,\n75, 100 percentiles)')
    ax2.legend(handles, labels, handler_map = {tuple : HandlerTupleVertical()},
               loc='lower right', fontsize='x-small', handlelength=1.6)


def plot_comparison_with_other_attenuation_studies():
    import seaborn as sns
    sns.set_style('ticks')
    fig = plt.figure(figsize=(168 / 25.4, 6.5))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222, sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(223, sharex=ax1)
    ax4 = fig.add_subplot(224, sharex=ax3, sharey=ax3)
    plot_comparison_CI_BFS(ax1, ax2)
    plot_comparison_California(ax3, ax4)
    ax3.set_ylim((-4.6, -1.4))
    _format_axes(ax1, ax2)
    _format_axes(ax3, ax4)
    kw = dict(xy=(0, 1), xytext=(-60, -10), xycoords='axes fraction',
              textcoords='offset points', fontsize='large',
              annotation_clip=False)
    ax1.annotate('a)', **kw)
    ax3.annotate('b)', **kw)
    plt.tight_layout(w_pad=1.5)
    fig.savefig(FIGOUT2 + 'comparison.pdf', bbox_inches='tight')
    plt.close(fig)


def plot_Mw_MR_comparison(MR=3.4):
    import seaborn as sns
    sns.set_style('ticks', {'axes.grid': True})
    with open(EVENTRESULTSJSON) as f:
        events = json.load(f)
    mag = np.asarray(events['MR'], dtype=np.float)
    Mw = np.asarray(events['Mw'], dtype=np.float)
    lon = np.asarray(events['lon'], dtype=np.float)
    sds = np.asarray(events['source_displacement_spectrum'], dtype=np.float)
    fig = plt.figure(figsize=(168 / 25.4, 4))
    rect1 = [0.1, 0.15, 0.4, 0.6]
    rect2 = [0.55, 0.15, 0.4, 0.6]
    rect3 = [0.1, 0.78, 0.4, 0.15]
    rect4 = [0.55, 0.78, 0.4, 0.15]
    ax1 = fig.add_axes(rect1)
    ax2 = fig.add_axes(rect2, sharex=ax1, sharey=ax1)
    ax3 = fig.add_axes(rect3)
    ax4 = fig.add_axes(rect4, sharex=ax3, sharey=ax3)
    ind3 = lon < -103  # west
    ind4 = lon > -103  # east
    kw = dict(bins=np.arange(1.45, 3.65, 0.1), alpha=0.75, rwidth=0.8)
    ax3.hist(mag[ind3][~np.isnan(Mw[ind3])], **kw)
    ax4.hist(mag[ind4][~np.isnan(Mw[ind4])], **kw)
    kw = dict(fliersize=3, linewidth=1)
    sns.boxplot(mag[ind3], Mw[ind3], ax=ax1, **kw)
    sns.boxplot(mag[ind4], Mw[ind4], ax=ax2, **kw)

    def cx(x):
        return 10 * (x - 1.5)

    x = np.array([1, 4])
    kw = dict(color='gray', alpha=0.5)
    ax1.plot(cx(x), 0.9 * x + 0.1, **kw)  # Mw SCSN
    ax1.plot(cx(x), 0.7 * x + 1.1, **kw)  # Mw NCSN
    ax1.plot(cx(x), 0.6 * x + 1.3, **kw)  # Mw UNR-NBE
    ax1.plot(cx(x), 0.7 * x + 0.8, **kw)  # Mw UUSS
    ax2.plot(cx(x), 0.6 * x + 0.8, **kw)  # Ml CERI
    ax2.plot(cx(x), 0.5 * x + 1.5, **kw)  # Ml LCSN
    kw = dict(textcoords='offset points', color='gray', alpha=0.8,
              size='x-small')
    ax1.annotate('4', (cx(1.5), 1.5), (0, 2), **kw)
    ax1.annotate('1', (cx(3.5), 0.7 * 3.5 + 1.1), (6, 0), **kw)
    ax1.annotate('2', (cx(3.5), 0.6 * 3.5 + 1.3), (6, 0), **kw)
    ax1.annotate('3', (cx(1.5), 0.7 * 1.5 + 0.8), (0, -7), **kw)
    ax2.annotate('6', (cx(3.5), 0.6 * 3.5 + 0.8), (6, 0), **kw)
    ax2.annotate('5', (cx(3.5), 0.7 * 3.5 + 0.8), (6, 0), **kw)
    t = '1 - Mw NCSN\n2 - Mw UNR-NBE\n3 - Mw UUSS\n4 - Mw SCSN'
    t2 = '5 - Ml LCSN\n6 - Ml CERI'
    kw.update(dict(ha='right', va='bottom'))
    ax1.annotate(t, (1, 0), (-5, 5), 'axes fraction', **kw)
    ax2.annotate(t2, (1, 0), (-5, 5), 'axes fraction', **kw)
    ax3.grid(False)
    ax4.grid(False)
    ax1.set_xticks(ax1.get_xticks()[::5])
    ax1.set_xticklabels(np.arange(1.5, 4, 0.5))
    ax1.set_ylim(1.5, 4.2)
    ax3.set_xlim(1.45, 3.55)
    ax1.set_ylabel('estimated Mw')
    ax1.set_xlabel('ANF MR')
    ax2.set_xlabel('ANF MR')
    ax3.set_yticks([0, 500, 1000, 1500])
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)
    sns.despine(ax=ax1)
    sns.despine(ax=ax2)
    sns.despine(ax=ax3)
    sns.despine(ax=ax4)
    kw = dict(xy=(1, 1), xytext=(-15, -15), xycoords='axes fraction',
              textcoords='offset points', ha='right')
    ax3.annotate('west', **kw)
    ax4.annotate('east', **kw)
    fig.savefig(FIGOUT2 + 'mags.pdf')
    fig2 = plt.figure(figsize=(5.5, 3))
    ax5 = fig2.add_subplot(121)
    ax6 = fig2.add_subplot(122, sharex=ax5, sharey=ax5)
    cond1 = np.logical_and(ind3, mag == MR)
    cond2 = np.logical_and(ind4, mag == MR)
    for i in range(np.count_nonzero(cond1)):
        ax5.loglog(FREQ, sds[cond1, :][i, :], 'k', alpha=0.1)
    for i in range(np.count_nonzero(cond2)):
        ax6.loglog(FREQ, sds[cond2, :][i, :], 'k', alpha=0.1)
    ax5.set_ylabel(r'$\omega$M (Nm)')
    ax5.set_xlabel('frequency (Hz)')
    ax6.set_xlabel('frequency (Hz)')
    ax5.set_xlim(1.5 / 1.2, 17 * 1.2)
    ax5.set_ylim(1e12 / 3, 1e15 * 3)
    ax5.set_xticks([1.5, 3, 6, 12])
    ax5.set_xticklabels([1.5, 3., 6., 12.])
    ax5.set_xticks([2.1, 4.2, 8.5, 17.], minor=True)
    plt.setp(ax6.get_yticklabels(), visible=False)
    ax5.tick_params(axis='y', which='minor', left='off', right='off')
    ax6.tick_params(axis='y', which='minor', left='off', right='off')
    ax5.annotate('west', **kw)
    ax6.annotate('east', **kw)
    plt.tight_layout(w_pad=1.5)
    sns.despine(ax=ax5)
    sns.despine(ax=ax6)
    fig2.savefig(FIGOUT2 + 'sds_MR%.1f.pdf' % MR, bbox_inches='tight')
    plt.close('all')

write_js()
plot_stations()
drop_stations()
print('Number of stations:', len(RES['stations']))
print('Number of stations with min 1 observation:', len(RES['nobs']))
print('Number of stations with less than 4 observations (freqi=2):',
      len(TODROP[2]))
print('Number of stations not used for computing the slope:',
      len(TODROP[None]))
calc_tifs()
plot_all_pngs()
for w in OBS + OBSDERIV:
    image_suppl(w)
image_paper()

plot_comparison_with_other_attenuation_studies()
plot_Mw_MR_comparison()
plot_mean()
