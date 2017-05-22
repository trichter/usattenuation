#!/usr/bin/env python
# Copyright 2017 Tom Eulenfeld, MIT license
"""Download waveforms and calculate results with Qopen package"""
from collections import defaultdict, OrderedDict
from copy import deepcopy
from glob import glob
import json
import logging

import matplotlib
matplotlib.use('AGG')

import numpy as np
from obspy import read_events, read_inventory
from obspy.clients.fdsn import Client
import shove
from qopen.core import ConfigJSONDecoder, get_freqs, run, sort_dict
from qopen.imaging import calc_dependent
from qopen.site import align_site_responses
from qopen.util import gerr, linear_fit

from util.cache import waveform_cache

OUT = '../results/'
TMP = '../tmp/'

EVENT_FNAME = '../data/events/individual/%s.xml'
EVENT_GLOB = EVENT_FNAME % '*'
STATIONS_FNAME = '../data/stations/us_array_bh.xml'
CACHE_PATH = '../data/waveforms/%s_%s_%ds_%ds.mseed'
CACHE_REQUEST_AGAIN = True
SHOVE_FNAME = 'lite://' + TMP + 'results.sqlite'
EVENTSJSON = TMP + '../data/events/us_events.json'
CONF_FNAME = 'conf_us.json'
CACHE_WIN = (-50, 200)

ORDER = ['title', 'author', 'license', 'description',
         'freqs', 'freq_bands', 'stations', 'Qsc', 'Qi', 'R', 'nobs', 'v0',
         'mean', 'error', 'slope', 'intercept']

IGNORE_STATIONS = None


class JsonShove(shove.Shove):

    def __getitem__(self, key):
        value = super(JsonShove, self).__getitem__(key)
        return json.loads(value)

    def __setitem__(self, key, value):
        value = json.dumps(value)
        super(JsonShove, self).__setitem__(key, value)

    def __enter__(self):
        return self

    def __exit__(self, type_, value, tb):
        self.close()


def calc_stuff(conf, eventids=None, sh=None):
    get_waveforms = waveform_cache(
        CACHE_PATH, CACHE_WIN, request_again=CACHE_REQUEST_AGAIN)(
                Client().get_waveforms)
    stations = read_inventory(STATIONS_FNAME, 'STATIONXML')
    print('warm up finished')
    if eventids is None:
        fnames = sorted(glob(EVENT_GLOB))
    else:
        fnames = [EVENT_FNAME % eid for eid in eventids]
    for fname in fnames:
        eventid = fname.split('/')[-1].replace('.xml', '')
        if sh is not None and eventid in sh:
            continue
        print('processing event %s...' % eventid)
        event = read_events(fname, 'QUAKEML')
        results = run(events=event, inventory=stations,
                      get_waveforms=get_waveforms, **conf)
        if sh is None:
            continue
        if results is None:
            sh[eventid] = None
            continue
        res = list(results['events'].values())[0]
        if res is None or all([g0 is None for g0 in res['g0']]):
            sh[eventid] = None
            continue
        for st in tuple(res['R'].keys()):
            if all([r is None for r in res['R'][st]]):
                res['R'].pop(st)
        sh[eventid] = res


def collect_stations(boreholes=True):
    borehole_stations = set()
    stations = read_inventory(STATIONS_FNAME, 'STATIONXML')
    coords = {}
    for net in stations.networks:
        for sta in net.stations:
            for cha in sta.channels:
                lat = cha.latitude or sta.latitude
                lon = cha.longitude or sta.longitude
                dep = cha.depth
                key = '%s.%s' % (net.code, sta.code)
                if dep > 0:
                    borehole_stations.add(key)
                if IGNORE_STATIONS and key in IGNORE_STATIONS:
                    continue
                coords[key] = (lat, lon)
    if boreholes:
        return borehole_stations
    else:
        return coords


BOREHOLE_STATIONS = collect_stations()


def eventresults2json(results):
    events = defaultdict(list)
    with open(EVENTSJSON) as f:
        events2 = json.load(f)
    for evid, val in sorted(results['events'].items(), key=lambda i: i[0]):
        if val is None:
            continue
        events['id'].append(evid)
        time, lat, lon, depth, mag, magtype = events2[evid]
        events['time'].append(time)
        events['lat'].append(lat)
        events['lon'].append(lon)
        events['MR'].append(mag)
        events['depth'].append(depth)
        for key in ('sds', 'Mw', 'M0', 'fit_error'):
            events[key].append(val[key] if key in val else None)
    order = ('id', 'time', 'lat', 'lon', 'depth', 'MR', 'Mw', 'M0',
             'fit_error', 'source_displacement_spectrum')
    # rename a key
    events['source_displacement_spectrum'] = events.pop('sds')
    with open(EVENTRESULTSJSON, 'w') as f:
        f.write(to_json(sort_dict(dict(events), order=order)))


GENERATE_TEST_FILE = False  # for tests of Qopen package
MAX_NOBS = 2


def collect_results(correct=True, events=False):
    with open(CONF_FNAME) as f:
        conf = json.load(f, cls=ConfigJSONDecoder)
    freq = get_freqs(**conf['freqs'])
    fbands = list(freq.values())
    freq = list(freq.keys())
    fa = np.asarray(freq)
    print('copy results from shove')
    r = {'freq': freq, 'events': {}, 'config': conf}
    if GENERATE_TEST_FILE:
        nobs_station = defaultdict(int)
    with JsonShove(SHOVE_FNAME, compress=True) as sh:
        msg = '%d in shove / %d events processed'
        print(msg % (len(sh), len(glob(EVENT_GLOB))))
        for evid in sh:
            if sh[evid] is None:
                continue
            if GENERATE_TEST_FILE:
                stations = sh[evid]['R'].keys()
                if min(nobs_station[s] for s in stations) >= MAX_NOBS:
                    continue
                for s in stations:
                    nobs_station[s] = nobs_station[s] + 1
            r['events'][evid] = sh[evid]
    if GENERATE_TEST_FILE:
        i = 6
        for evid in r['events']:
            for sta in r['events'][evid]['R']:
                r['events'][evid]['R'][sta] = \
                    r['events'][evid]['R'][sta][i:i + 1]
            r['events'][evid] = {'W': r['events'][evid]['W'][i:i + 1],
                                 'R': r['events'][evid]['R']}
        r['freq'] = r['freq'][i:i + 1]
        del r['config']
        with open(TMP + 'usarray_dataset.json', 'wb') as f:
            json.dump(r, f)
        return
    if IGNORE_STATIONS is not None:
        print('remove stations %s' % IGNORE_STATIONS)
        sh = r['events']
        for evid in sh.keys():
            for station in sh[evid]['R'].keys()[:]:
                if station in IGNORE_STATIONS:
                    sh[evid]['R'].pop(station)
                    if len(sh[evid]['R']) == 0:
                        sh.pop(evid)

    Nres = len(r['events'])
    print("number of events with results: %d" % Nres)
    rorig = r
    if correct:
        rorig = deepcopy(r)
        print('correct sensitivities')
        # r = align_site_responses(r, station=BOREHOLE_STATIONS, response=0.25)
        r = align_site_responses(r)
    if events:
        print('collect magnitudes')
        eventresults2json(r)
    print('collect results')
    data = defaultdict(lambda: defaultdict(lambda: [[] for _ in freq]))
    data['v0'] = defaultdict(list)
    results = defaultdict(lambda: defaultdict(dict))
    sh = r['events']
    sh2 = rorig['events']
    for evid in sh:
        for station in sh2[evid]['R']:
            v0 = sh2[evid]['v0']
            data['v0'][station].append(v0)
            for i in range(len(sh2[evid]['R'][station])):
                R1 = sh[evid]['R'][station][i]
                R2 = sh2[evid]['R'][station][i]
                b = sh[evid]['b'][i]
                g0 = sh[evid]['g0'][i]
                if R2 is not None:
                    if R1 is not None:
                        data['R'][station][i].append(R1)
                    Qi = calc_dependent('Qi', b, freq[i], v0)
                    data['Qi'][station][i].append(Qi)
                    Qsc = calc_dependent('Qsc', g0, freq[i], v0)
                    data['Qsc'][station][i].append(Qsc)
    print('calc mean, err and slope')
    for what in ('v0', 'R', 'Qsc', 'Qi'):
        for station in data[what]:
            d = data[what][station]
            if what == 'v0':
                results[what][station] = np.mean(d)
                continue
            mean, err, nobs = [], [], []
            for i in range(len(freq)):
                da = np.asarray(d[i], dtype=np.float)
                m, e1, e2 = gerr(da)
                e = np.log10(m) - np.log10(e1)
                m = np.log10(m)
                mean.append(None if not np.isscalar(m) or np.isnan(m) else m)
                err.append(None if not np.isscalar(e) or np.isnan(e) else e)
                nobs.append(np.count_nonzero(np.logical_not(np.isnan(da))))
            results[what][station]['mean'] = mean
            results[what][station]['error'] = err
            if what == 'Qsc':
                results['nobs'][station] = nobs
            mean = np.asarray(mean, dtype=np.float)
            ind = np.logical_not(np.isnan(mean))
            if what == 'R':
                continue
            elif np.count_nonzero(ind) < 3:
                slope, intercept = None, None
            else:
                slope, intercept = linear_fit(mean[ind], np.log10(fa[ind]))
            results[what][station]['slope'] = slope
            results[what][station]['intercept'] = intercept
    print('convert to OrderedDict and save json')
    for key in results:
        if key not in ('v0', 'nobs'):
            for sta in results[key]:
                results[key][sta] = sort_dict(results[key][sta], order=ORDER)
        results[key] = OrderedDict(sorted(results[key].items()))
    results['freqs'] = freq
    results['freq_bands'] = fbands
    results['stations'] = collect_stations(False)  # get coordinates
    results['title'] = TITLE
    results['description'] = DESC
    results['author'] = 'Tom Eulenfeld'
    # results['copyright'] = COPYR
    results['license'] = LICENSE
    results = sort_dict(results, order=ORDER)
    with open(RESULTSJSON, 'w') as f:
        f.write(to_json(results, nl_after_str=True))


# http://stackoverflow.com/a/10420059
INDENT = 4
SPACE = " "
NEWLINE = "\n"


def json_collection(o, tweak=None, tweak2=False, level=0, nl_after_str=False):
    if not isinstance(o, dict):
        tweak2 = not (len(o) > 0 and isinstance(o[0], str) and
                      nl_after_str)
    op = '{' if isinstance(o, dict) else '['
    cl = '}' if isinstance(o, dict) else ']'
    nl = "" if tweak2 else NEWLINE
    nl2 = " " if tweak2 else NEWLINE
    ind = 0 if tweak2 else INDENT
    ret = op + nl
    comma = ""
    for k in o:
        ret += comma
        comma = "," + nl2
        ret += SPACE * ind * (level + 1)
        if isinstance(o, dict):
            v = o[k]
            ret += '"' + str(k) + '":' + SPACE
            tweak2 = '.' in k
            tweak = ('%.4g' if (isinstance(v, float) and '.' in k)
                     else tweak)  # v0
            tweak = '%.4f' if isinstance(v, tuple) else tweak  # coordinates
            k = v
        ret += to_json(k, level + 1, tweak=tweak, tweak2=tweak2,
                       nl_after_str=nl_after_str)
    ret += nl + SPACE * ind * level + cl
    return ret


def to_json(o, level=0, tweak=None, tweak2=False, nl_after_str=False):
    ret = ""
    if isinstance(o, (dict, list, tuple)):
        ret += json_collection(o, tweak=tweak, tweak2=tweak2, level=level,
                               nl_after_str=nl_after_str)
    elif isinstance(o, str):
        ret += '"' + o + '"'
    elif isinstance(o, bool):
        ret += "true" if o else "false"
    elif isinstance(o, int):
        ret += str(o)
    elif isinstance(o, float) and tweak is None:
        ret += '%.3f' % o  # frequencies, values
    elif isinstance(o, float):
        ret += tweak % o
    elif o is None:
        ret += 'null'
    else:
        msg = "Unknown type '%s' for json serialization" % str(type(o))
        raise TypeError(msg)
    return ret


TITLE = 'Results of Qopen script for US array'
LICENSE = ('Creative Commons Attribution 4.0 International License (CC BY 4.0)'
           ', see http://creativecommons.org/licenses/by/4.0/')
DESC = [
    "It follow a short description of the available keys in the results "
    "dictionary.",
    "freqs:       central frequencies",
    "freq_bands: corner frequencies",
    "stations:   all used stations of US array with coordinates "
    "(latitude, longitude)",
    "Qi:         quality factor for intrinsic attenuation in the form "
    "-log_10(Q_i) for each station "
    "(fields 'mean', 'error', 'slope', 'intercept')",
    "Qsc:        quality factor for scattering attenuation in the form "
    "-log_10(Q_sc) for each station "
    "(fields 'mean', 'error', 'slope', 'intercept')",
    "R:          energy site amplification factors in the form log_10(R) for "
    "each station (fields 'mean', 'error')",
    "nobs:       number of observation for each station",
    "v0:         mean S-wave velocity determined from picks for each station",
    "The fields 'mean' (robust geometric mean) and 'error' "
    "(median absolute deviation) of Qi, Qsc, R contain one value for each "
    "frequency band. Same is the case for nobs.",
    "The fields 'slope' and 'intercept' (at 1 Hz) of Qi and Qsc describe the "
    "frequency dependence."
]

with open(CONF_FNAME) as f:
    CONF = json.load(f, cls=ConfigJSONDecoder)
RESULTSJSON = OUT + 'results.json'
EVENTRESULTSJSON = OUT + 'eventresults.json'


if __name__ == '__main__':
    with JsonShove(SHOVE_FNAME, compress=True) as sh:
        calc_stuff(CONF, sh=sh)
    logging.basicConfig(level=logging.DEBUG)
    collect_results(correct=True, events=True)
    RESULTSJSON = OUT + 'results_uncorrected.json'
    EVENTRESULTSJSON = OUT + 'eventresults_uncorrected.json'
    collect_results(correct=False, events=True)
    msg = '%d borehole stations: %s'
    print(msg % (len(BOREHOLE_STATIONS), ' '.join(sorted(BOREHOLE_STATIONS))))