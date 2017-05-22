#!/usr/bin/env python
# Copyright 2017 Tom Eulenfeld, MIT license
"""
Compare estimated moment magnitudes with reference catalogs

You will need to create files with the magnitudes from reference catalogs
(see comments in this script or ask me)
"""
from copy import deepcopy
import json
import numpy as np
from obspy import UTCDateTime as UTC
import re


with open('../results/eventresults.json') as f:
    ev = json.load(f)
with open('../results/eventresults_uncorrected.json') as f:
    ev2 = json.load(f)
ev = [dict(zip(ev, t)) for t in zip(*ev.values())]
ev2 = [dict(zip(ev2, t)) for t in zip(*ev2.values())]


def dic(a, des, MR):
    a = deepcopy(a)
    a['time'] = a['time'][:19].replace('T', ' ')
    a['des'] = des
    a['MR'] = MR
    a['depth'] = a['depth'] / 1000
    for b in ev2:
        if a['id'] == b['id']:
            a['Mw_uncorrected'] = b['Mw']
    return a


def iter_matching_event(lat, lon, ddeg, time):
    for a in ev:
        if (abs(a['lat'] - lat) < ddeg and abs(a['lon'] - lon) < ddeg and
                abs(time - UTC(a['time'])) < 5):
            yield a


# This prints information for the two events for which envelope fits are
# plotted in us_plot_fits.py script
for a in ev:
    if a['id'] in ('2511938', '4808371'):
        print('envelope fits for event', a, '\n')

# http://www.ncedc.org/ncedc/catalog-search.html
# -> Mechanism Catalog
# -> Moment Tensors in CMT Standard Format
with open('../data/events/refcatalogs/berkeley.txt') as f:
    text1 = f.read()
# http://www.eas.slu.edu/eqc/eqc_mt/MECH.NA/
# copy & paste from links
with open('../data/events/refcatalogs/slu.txt') as f:
    text2 = f.read()
# https://service.iris.edu/fdsnws/event/1/query?starttime=2005-01-01&
# endtime=2016-01-01&maxmag=4&magtype=MW&catalog=ISC&orderby=time&
# format=text&maxlat=51&minlon=-132&maxlon=-56&minlat=23&nodata=404
with open('../data/events/refcatalogs/iris_fdsn.txt') as f:
    text3 = f.read()

events = []

# Event Id= 21437727
# Date: 2005/02/05   Time: 18:43:30.3100 GMT
# Lat=  37.40  Lon= -121.48  Depth=   7.20
# Moment Tensor Depth=  11.00
# Expo= +22   -0.193150  -2.045900   2.239000   0.137150   0.253480   0.837600
# Md=  4.25  Ml=  4.42  Mw=  4.18  Scalar Moment= 2.317000E+22
# Fault plane:  strike= 56  dip=88  slip=   6
# Fault plane:  strike=326  dip=84  slip= 178
pattern = (r'Date:\s*([\d/]+)\s*Time:\s*([\d:.]+).*?'
           r'Lat=\s*([.\d-]+)\s*Lon=\s*([.\d-]+).+?Mw=\s*([\d.]+)')

for match in re.findall(pattern, text1, flags=re.DOTALL):
    date, time, lat, lon, mag = match
    lat = float(lat)
    lon = float(lon)
    mag = float(mag)
    utc = UTC(date.replace('/', '_') + ' ' + time)
    for a in iter_matching_event(lat, lon, 0.5, utc):
        v = (a['time'], 0, dic(a, 'Berkeley', mag))
        events.append(v)


# 2005/01/28 22:37:07 34.71 -111.00 3 4.0 Arizona Y 20050128223707 Mechanism
for line in text2.split('\n'):
    if line.startswith('#'):
        continue
    date, time, lat, lon, dep, mag, *bla = line.split()
    utc = UTC(date.replace('/', '_') + ' ' + time)
    lat = float(lat)
    try:
        lon = float(lon)
    except:
        lon = float(lon + dep)
        mag = bla[0]
    try:
        mag = float(mag)
    except:
        mag = float(mag[:-1])
    for a in iter_matching_event(lat, lon, 0.5, utc):
        v = (a['time'], 1, dic(a, 'SLU', mag))
        events.append(v)

# EventID | Time | Latitude | Longitude | Depth/km | Author | Catalog |
# Contributor | ContributorID | MagType | Magnitude | MagAuthor |
# EventLocationName
for line in text3.split('\n'):
    if line.startswith('#'):
        continue
    id_, time, lat, lon, dep, _, cat, _, _, mt, mag, mc, _ = line.split('|')
    if 'mw' not in mt.lower():
        continue
    utc = UTC(time)
    lat = float(lat)
    lon = float(lon)
    mag = float(mag)
    for a in iter_matching_event(lat, lon, 0.5, utc):
        if mc == '':
            print('ignore:', line, '\n')
            continue
        v = (a['time'], 2, dic(a, cat + ' ' + mc, mag))
        events.append(v)

print('Table for publication')
last_id = None
for event in sorted(events):
    s = (r'{id} & {time} & {lat:.3f} & {lon:.3f} & {depth:.1f} & '
         r'{des} & {MR} & {Mw_uncorrected:.2f} & {Mw:.2f}\\')
    # ignore non-unique events
    if last_id is not None and event[0] == last_id:
        continue
    last_id = event[0]
    print(s.format(**event[2]))

print('number of earthuqakes', len(events), '\n')


def residual(s, dif=0):
    mags = [(a[2][s] + dif, a[2]['MR']) for a in events]
    mags = np.array(mags)
    print(s, 'shifted by', dif)
    print('mean signed deviation', np.mean(mags[:, 0] - mags[:, 1]))
    print('   mean abs deviation', np.mean(np.abs(mags[:, 0] - mags[:, 1])))
    print('                 stdv',
          (np.mean((mags[:, 0] - mags[:, 1]) ** 2)) ** 0.5)


residual('Mw_uncorrected')
residual('Mw_uncorrected', -0.077)
residual('Mw')
residual('Mw', -0.2808)
