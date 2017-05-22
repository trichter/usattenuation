#!/usr/bin/env python
# Copyright 2017 Tom Eulenfeld, MIT license
"""Download station and event metadata from IRIS"""

from collections import OrderedDict
import json

from obspy import read_events
from obspy.core import UTCDateTime as UTC
from obspy.core.event import ResourceIdentifier
from obspy.clients.fdsn import Client as FSDNClient

EVENT_FNAME = '../data/events/us_events'
EVENT_PICKS_FNAME = '../data/events/individual/'
STATION_FNAME = '../data/stations/us_array'


def _extract_eventid(event):
    return str(event.resource_id).split('=')[-1]


def get_events(individual=True):
    endtime = UTC('2015-09-01')
    kw = {'starttime': UTC('2005-04-01'), 'endtime': endtime,
          'minmagnitude': 1.5, 'maxmagnitude': 3.5,
          'minlatitude': 25, 'maxlatitude': 49.5,
          'minlongitude': -124, 'maxlongitude': -67,
          'maxdepth': 40,
          'includearrivals': individual,
          'catalog': 'ANF'}
    client = FSDNClient()
    if individual:
        for t in range(6, 17):
            kw['endtime'] = UTC('20%02d-01-01' % t) if t < 16 else endtime
            events = client.get_events(**kw)
            for ev in events:
                id_ = _extract_eventid(ev)
                ev.resource_id = ResourceIdentifier(id_)
                ev.write(EVENT_PICKS_FNAME + id_ + '.xml', 'QUAKEML')
            print('fetched events of year %d' % kw['starttime'].year)
            kw['starttime'] = kw['endtime']
    else:
        events = client.get_events(**kw)
        events.plot(projection='local', outfile=EVENT_FNAME + '.png')
        for ev in events:
            id_ = _extract_eventid(ev)
            ev.resource_id = ResourceIdentifier(id_)
        events.write(EVENT_FNAME + '.xml', 'QUAKEML')
        return events


def merge_stations(s1, s2, verbose=False):
    """
    Merge stations from s2 into s1

    Ignore stations which do already exist in s1.
    """
    for net2 in s2.networks:
        if not any(net1.code == net2.code for net1 in s1.networks):
            s1.networks.append(net2)
            continue
        for net1 in s1.networks:
            if net1.code == net2.code:
                break
        for sta2 in net2.stations:
            if not any(sta1.code == sta2.code for sta1 in net1.stations):
                net1.stations.append(sta2)
                if verbose:
                    print(sta2)


def get_stations(level='station'):
    fname = STATION_FNAME + '_bh' * (level == 'channel') + '.xml'
    kw = {'network': '_US-REF', 'location': '', 'channel': 'BH?',
          'starttime': UTC('2005-04-01'), 'endtime': UTC('2015-09-01'),
          'minlatitude': 20, 'maxlatitude': 50,
          'minlongitude': -130, 'maxlongitude': -65,
          'level': level}
    client = FSDNClient()
    stations = client.get_stations(**kw)
    print('num stations:', len(stations.get_contents()['stations']))
    kw['network'] = '_US-TA'
    stations2 = client.get_stations(**kw)
    merge_stations(stations, stations2)
    print('num stations:', len(stations.get_contents()['stations']))
    # use HH channel only if BH channel not availlable
    kw['network'] = '_US-REF'
    kw['channel'] = 'HH?'
    stations2 = client.get_stations(**kw)
    merge_stations(stations, stations2)
    print('num stations:', len(stations.get_contents()['stations']))
    kw['network'] = '_US-TA'
    stations2 = client.get_stations(**kw)
    merge_stations(stations, stations2)
    print('num stations:', len(stations.get_contents()['stations']))
    # get stations with loc code 00 if no other station
    kw['network'] = '_US-REF'
    kw['channel'] = 'BH?'
    kw['location'] = '00'
    stations2 = client.get_stations(**kw)
    merge_stations(stations, stations2)
    print('num stations:', len(stations.get_contents()['stations']))
    kw['network'] = '_US-TA'
    stations2 = client.get_stations(**kw)
    merge_stations(stations, stations2)
    print('final num stations:', len(stations.get_contents()['stations']))
    stations.write(fname, 'STATIONXML')
    if level == 'station':
        stations.plot(projection='local', size=1,
                      outfile=STATION_FNAME + '.pdf')
    return stations


def events2json():
    events_xml = read_events(EVENT_FNAME + '.xml')
    events = OrderedDict()
    for event in events_xml[::-1]:
        evid = str(event.resource_id).split('/')[-1]
        ori = event.preferred_origin() or event.origins[0]
        mag = event.preferred_magnitude() or event.magnitudes[0]
        events[evid] = (str(ori.time), ori.latitude, ori.longitude, ori.depth,
                        mag.mag, mag.magnitude_type)
    with open(EVENT_FNAME + '.json', 'w') as f:
        json.dump(events, f)


stations = get_stations()
print('fetched stations')
stations2 = get_stations(level='channel')
print('fetched stations with channels')
print('\nEvents and picks of the ANF event catalog were downloaded via '
      'IRIS FDSN webservice. Unfortunately, IRIS is shutting down its event '
      'service '
      '(see e.g. http://ds.iris.edu/message-center/thread/1851/#m-3321).'
      ' The ANF catalog is available at http://anf.ucsd.edu/tools/events/ '
      'in CSS format but not in QuakeML format.\n'
      'The event data in QuakeML (90 MB of packed XML files) is needed for '
      'the us_qopen.py script. If needed, download the event metadata from '
      'https://github.com/trichter/misc/raw/master/data/us_events.7z '
      'and put the contents of the archive into the data/events folder.')
#events = get_events(individual=False)
#print('fetched all events without picks')
#events2json()
#print('created events json file')
#events = get_events()
#print('fetched all events')
