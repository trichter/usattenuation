# Copyright 2017 Tom Eulenfeld, MIT license
import logging
import functools
import os.path

import numpy as np
from obspy import read
from obspy.core import Stream, Trace


try:
    cache = functools.lru_cache(maxsize=512)
except AttributeError:
    def cache(f, ignore=None):
        """Caching decorator"""
        cache = {}

        @functools.wraps(f)
        def _f(*args, **kwargs):
            kws = tuple(sorted(kwargs.items()))
            key = args + kws
            try:
                return cache[key]
            except KeyError:
                cache[key] = result = f(*args, **kwargs)
                return result

        def cache_clear():
            _f.cache = {}
        _f.cache = cache
        _f.cache_clear = cache_clear
        return _f


def get_eventid(event):
    """Event id from event"""
    return str(event.resource_id).split('/')[-1]


def get_origin(event):
    """Preferred or first origin from event"""
    return event.preferred_origin() or event.origins[0]


def waveform_cache(path, win=None, format='MSEED', request_again=False):
    """Return caching decorator for waveforms
    gw_orig: original get_waveforms function
    path: filename with 4 or 3 '%s' (event mode or continuous mode)
      ('%s', will be filled with evid, seedid, win[0], win[1] in event mode
       '%s', will be filled with seedid, starttime, endtime in continous mode)
    win: length 2 tuple with window around origin time (event mode) or None
      (continuous mode)
    format: format to save waveforms in
    request_again: request waveforms again if no waveforms were obtained in
      previous runs.
    """
    assert win is None or len(win) == 2
    assert path.count('%') == 3 + (win is not None)

    def waveform_cache_decorator(gw_orig):

        @functools.wraps(gw_orig)
        def wrapper(network, station, location, channel, starttime, endtime,
                    event=None):
            assert (win is None) == (event is None)
            if event is None:
                st = starttime
                et = endtime
            else:
                etime = get_origin(event).time
                evid = get_eventid(event)
                st = etime + win[0]
                et = etime + win[1]
                if starttime < st or endtime > et:
                    log = logging.getLogger('waveform_cache')
                    log.error('Window has to be inside (%.1fs, %.1fs)', *win)
                    raise ValueError
            seedid = '.'.join((network, station, location, channel))
            if win is None:
                fname = path % (seedid, str(st)[:19], str(et)[:19])
            else:
                assert evid is not None
                fname = path % (evid, seedid, win[0], win[1])
            if os.path.exists(fname):
                stream = read(fname, format)
            if (not os.path.exists(fname) or
                    (request_again and len(stream) < 3)):
                args = (network, station, location, channel, st, et)
                try:
                    stream = gw_orig(*args)
                    stream.merge(method=1, fill_value='interpolate',
                                 interpolation_samples=-1)
                    if len(stream) == 0:
                        raise
                except (KeyboardInterrupt, SystemExit):
                    raise
                except Exception as ex:
                    stream = Stream(traces=[Trace(data=np.array([0.]))])
                    log = logging.getLogger('waveform_cache')
                    msg = 'channel %s: error while retrieving data: %s'
                    log.info(msg, seedid, ex)
                stream.write(fname, format)
            if len(stream) == 1 and len(stream[0]) == 1:
                raise ValueError('No data available')
            stream.trim(starttime, endtime)
            return stream
        return wrapper
    return waveform_cache_decorator
