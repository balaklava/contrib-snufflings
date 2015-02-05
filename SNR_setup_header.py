import numpy as num
from pyrocko.guts import Object, String, List, Float, \
			 Int, Tuple, Timestamp, load
from pyrocko.guts_array import Array


class Taper(Object):
    type = String.T()
    deltat_taper = String.T()


class Filter(Object):
    type = String.T()
    bandwidth = Float.T()


class SNRSetup(Object):
    number_of_events = Int.T()
    stations = List.T(String.T())
    components = List.T(String.T())
    sampling_rate = Float.T()
    time_window_length = Float.T()

    taper = Taper.T(default=Taper.D(type="trace.CosTaper", 
  				deltat_taper="0.1*WindowLength"))
    
    filter = Filter.T(default=Filter.D(type="obspy.signal.konnoohmachismoothing",
				bandwidth=40.0))
    
    results_columns = String.T()


class Results(Object):
    setup = SNRSetup.T()
    values = Array.T(dtype=num.float, shape=(None, 3), serialize_as='table')

