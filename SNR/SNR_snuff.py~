import os
import sys
import glob
import math
import time
import numpy
from pyrocko import util
from scipy import stats
from pyrocko import io
from pyrocko import trace
from pyrocko.snuffling import Snuffling
from obspy.core.utcdatetime import UTCDateTime
from obspy.signal.konnoohmachismoothing import konnoOhmachiSmoothingWindow as konno
from obspy.signal.konnoohmachismoothing import konnoOhmachiSmoothing as konno_spec
from obspy.signal.konnoohmachismoothing import calculateSmoothingMatrix as konno_mat


print "hallo"

def sampling_rate(self,tr_list):
    '''
        Checks if the sampling rate for all stations is the same.
    '''
    for i in range(1,len(tr_list)): 
        if tr_list[i].deltat != tr_list[0].deltat:
            self.fail('The selected stations have a different sampling rate!')



def swap_PSarrivals(p_arrival_time, s_arrival_time):
    '''
        Checks if the order of the P and S arrival is correct.
        If the order is incorrect, then swap the postition of the phase arrivals. 
    '''
    for i in range(len(p_arrival_time)):
        if p_arrival_time[i] > s_arrival_time[i]:
            p_arrival_time[i],s_arrival_time[i] = s_arrival_time[i],p_arrival_time[i]

def readMarker(self,phase_markers,event_marker):
    '''
        Creates the dictionary of phases for one event. The key is nslc - network, station, location, channel.
    '''

    ArrivalDict = {}
    keys_stations = []
    print "phase markers", phase_markers
    for i in range(len(phase_markers)):
        tmp = list(phase_markers[i].nslc_ids)
        print tmp
        sta_channel = tuple(tmp[0]) 
        keys_stations.extend([sta_channel])
    keys_stations = tuple(set(keys_stations)) 
    
    #set the keys for dictionary
    for key in keys_stations:
        ArrivalDict[key] = []
        print key

    #filling the dictionary
    for key in keys_stations:
        for i in range(len(phase_markers)):
            tmp = list(phase_markers[i].nslc_ids)
            sta_channel = tuple(tmp[0])
            if key == sta_channel:
                ArrivalDict[key].append((phase_markers[i].tmin))

    return ArrivalDict

class wave_data:
    '''
        Initialises the class for the P and S arrivals.
    '''
    def __init__(self):
    
        self.arrival_time    = None
        self.arrival_shift   = 0.0
        self.start_time      = None
        self.end_time        = None
        self.choped_trace    = []
        self.win_time        = None
        self.win_ydata       = None
        self.fftx            = []
        self.ffty            = []
        self.station         = []        
        self.channel         = []    



def WindowTime(evm, p, s, tr_list, window_length, ArrivalDict):
    '''
        Cuts the time windows for each trace according to the P and S arrival time.
    '''

    p.arrival_time  = numpy.zeros(len(tr_list))
    p.end_time      = numpy.zeros(len(tr_list))
    p.start_time    = numpy.zeros(len(tr_list))

    s.arrival_time  = numpy.zeros(len(tr_list))
    s.end_time      = numpy.zeros(len(tr_list))
    s.start_time    = numpy.zeros(len(tr_list))

    #for tr in tr_list:
    #    print "nslc:", tr.nslc_id, tr.nslc_id[0] + tr.nslc_id[1] + tr.nslc_id[2] + tr.nslc_id[3]

    for i in range(len(tr_list)):
        for key in ArrivalDict:
            nslc = tr_list[i].nslc_id
            if key == tr_list[i].nslc_id:
                for ind in range(len(ArrivalDict[key])):
                    if tr_list[i].tmin < evm.tmin and evm.tmin < tr_list[i].tmax: 
                        if ind % 2 == 0: 
                            p.arrival_time[i] = ArrivalDict[key][ind]
                        else:    
                            s.arrival_time[i] = ArrivalDict[key][ind]

                        

    swap_PSarrivals(p.arrival_time, s.arrival_time)

    p.arrival_shift = 0.1*window_length
    p.end_time      = p.arrival_time - p.arrival_shift;
    p.start_time    = p.end_time - window_length; 
        
    s.arrival_shift = 0.1*window_length;
    s.start_time    = s.arrival_time - s.arrival_shift;
    s.end_time      = s.start_time + window_length;

    print ArrivalDict

def stations_swap(sta_list, selected_sta1, selected_sta2):
    ''' 
        Orders the stations for the ratio according to the order given by a user.
    '''

    if sta_list[0] != selected_sta1:
        sta_list[0], sta_list[1] = sta_list[1], sta_list[0] 

def ProcessTrace(self, evm, p, s, tr_list, selected_sta1, selected_sta2,list_of_channels):
    '''
        Calculates SNR between two stations.
    '''
    

    # to cut the data using determined windows
    for i in range(len(tr_list)): 
        if tr_list[i].tmin < evm.tmin and evm.tmin < tr_list[i].tmax: 
            print "the trace item and station", i, tr_list[i].station
            print "start time p and s:", p.start_time[i], s.start_time[i]
            
            try:
                tmp_p = tr_list[i].chop(p.start_time[i], p.end_time[i], inplace=False, want_incomplete=False)
                tmp_s = tr_list[i].chop(s.start_time[i], s.end_time[i], inplace=False, want_incomplete=False)
                
                n = min(tmp_p.data_len(), tmp_s.data_len())
                tmp_p.set_ydata(tmp_p.get_ydata()[:n])
                tmp_s.set_ydata(tmp_s.get_ydata()[:n])
                

                p.choped_trace.append(tmp_p) 
                s.choped_trace.append(tmp_s) 
                p.station.append(tmp_p.station) 
                s.station.append(tmp_s.station)
                p.channel.append(tmp_p.channel) 
                s.channel.append(tmp_s.channel)
            except (trace.NoData):
                self.fail("My Dear Snufflinger! Some of the traces can not be cropped.\nThere might be two main reasons for it:\n1. The length of the time window is too large (decrease the time window);\n2. There is a gap in the trace data for the event " + util.time_to_str(evm.tmin, format='%Y-%m-%d %H:%M:%S.3FRAC') + " (remove this event)");
                print "hehe hehe ---- the is no data!"
           

    #print 'length of the chopped traces:', len(p.choped_trace[0].ydata), len(p.choped_trace[1].ydata), len(p.choped_trace[2].ydata),len(p.choped_trace[3].ydata),len(p.choped_trace[4].ydata),len(p.choped_trace[5].ydata)
    #print "the length of the p.choped_trace", len(p.choped_trace)

    # DC offset
    for i in range(len(p.choped_trace)):
        p.choped_trace[i].ydata -= numpy.mean(p.choped_trace[i].ydata)
        s.choped_trace[i].ydata -= numpy.mean(s.choped_trace[i].ydata)
    # taper

    for i in range(len(p.choped_trace)):
        deltat_taper = 0.1 * (p.choped_trace[i].tmax - p.choped_trace[i].tmin)
        taper_p = trace.CosTaper(p.choped_trace[i].tmin, p.choped_trace[i].tmin + deltat_taper, p.choped_trace[i].tmax - deltat_taper, p.choped_trace[i].tmax)
        taper_s = trace.CosTaper(s.choped_trace[i].tmin, s.choped_trace[i].tmin + deltat_taper, s.choped_trace[i].tmax - deltat_taper, s.choped_trace[i].tmax)
        p.choped_trace[i].taper(taper_p, inplace=True)
        s.choped_trace[i].taper(taper_s, inplace=True)

    #FFT and the amplitude spectra (can be done together with tapering, if to use the parameter tfade)
    bandwidth = 40
    print "-----BANDWIDTH"
    for i in range(len(p.choped_trace)):
        tmp_x_p, tmp_y_p = p.choped_trace[i].spectrum( pad_to_pow2=False, tfade=None)
        tmp_x_s, tmp_y_s = s.choped_trace[i].spectrum( pad_to_pow2=False, tfade=None)

        x_p_soothing_mat = konno_mat(tmp_x_p, bandwidth, normalize=False)
        x_s_soothing_mat = konno_mat(tmp_x_s, bandwidth, normalize=False)

        tmp2_y_p = numpy.dot(abs(tmp_y_p), x_p_soothing_mat)
        tmp2_y_s = numpy.dot(abs(tmp_y_s), x_s_soothing_mat)

        p.fftx.append(tmp_x_p), p.ffty.append(tmp2_y_p)
        s.fftx.append(tmp_x_s), s.ffty.append(tmp2_y_s)

    #average the ratios for each stations over the components
    #marking the stations
    sta = []
    ista1 = [] #borehole
    ista2 = [] #surface
    sta = list(set(p.station))
    print ">>>>>>>>> station list: ", sta
    stations_swap(sta, selected_sta1, selected_sta2)
    sta_tmp = list(p.station)

    for i in range(len(p.station)):
        if p.station[i] == sta[0]:
            ista1.extend([i])
        else:
            ista2.extend([i])

    ffty_sta1 = numpy.zeros(len(p.fftx[0]))
    ffty_sta2 = numpy.zeros(len(p.fftx[0]))
    ffty_sta12 = numpy.zeros(len(p.fftx[0]))
    

    if len(list_of_channels) == 1:
        ffty_sta1 = s.ffty[ista1[0]]/p.ffty[ista1[0]]
        ffty_sta2 = s.ffty[ista2[0]]/p.ffty[ista2[0]]
    else:
        ffty_sta1 = stats.gmean((s.ffty[ista1[0]]/p.ffty[ista1[0]], s.ffty[ista1[1]]/p.ffty[ista1[1]], s.ffty[ista1[2]]/p.ffty[ista1[2]]))
        ffty_sta2 = stats.gmean((s.ffty[ista2[0]]/p.ffty[ista2[0]], s.ffty[ista2[1]]/p.ffty[ista2[1]], s.ffty[ista2[2]]/p.ffty[ista2[2]]))  

    
        
    ffty_sta12 = ffty_sta1/ffty_sta2

    #print "test the lesngth of the ffty_sta12", s.ffty[ista1[0]], p.ffty[ista1[0]]
    
    #log10_ffty_sta12 = [ math.log10(x) for x in ffty_sta12]
    log10_ffty_sta12 = numpy.log10(ffty_sta12)
   
    
    return log10_ffty_sta12, p.fftx[0]

def std_s(i_event, std_SNR, Nevents):
    data_fn = "std_frequency_dependent_" + str(Nevents) + "events.txt" 
    #print "######################", i_event  
    if  i_event == 0:
        file = open(data_fn, "w")
    else:
        file = open(data_fn, "a")
        std_SNR = 20.0*numpy.std(std_SNR[0:i_event,:],axis=0)
        array_length = len(std_SNR)
        results_shaped_for_file = numpy.insert(std_SNR, 0, i_event+1)
        #print numpy.shape(results_shaped_for_file)
        #print numpy.shape(numpy.reshape(results_shaped_for_file,(1,array_length+1)))
        
        numpy.savetxt(file,numpy.reshape(results_shaped_for_file,(1,array_length+1)))
        #print i_event, numpy.shape(std_SNR)
    file.close()


