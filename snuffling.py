from pyrocko.snuffling import Snuffling, Choice, Param
from pyrocko.model import Event
from SNR_snuff import *
from SNR_setup_header import SNRSetup, Results, load
import numpy
from pyrocko import util
from matplotlib import *
import SNR_snuff

reload(SNR_snuff)


class SNR(Snuffling):

    '''
    Signal to noise ratio (SNR) calculation as a ration in decibels for two stations.

                SNR(station 1)
    SNR = --------------------, 
                SNR(station 2)

    where SNR(station 1) = signal/noise.


    USER GUIDE:
        0. Download the traces;

        1. To select two DIFFERENT stations from the list of stations given in the menu ("Select Station (borehole)" and "Select Station (surface)" ); 

        2. There are two possibilities to to calculate the SNR in the menu:
            a) "Active Event Marker (one event)" for SNR using a single event;
            b) "All Events Markers (several events)" for the SNR of the several events (all events on the traces);
        Please, select one of the options in the "Number of events" menu.
        Obtain markers by picking the Events and Phases or by downloading created earlier markers from the file. 
        One event marker should correspond to 6 phase markers (P and S markers for three components of each selected station).
        If option "a" is chosen activate the event and associated with it phase markers by clicking first on the "Event" sign and pressing "e".
        If option "b" is selected you do not need to activate any markers, all the markers will be used in the calculation. 

        3. Press "Run" to start the calculations;

        4. Press "Save" to save obtaned results in the file (first column - frequencies, second column - SNR). 

    '''


    def setup(self):
        '''
            Creates a GUI menu.
        '''

        self.set_name('SNR - signal to noise ratio')
        self.add_parameter(Choice('Select Station (borehole)','selected_sta1', '<empty>', ['<empty>']))
        self.add_parameter(Choice('Select Station (surface)','selected_sta2', '<empty>', ['<empty>']))
        self.add_parameter(Choice('Channels', 'selected_channel', '<empty>', ['<empty>']))
        self.add_parameter(Choice('Number of Events','selected_marker_type', 
				'Active Event Markers (one event)', 
				['Active Event Markers (one event)', 
				'All Events Markers (several events)']))
        self.add_parameter(Param('Time Window [s]', 'window_length', 1., 0.08, 1200.0))    # for a Signal/Noise Selection [s]'
        self.add_parameter(Choice('Output file header format', 'output_header', 'yaml (full)', 
				('yaml (full)', 'alternative', 'without header')))
        
        self.add_trigger('Save', self.save)
        self.current_stuff = None
        self.set_live_update(False)


    def save(self):
        '''
            Saves obtatined calculation into a file.
        '''

        if not self.current_stuff:
            self.fail('Nothing to save.')

        SNR, freq, std_SNR, sampling_rate, i_tmin = self.current_stuff

        if self.selected_marker_type == 'All Events Markers (several events)':
            dir = 'SNR_' + self.selected_sta1 + '_' + self.selected_sta2 + \
		 '_'  + 'Nevents'+ str(len(self.eventsMarkers)) + '_' + \
		 util.time_to_str(self.eventsMarkers[i_tmin].tmin, format='%Y-%m-%d_%H-%M-%S') + '.dat' 
        else: 
            dir = 'SNR_' + self.selected_sta1 + '_' + self.selected_sta2 + \
		 '_'  + str(self.eventsMarkers[0].get_event().name) + '_' + \
		 util.time_to_str(self.eventsMarkers[0].tmin, format='%Y-%m-%d_%H-%M-%S') + '.dat' 

        data_fn = self.output_filename(caption='Save Data', dir = dir)        
        file = open(data_fn, "w")
        if self.output_header == "yaml (full)":
            setup = SNRSetup(
                number_of_events = len(self.eventsMarkers),
                stations = [self.selected_sta1, self.selected_sta2],
                sampling_rate = sampling_rate,
                components = [self.selected_channel],
                results_columns = "1.frequency [Hz]     2.SNR     3.std")
            results_shaped_for_file=zip( freq, SNR, std_SNR)
            results = Results(setup=setup, values=numpy.array(results_shaped_for_file,dtype=numpy.float))
            s = results.dump()
            file.write("%s" % s)
        elif self.output_header == 'alternative':
            file.write("# Signal to Noise Ratio improvment for " + \
			str(self.selected_sta1) + "/" + str(self.selected_sta2) + " stations.\n")
            file.write("# Station components:       " + str(self.selected_channel) + "\n" )
            file.write("# Station sampling rate:    " + str(sampling_rate) + "\n")
            file.write("# number of events:         " + str(len(self.eventsMarkers)) + "\n")
            file.write("# \n")
            file.write("# freq\tSNR\tstd\n" )
            for i in range(len(freq)):
                file.write(("%e %e %e\n" % (freq[i], SNR[i], std_SNR[i])))
        else:
            for i in range(len(freq)):
                file.write(("%e %e %e\n" % (freq[i], SNR[i], std_SNR[i])))
        
        file.close()


    def get_all_event_and_phase_markers(self):
        ''' 
            Gets and separates all event and phase markers, 
		if an active event is not marked.
        ''' 

        all_markers = self.get_markers() 
        evm_all = []
        phm_all = []
        for mk in all_markers:
            if str(type(mk)) == "<class 'pyrocko.gui_util.EventMarker'>":
                evm_all.append(mk)
            elif str(type(mk)) == "<class 'pyrocko.gui_util.PhaseMarker'>":
                phm_all.append(mk)
            else:
                self.fail('This marker type is not acceptable in the calculations!')
        return evm_all, phm_all


    def marker_upload(self):
        ''' 
            Gets phase and event markers according to the menu selection.
        '''

        if self.selected_marker_type == 'All Events Markers (several events)':
            evm, phm = self.get_all_event_and_phase_markers()
        else:
            tmp_evm, phm = self.get_active_event_and_phase_markers()
            evm = []
            evm.append(tmp_evm)
        return evm, phm


    def event_marker_dict(self,evm, phm):
        '''
            Creates a dictionary, where each event marker 
		will be assotiated with six P and S phases.  
        '''

        evm = sorted(evm, key=lambda x: x.tmin)
        EventDict = {}
        for i in range(len(evm)):
            EventDict[evm[i]] = []
            for ph in phm:
                if (i < len(evm)-1):
                    if (evm[i].tmin <= ph.tmin and ph.tmin < evm[i+1].tmin):
                        EventDict[evm[i]].append((ph))
                else:
                    if (evm[i].tmin <= ph.tmin):
                        EventDict[evm[i]].append((ph))
        return EventDict


    def auto_plotting(self, evm, freq, SNR, std_SNR, sampling_rate):
        
        fframe = self.figure_frame()
        fig = fframe.gcf()
        ax = fig.add_subplot(1, 1, 1)

        ax.tick_params(axis='both',which='minor')
        ax.set_xlim(0,sampling_rate/2)
        
        if self.selected_channel == "*":
            title = 'Stations: ' + self.selected_sta1 +'/'+ self.selected_sta2 \
		 + ';    Number of events: ' + str(len(evm)) + ';  Components: all'
        else:
            title = 'Stations: ' + self.selected_sta1 +'/'+ self.selected_sta2 \
		 + ';    Number of events: ' + str(len(evm)) + ';  Components: '\
		 + self.selected_channel

        zero_line = numpy.zeros(len(freq))
        ax.plot(freq, zero_line, ':k')
        ax.plot(freq, SNR, 'b', lw=2)
        if len(evm) > 0:
            ax.fill_between(freq, SNR+std_SNR, SNR-std_SNR, facecolor='blue', alpha=0.2)

        ax.set_title(title)
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel('S-signal ratio / noise ratio [dB]')
        ax.minorticks_on()

        fig.canvas.draw()

        
    def event_with_tmin(self):
        evm_min_time = 100000000000000000000
        i_tmin = 0
        evm = self.eventsMarkers
        for i in range(0, len(evm)):
            if evm_min_time > evm[i].tmin:
                evm_min_time = evm[i].tmin
                i_tmin = i
        return i_tmin


    def time_range_chopper(self, evm):
        self.window_length
        chopper_tmin = evm.tmin - 0.01
        chopper_tmax = evm.tmin + 180.0

        return chopper_tmin, chopper_tmax


    def call(self):
        '''
           Starts the calculations by the button "Run" clicking.
        '''
        
        pile = self.get_pile()
        chans = list(pile.channels)

        # checks if the selected stations are different
        if self.selected_sta1 == self.selected_sta2:
            self.fail('Selected stations should be different!')

        if self.selected_channel == '*':
            list_of_channels = chans
        else:
            list_of_channels = [self.selected_channel]

        evm, phm = self.marker_upload()

        EventDict = self.event_marker_dict(evm, phm)

        # iterations over the events#
        for i in range(len(evm)):
            arr_dict = readMarker(self,EventDict[evm[i]],evm[i])     
            chopper_tmin, chopper_tmax = self.time_range_chopper(evm[i])

            for tr_list in pile.chopper(tmin=chopper_tmin,
                                        tmax=chopper_tmax,
                                        trace_selector=lambda tr:  (tr.station in [self.selected_sta1, self.selected_sta2]) and (tr.channel in list_of_channels)):

                sampling_rate(self,tr_list)

                p = wave_data()
                s = wave_data()
                WindowTime(evm[i],p, s, tr_list,self.window_length,arr_dict)
                (SNR_tmp, freq) = ProcessTrace(self,evm[i],p,s,tr_list, 
					self.selected_sta1, self.selected_sta2, list_of_channels)
                
            if i == 0:
                SNR = numpy.array(SNR_tmp)
                std_SNR = numpy.zeros((len(evm),len(freq))) 
            else:
                SNR = SNR+numpy.array(SNR_tmp)
            std_SNR[i,:] = SNR_tmp 

            std_s(i, std_SNR, len(evm))

        std_SNR = 20.0*numpy.std(std_SNR,axis=0)
        SNR = 20.0*SNR/(len(evm)*1.0)

        sampl_rate = int(1/tr_list[0].deltat)
        self.auto_plotting(evm, freq, SNR, std_SNR, sampl_rate)
      
        # save results into a file  
        self.eventsMarkers = evm
        i_tmin = self.event_with_tmin()
        self.current_stuff = (SNR, freq, std_SNR, sampl_rate, i_tmin)

        print "The calculations are done, my Dear Snufflinger!"


    def panel_visibility_changed(self, bool):
        if bool:
            self.adjust_controls()
            self.enable_pile_changed_notifications()
        else:
            self.disable_pile_changed_notifications()


    def pile_changed(self):
        self.adjust_controls()


    def adjust_controls(self):
        p = self.get_pile()
        stas = sorted(list(p.stations))
        chans = list(p.channels)
        if not stas:
            stas = ['<empty>']
        if not chans:
            chans = ['<empty>']

        chans.insert( 0, '*')
        self.set_parameter_choices('selected_sta1', stas)
        self.set_parameter_choices('selected_sta2', stas)
        self.set_parameter_choices('selected_channel', chans)

        
def __snufflings__():
    return [SNR()]

