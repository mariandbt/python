# S2 module to create bb events

from import_modules import *

import set_up as setup

import random
from pint import UnitRegistry

unit = UnitRegistry()
class InvalidUnitInputError(Exception):
    pass

class nexusEvent:
    def __init__(self, nexus_event_path, event_id):

        self.EventID    = event_id

        hits            = pd.read_hdf(nexus_event_path, "/MC/hits")
        hits            = hits.query(f'event_id == {self.EventID}')

        self.HitsX      = np.array(hits.x) *unit.mm
        self.HitsY      = np.array(hits.y) *unit.mm
        self.HitsZ      = np.array(hits.z) *unit.mm
        self.HitsTime   = np.array(hits.time) *unit.ns
        self.HitsEnergy = (np.array(hits.energy) *unit.MeV).to(unit.eV)

        self.PrimaryElectronX   = self.HitsX[self.HitsTime == self.HitsTime.min()]
        self.PrimaryElectronY   = self.HitsY[self.HitsTime == self.HitsTime.min()]
        self.PrimaryElectronR   = np.sqrt(self.PrimaryElectronX**2 + self.PrimaryElectronY**2)


    def AddDriftAndDiffusion(self, TPC, event_type = 'bb0nu'):

        hit_e   = self.HitsEnergy

        t_hit   = self.HitsTime # [ns]
        x_hit   = self.HitsX # [mm]
        y_hit   = self.HitsY # [mm]
        z_hit   = self.HitsZ # [mm]
        z_EL    = TPC.StartELPositionZ # [mm] Z position towards which the e⁻ drift 

        # Drift
        z_drift = np.fabs(z_EL - z_hit) # [mm]
        v_drift = TPC.ActiveDriftVelocity
        t_drift = (z_drift/v_drift).to(unit.ns)

        t_measurement = t_hit + t_drift # After drift

        recombi_term  = 1 - TPC.RecombinationFactor
        lifetime_term = np.exp(-t_drift/TPC.ElectronLifetime)

        if event_type == 'bb0nu':
            ionization_energy = TPC.XeIonization
        if event_type == 'KrCalibration':
            ionization_energy = TPC.KrIonization

        n_ie = hit_e*recombi_term*lifetime_term/ionization_energy

        self.NIonElectrons  = np.vectorize(int)(n_ie.magnitude)
        self.DriftTime      = t_drift

        # Diffusion
        drift_trans_diff = TPC.ActiveTransDiffusion
        drift_long_diff  = TPC.ActiveLongDiffusion

        sigma_trans = (drift_trans_diff*(z_drift**.5)).to(unit.mm) # [mm]
        sigma_long  = (drift_long_diff*(z_drift**.5)/v_drift).to(unit.ns) # [ns]

        self.TransDiff      = sigma_trans
        self.LongDiff       = sigma_long

        t_diff = np.array([])
        for t, std, n in zip(t_measurement, sigma_long, n_ie):
            mu      = t.magnitude
            sigma   = std.magnitude
            n       = int(n.magnitude)

            gaussian_diff = np.random.normal(mu, sigma, size = n)*t.units

            t_diff = np.concatenate((t_diff, gaussian_diff))
        t_diff = t_diff.to(unit.ns)

        x_diff = np.array([])
        for x, std, n in zip(x_hit, sigma_trans, n_ie):
            mu    = x.magnitude
            sigma = std.magnitude
            n     = int(n.magnitude)

            gaussian_diff = np.random.normal(mu, sigma, size = n)*x.units

            x_diff = np.concatenate((x_diff, gaussian_diff))
        x_diff = x_diff.to(unit.mm)

        y_diff = np.array([])
        for y, std, n in zip(y_hit, sigma_trans, n_ie):
            mu    = y.magnitude
            sigma = std.magnitude
            n     = int(n.magnitude)

            gaussian_diff = np.random.normal(mu, sigma, size = n)*y.units

            y_diff = np.concatenate((y_diff, gaussian_diff))
        y_diff = y_diff.to(unit.mm)

        total_n_ie  = self.NIonElectrons.sum()

        self.ElectronsIDs               = 10**(len(f'{total_n_ie}'))*self.EventID + np.arange(total_n_ie)
        self.ElectronsMeasurementTime   = t_diff.to(unit.ns)
        self.ElectronsFinalX            = x_diff.to(unit.mm)
        self.ElectronsFinalY            = y_diff.to(unit.mm)
        self.ElectronsFinalZ            = (z_EL*np.ones_like(self.ElectronsFinalY)).to(unit.mm)
        self.ElectronsFinalR            = np.sqrt(self.ElectronsFinalX**2 + self.ElectronsFinalY**2).to(unit.mm)
        self.ElectronsFinalAlpha        = np.arctan2(self.ElectronsFinalY, self.ElectronsFinalX).to(unit.rad)

        



class FiberBarrelTPC:
    def __init__(self, 
                 inner_rad          = 0 *unit.mm, 
                 outer_rad          = 500 *unit.mm, 
                 initial_phi        = 0*unit.rad, 
                 delta_phi          = 2*math.pi *unit.rad, 
                 length             = 3 *unit.m,
                 n_panels           = 18
                 ):

        # We unify the units for later calculations
        self.InnerRad    = inner_rad.to(unit.mm)
        self.OuterRad    = outer_rad.to(unit.mm) 
        self.InitialPhi  = initial_phi.to(unit.rad) 
        self.DeltaPhi    = delta_phi.to(unit.rad)   
        self.Length      = length.to(unit.mm)

        self.NPanels            = n_panels
        self.DeltaTheta         = 2*math.pi/self.NPanels  *unit.rad      

        self.bb0nuEnergyInXe  = (2458 *unit.keV).to(unit.eV)
        self.XeIonization     = (12.13 *unit.eV).to(unit.eV)

        self.KrDecayEnergy    = (41.5 *unit.keV).to(unit.eV)
        self.KrIonization     = (13.9996 *unit.keV).to(unit.eV)

    def SetActiveDriftVelocity(self, drift_velocity = 1e-3 *unit.mm/unit.ns):
        self.ActiveDriftVelocity = drift_velocity.to(unit.mm/unit.ns)

    def SetActiveTransDiffusion(self, drift_trans_diff = 1. *unit.mm/(unit.cm**.5)):
        # value reference: nexus simulation
        self.ActiveTransDiffusion = drift_trans_diff.to(unit.mm/(unit.mm**.5))

    def SetActiveLongDiffusion(self, drift_long_diff = .3 *unit.mm/(unit.cm**.5)):
        # value reference: nexus simulation
        self.ActiveLongDiffusion = drift_long_diff.to(unit.mm/(unit.mm**.5))

    def SetRecombinationFactor(self, recombi_factor = 0.026):
        # value reference: (https://iopscience.iop.org/article/10.1088/1748-0221/10/03/P03025/pdf)
        self.RecombinationFactor = recombi_factor

    def SetElectronLifetime(self, lifetime = 1e7 *unit.ns):
        self.ElectronLifetime = (lifetime).to(unit.ns)

    def SetEL(self, v_drift_EL = 2.5 *unit.mm/unit.us):
        self.StartELPositionZ   = 0 *unit.mm # start of EL
        self.HalfWidthEL        = 5.1 *unit.mm
        self.DriftVelocityEL    = v_drift_EL.to(unit.mm/unit.ns)


    def SetSensors(self, sns_path):
        sns_positions, _        = setup.read_fiber_sens(sns_path)
        self.NSensors           = len(sns_positions) # update number of sensors
        self.SensorsPerPanel    = int(self.NSensors/self.NPanels) 

        self.SensorsPosID   = np.arange((-self.NSensors/2), (self.NSensors/2))
        self.SensorsIDs     = np.array(sns_positions.sensor_id)
        self.SensorX        = np.array(sns_positions.x) *unit.mm
        self.SensorY        = np.array(sns_positions.y) *unit.mm
        self.SensorTheta    = np.round(np.arctan2(self.SensorY, self.SensorX), 3)

    def SetDefaults(self, sns_path):
        self.SetActiveDriftVelocity()
        self.SetActiveLongDiffusion()
        self.SetActiveTransDiffusion()
        self.SetRecombinationFactor()
        self.SetElectronLifetime()
        self.SetEL()
        self.SetSensors(sns_path)

        print('Default parameters set succesfully :) TPC ready to use')


class s2Table:
    def __init__(self, s2_table_path):

        self.FilePath   = s2_table_path

        s2_table    = setup.read_s2_table(self.FilePath)
        s2_tab      = s2_table['sens_200']

        self.SensorIDs  = s2_table.keys()

        self.NBinsX     = len(s2_tab.bin_initial_x.unique())
        self.NBinsY     = len(s2_tab.bin_initial_y.unique())

        self.WidthBinX  = (s2_tab.bin_final_x - s2_tab.bin_initial_x)[0] *unit.mm
        self.WidthBinY  = (s2_tab.bin_final_y - s2_tab.bin_initial_y)[0] *unit.mm
        
        self.MinimumX   = s2_tab.bin_initial_x.min() *unit.mm
        self.MinimumY   = s2_tab.bin_initial_y.min() *unit.mm

        self.MaximumX   = s2_tab.bin_initial_x.max() *unit.mm
        self.MaximumY   = s2_tab.bin_initial_y.max() *unit.mm

    def BuildS2TablesDict(self):

        s2_table    = setup.read_s2_table(self.FilePath)
        s2_tab_dict = {}

        for sens_id in self.SensorIDs:

            s2_tab = s2_table[sens_id]
            s2_matrix = s2_tab.s2.to_numpy().reshape(self.NBinsX, self.NBinsY)

            s2_tab_dict[sens_id] = s2_matrix
        
        self.S2TablesDict   = s2_tab_dict
    


def FindRotation(electron_final_alpha, TPC):

    dtheta = TPC.DeltaTheta.magnitude
    
    rot = np.where(electron_final_alpha > 0, 
                   int((electron_final_alpha + dtheta/2)/dtheta),
                   int((electron_final_alpha - dtheta/2)/dtheta)
                  )
    return rot
FindRotation = np.vectorize(FindRotation) # Vectorize the function

def BuildSensorsDict(TPC):
    
    x_dict      = dict(zip(TPC.SensorsIDs, TPC.SensorX.magnitude))
    y_dict      = dict(zip(TPC.SensorsIDs, TPC.SensorY.magnitude))

    sensor_id_sorted_by_theta   = TPC.SensorsIDs[np.argsort(TPC.SensorTheta.magnitude)]
    pos_to_id_dict              = dict(zip(TPC.SensorsPosID,  sensor_id_sorted_by_theta))
    id_to_pos_dict              = dict(zip(sensor_id_sorted_by_theta,  TPC.SensorsPosID))

    dict_names = ('x_dict', 'y_dict', 'pos_to_id_dict', 'id_to_pos_dict')
    dicts = (x_dict, y_dict, pos_to_id_dict, id_to_pos_dict)

    setup.create_or_update_global_variable(globals(), dict_names, dicts, verbose = False)

def FindSensor(sensor_id, rotation, TPC):

    pos     = id_to_pos_dict[sensor_id]

    new_pos = pos - rotation*TPC.SensorsPerPanel

    if (new_pos == TPC.NSensors/2):
        new_pos = -TPC.NSensors/2

    if (new_pos > TPC.NSensors/2):
        new_pos = new_pos%(TPC.NSensors/2)

    if (new_pos < -TPC.NSensors/2):
        new_pos = new_pos%(-TPC.NSensors/2)

    new_sens_id = pos_to_id_dict[new_pos]

    return new_sens_id
FindSensor = np.vectorize(FindSensor) # Vectorize the function

def FindS2(sensor_id, nexusEvent, s2Table, TPC):

    rotation    = FindRotation(nexusEvent.ElectronsFinalAlpha.magnitude, TPC)

    new_alpha   = nexusEvent.ElectronsFinalAlpha - rotation * TPC.DeltaTheta
    new_xx      = nexusEvent.ElectronsFinalR * np.cos(new_alpha)
    new_yy      = nexusEvent.ElectronsFinalR * np.sin(new_alpha)

    new_sens_id     = FindSensor(sensor_id, rotation, TPC)
    
    x_bin = ((new_xx - s2Table.MinimumX)//(s2Table.WidthBinX)).magnitude.astype(int)
    y_bin = ((new_yy - s2Table.MinimumY)//(s2Table.WidthBinY)).magnitude.astype(int)
    
    zero_signal_conditions = (new_xx > s2Table.MaximumX) | \
                             (new_xx < s2Table.MinimumX) | \
                             (new_yy > s2Table.MaximumY) | \
                             (new_yy < s2Table.MinimumY)
    
    # Initialize the output array with zeros
    s2_signal = np.zeros_like(new_xx)

    # Apply the conditions
    valid_indices = ~zero_signal_conditions
    
    # For the valid indices, get the corresponding values from the s2Table.S2TablesDict
    signal_values = np.array([s2Table.S2TablesDict[f'sens_{sid}'][xb, yb] 
                              for sid, xb, yb 
                              in zip(new_sens_id[valid_indices], x_bin[valid_indices], y_bin[valid_indices])])

    # Generate Poisson-distributed random numbers for these values
    poisson_fluctuations = np.random.poisson(signal_values)

    # Place the Poisson values into the s2_signal array
    s2_signal[valid_indices] = poisson_fluctuations

    return s2_signal

def ResponseSiPM(q_in_pes, t, t0, tau):

    """
    NOTE: units of t, t0, tau and (tau * rise_time) must be the same
    """
    # SiPM response parameters
    # tau = 200   # [ns] Decay time constant
    rise_time = 1 # Rise time constant

    rise_term = 1 - np.exp(-(t - t0) / (tau * rise_time))
    decay_term = np.exp(-(t - t0) / tau)

    signal = (rise_term * decay_term)
    signal[t<t0] = 0
    signal_area = np.trapz(x = t, y = signal) or 1

    normalized_signal = q_in_pes*signal/signal_area

    return normalized_signal

class s2Signal:
    def __init__(self, s2Table, TPC, nexusEvent):

        s2Table.BuildS2TablesDict()
        BuildSensorsDict(TPC)
        nexusEvent.AddDriftAndDiffusion(TPC)

        if nexusEvent.NIonElectrons.sum() > 0:

            # BuildElectronsDict(nexusEvent, TPC)

            time_data   = nexusEvent.ElectronsMeasurementTime
            t_delay     = (TPC.Length/2 + TPC.HalfWidthEL)/TPC.DriftVelocityEL 
            t_delay     = t_delay.to(unit.ns)
            time_data   = (time_data + t_delay).to(unit.ns) # [ns]
            t_values    = time_data.magnitude.astype(np.float32)

            prim_e_r                = nexusEvent.PrimaryElectronR.to(unit.mm)
            self.PrimaryElectronsR  = prim_e_r.magnitude.astype(np.float32) # [mm]
            self.Time               = t_values # [ns]
            self.SensorResponse     = {}

            for jj, sens_id in enumerate(TPC.SensorsIDs[:]):

                if (((jj+1)%1 == 0) or jj == 0):
                    processing_message = f'Processing sensor {jj+1}/{TPC.NSensors}...'
                    print(processing_message + ' '*len(processing_message), end = '\r')

                s2_data     = FindS2(sens_id, nexusEvent, s2Table, TPC)
                s2_values   = s2_data.astype(np.float32) # [pes]

                self.SensorResponse[sens_id]    = s2_values # [pes]
            
            print(f'{processing_message} Signal created succesfully :)', end = '\r')


    def AddShapinAndSamplin(self, shapin_tau = 155 *unit.ns, samplin_rate = 25 *unit.ns, t_binin = 0.1 *unit.ns):

        shapin_tau      = shapin_tau.to(unit.ns).magnitude
        samplin_rate    = samplin_rate.to(unit.ns).magnitude
        t_binin         = t_binin.to(unit.ns).magnitude

        sensor_keys = list(self.SensorResponse.keys())
        n_sensors   = len(sensor_keys)

        time_data   = self.Time # [ns]
        prim_e_r    = self.PrimaryElectronsR # [mm]

        # for the s2 as deltas
        tail_in_ns  = shapin_tau*4 # [ns]
        bin_edges   = np.arange(time_data.min(), time_data.max() + tail_in_ns, t_binin)

        # for the Shaping
        bin_means               = (bin_edges[:-1] + bin_edges[1:])/2
        generic_sipm_response   = ResponseSiPM(1, bin_means, bin_means.mean(), shapin_tau)

        # for the Samplin
        samplin_step = int(samplin_rate//t_binin)

        for jj, sensor in enumerate(sensor_keys[:]):

            if (((jj+1)%1 == 0) or jj == 0):
                processing_message = f'Processing sensor {jj+1}/{n_sensors}...'
                print(processing_message + ' '*2*len(processing_message), end = '\r')

            s2_data     = self.SensorResponse[sensor] # [pes]

            # s2 as deltas
            s2_deltas, _    = np.histogram(time_data, bins=bin_edges, weights = s2_data)
            print(f'{processing_message} deltas DONE!' + ' '*10, end = '\r')

            # Shaping
            s2_data_shaped          = np.convolve(s2_deltas, generic_sipm_response, mode='same')
            print(f'{processing_message} shapin DONE!' + ' '*10, end = '\r')

            # Sample: Filter data
            s2_data_shaped_sampled  = s2_data_shaped[::samplin_step]
            s2_values               = s2_data_shaped_sampled
            print(f'{processing_message} samplin DONE!' + ' '*10, end = '\r')

            self.SignalShapedSampled    = s2_values # [pes]
            self.SamplinRate            = samplin_rate *unit.ns # [ns]
            self.ShapinDecayConstant    = shapin_tau *unit.ns # [ns]

        print(f'{processing_message} Shapin and samplin done succesfully :)' + ' '*10, end = '\r')






# class s2Signal:
#     def __init__(s2Table, nexusEventDict, TPC):

#         s2Table.BuildS2TablesDict()
#         BuildSensorsDict(TPC)

#         for event in nexusEventDict.keys():

#             nexusEvent    = nexusEventDict[event]
#             nexusEvent.AddDriftAndDiffusion(TPC)

#             if len(nexusEvent.NIonElectrons.sum()) > 0:

#                 # BuildElectronsDict(nexusEvent, TPC)

#                 time_data   = nexusEvent.ElectronsMeasurementTime
#                 t_delay     = (TPC.Length/2 + TPC.HalfWidthEL)/TPC.DriftVelocityEL 
#                 t_delay     = t_delay.to(unit.ns)
#                 time_data   = (time_data + t_delay).to(unit.ns) # [ns]
#                 t_values    = time_data.magnitude

#                 prim_e_r    = nexusEvent.PrimaryElectronR.to(unit.mm).magnitude
#                 event_data  = {}
#                 event_data['prim_e_r_in_mm']   = prim_e_r # [mm]
# ################ chek cómo guardar la info ###########################
#                 event_data['signal']   = {} 
# ################ chek cómo guardar la info ###########################

#                 for jj, sens_id in enumerate(TPC.SensorIDs[:]):

#                     if (((jj+1)%1 == 0) or jj == 0):
#                         print(f'Sensor {jj+1}/{TPC.NSensors}; Event {event+1}/{n_bb_events_per_file}; File {ii+1}/{n_bb_files}')

#                     s2_data     = FindS2(sens_id, nexusEvent, s2_table, TPC)
#                     s2_values   = s2_data # [pes]

#                     sensor_data = {}
#                     sensor_data['time_in_ns']       = t_values # [ns]
#                     sensor_data['s2_in_pes']        = s2_values # [pes]
                        





#         file_index = 10**(int(math.log10(n_bb_events_per_file)) + 1)

#         # Open the HDF5 file in write mode
#         with h5py.File(output_file_path, 'w') as file:

#             for ii, bb_file_path in enumerate(list_of_bb_file_paths):

#                 bb_sns_pos, _ = setup.read_fiber_sens(sns_path)


#                 start = 0

#                 for event in nexusEventDict.keys():

#                     nexusEvent    = nexusEventDict[event]

#                     if len(bb_ie) > 0:
#                         event_group = file.create_group(str(file_index*ii + event))

#                         BuildElectronsDict(nexusEvent, TPC)

#                         z_ie = np.array(list(zz_dict.values()), dtype=np.float32) # [mm]
#                         time_data = np.array(list(tt_dict.values()), dtype=np.float32) # [ns]

#                         t_delay = (z_ie - z_half_EL)/v_drift_EL # [ns]

#                         time_data   = time_data + t_delay # [ns]
#                         t_values    = time_data

#                         # bin_edges = np.arange(time_data.min() - bin_width, time_data.max() + 2*bin_width, bin_width)
#                         # t_values = (bin_edges[:-1] + bin_edges[1:])/2

#                         prim_e_x = prim_e.initial_x.values[0] # [mm]
#                         prim_e_y = prim_e.initial_y.values[0] # [mm]
#                         prim_e_r = np.sqrt(prim_e_x**2 + prim_e_y**2) # [mm]

#                         for jj, sens_id in enumerate(bb_sns_pos.sensor_id[:]):

#                             sensor_group = event_group.create_group(f'sens_{sens_id}')

#                             if (((jj+1)%1 == 0) or jj == 0):
#                                 print(f'Sensor {jj+1}/{n_sensors}; Event {event+1}/{n_bb_events_per_file}; File {ii+1}/{n_bb_files}')

#                             # table_id = f'sens_{sens_id}'

#                             s2_data     = find_s2(sens_id, bb_ie['particle_id'])
#                             s2_values   = s2_data
                            
#                             # Integrate: Create a histogram
#                             # s2_values, _ = np.histogram(time_data,
#                             #                             bins=bin_edges,
#                             #                             weights = s2_data_shaped)
#                             #                             # weights = s2_data)


#                             # Create the dictionary after the loop
#                             print('s2_in_pes max = ', s2_values.max(), 'prim_e_r = ', prim_e_r, 'samplin_rate_in_ns = ', samplin_rate_in_ns)
#                             sensor_data = {}
#                             sensor_data['time_in_ns']              = t_values # [ns]
#                             sensor_data['s2_in_pes']               = s2_values # [pes]
#                             sensor_data['prim_e_r_in_mm']          = prim_e_r # [mm]
#                             # sensor_data['bin_width_in_ns']       = samplin_rate_in_ns # [ns]
#                             # sensor_data['samplin_rate_in_ns']    = samplin_rate_in_ns # [ns]
#                             # sensor_data['geant4_t_binin_in_ns']  = geant4_t_binin_in_ns # [ns]
#                             # sensor_data['bin_width_in_ns']       = bin_width # [ns]


#                             for data_key, values in sensor_data.items():
#                                 sensor_group.create_dataset(data_key, data=values)

#         print('Done! s2 signal created :)')



    # def ShapinAndSamplin(self):
    # def CreateS2SignalFile(self):

        




