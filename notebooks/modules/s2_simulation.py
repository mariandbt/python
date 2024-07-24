# S2 module to create bb events

from import_modules import *
import set_up as setup

unit = UnitRegistry()
class InvalidUnitInputError(Exception):
    pass

class nexusEvent:
    def __init__(self, nexus_event_path, event_id, event_type = 'bb0nu'):

        self.EventType  = event_type
        self.EventID    = event_id

        hits            = pd.read_hdf(nexus_event_path, "/MC/hits")
        hits            = hits.query(f'event_id == {self.EventID}')

        if hits.empty:
            print('This event is empty... Try another, this one won\'t take you far :(')
            self.Empty  = True
        
        else:
            self.Empty      = False
            self.HitsX      = np.array(hits.x) *unit.mm
            self.HitsY      = np.array(hits.y) *unit.mm
            self.HitsZ      = np.array(hits.z) *unit.mm
            self.HitsTime   = np.array(hits.time) *unit.ns
            self.HitsEnergy = (np.array(hits.energy) *unit.MeV).to(unit.eV)

            self.PrimaryElectronX   = self.HitsX[self.HitsTime == self.HitsTime.min()]
            self.PrimaryElectronY   = self.HitsY[self.HitsTime == self.HitsTime.min()]
            self.PrimaryElectronR   = np.sqrt(self.PrimaryElectronX**2 + self.PrimaryElectronY**2)


    def AddDriftAndDiffusion(self, TPC):

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

        recombi_prob    = 1 - TPC.RecombinationFactor
        lifetime_prob   = np.exp(-t_drift/TPC.ElectronLifetime).magnitude

        ionization_energy   = TPC.HPGXeIonization
        fano_factor         = TPC.HPGXeFano

        n_ie    = (hit_e/ionization_energy).magnitude
        sigma   = np.sqrt(n_ie*fano_factor)

        # Fano factor correction
        n_ie    = np.random.normal(loc = n_ie*recombi_prob, scale = sigma) # recombination
        n_ie    = np.random.normal(loc = n_ie*lifetime_prob, scale = sigma) # lifetime correction

        n_ie[n_ie < 0]  = 0
        n_ie            = np.vectorize(int)(n_ie)


        self.NIonElectrons  = n_ie
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

            gaussian_diff = np.random.normal(mu, sigma, size = n)*t.units

            t_diff = np.concatenate((t_diff, gaussian_diff))
        t_diff = t_diff.to(unit.ns)

        x_diff = np.array([])
        for x, std, n in zip(x_hit, sigma_trans, n_ie):
            mu    = x.magnitude
            sigma = std.magnitude

            gaussian_diff = np.random.normal(mu, sigma, size = n)*x.units

            x_diff = np.concatenate((x_diff, gaussian_diff))
        x_diff = x_diff.to(unit.mm)

        y_diff = np.array([])
        for y, std, n in zip(y_hit, sigma_trans, n_ie):
            mu    = y.magnitude
            sigma = std.magnitude

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


class HPGXeTPC:
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

        self.bb0nuEnergyInXe    = (2458 *unit.keV).to(unit.eV) # ref: (https://arxiv.org/abs/1804.01780)
        self.KrDecayEnergy      = (41.5 *unit.keV).to(unit.eV) # ref: (https://arxiv.org/abs/1804.01780)

        self.HPGXeIonization    = (22 *unit.eV).to(unit.eV) # ref: (https://doi.org/10.1016/j.nima.2009.10.076)
        self.HPGXeFano          = 0.15 # ref: (https://doi.org/10.1016/j.nima.2009.10.076)

    def SetActiveDriftVelocity(self, drift_velocity = 1 *unit.mm/unit.us):
        # ref: nexus simulation
        self.ActiveDriftVelocity = drift_velocity.to(unit.mm/unit.ns)

    def SetActiveTransDiffusion(self, drift_trans_diff = 1. *unit.mm/(unit.cm**.5)):
        # ref: nexus simulation
        self.ActiveTransDiffusion = drift_trans_diff.to(unit.mm/(unit.mm**.5))

    def SetActiveLongDiffusion(self, drift_long_diff = .3 *unit.mm/(unit.cm**.5)):
        # ref: nexus simulation
        self.ActiveLongDiffusion = drift_long_diff.to(unit.mm/(unit.mm**.5))

    def SetRecombinationFactor(self, recombi_factor = 0.026):
        # ref: (https://iopscience.iop.org/article/10.1088/1748-0221/10/03/P03025/pdf)
        self.RecombinationFactor = recombi_factor

    def SetElectronLifetime(self, lifetime = 1000 *unit.ms):
        # ref: nexus default
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
    def __init__(self, s2_table_path, light_yield = 1050, n_ie = 500*100e3):

        self.FilePath       = s2_table_path
        self.LightYield     = light_yield # [ph/ie⁻]
        self.Nie            = n_ie # number of ie⁻ used to create the table

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
    


def FindRotation(Electron_Final_Alpha, TPC):

    dtheta = TPC.DeltaTheta.magnitude
    
    rot = np.where(Electron_Final_Alpha > 0, 
                   int((Electron_Final_Alpha + dtheta/2)/dtheta),
                   int((Electron_Final_Alpha - dtheta/2)/dtheta)
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

def ResponseSiPM(q_in_pes, t, t0, rise_time, decay_time):

    """
    NOTE: units of t, t0, rise_time and decay_time must be the same
    """
    rise_term = 1 - np.exp(-(t - t0) / rise_time)
    decay_term = np.exp(-(t - t0) / decay_time)

    signal = (rise_term * decay_term)
    signal[t<t0] = 0
    signal_area = np.trapz(x = t, y = signal) or 1

    normalized_signal = q_in_pes*signal/signal_area

    return normalized_signal


def ConvertTomA(Waveform_in_pes_per_ns):

    sipm_gain = 4e6 # e/pes
    e_charge = 1.6e-19 # C/e
    
    pes_to_C = sipm_gain*e_charge # C/pes
    ns_to_s = 1e-9 # s/ns
    
    Waveform_in_A = Waveform_in_pes_per_ns*(pes_to_C/ns_to_s) # A
    Waveform_in_mA = Waveform_in_A*1e3 # mA

    return Waveform_in_mA
ConvertTomA = np.vectorize(ConvertTomA) # Vectorize the function


def ConvertTomV(Waveform_in_pes_per_ns, impedance_in_ohm = 50):

    Waveform_in_mV = ConvertTomA(Waveform_in_pes_per_ns)*impedance_in_ohm # mV

    return Waveform_in_mV
ConvertTomV = np.vectorize(ConvertTomV) # Vectorize the function



class s2Signal:
    def __init__(self, s2Table, TPC, nexusEvent):

        s2Table.BuildS2TablesDict()
        BuildSensorsDict(TPC)
        nexusEvent.AddDriftAndDiffusion(TPC)

        self.EventType          = nexusEvent.EventType
        self.EventID            = nexusEvent.EventID
        self.SensorResponse     = {}
        self.TotalNie           = nexusEvent.NIonElectrons.sum()

        prim_e_r                = nexusEvent.PrimaryElectronR.to(unit.mm)
        self.PrimaryElectronsR  = prim_e_r.magnitude.astype(np.float32)[0] # [mm]

        if nexusEvent.NIonElectrons.sum() > 0:

            time_data   = nexusEvent.ElectronsMeasurementTime
            t_delay     = (nexusEvent.ElectronsFinalZ + TPC.HalfWidthEL)/TPC.DriftVelocityEL 
            t_delay     = t_delay.to(unit.ns)
            time_data   = (time_data + t_delay).to(unit.ns) # [ns]
            t_values    = time_data.magnitude.astype(np.float32)

            self.Time   = t_values *unit.ns# [ns]

            for jj, sens_id in enumerate(TPC.SensorsIDs[:]):

                if (((jj+1)%1 == 0) or jj == 0):
                    processing_message = f'Creating signal (Processing sensor {jj+1}/{TPC.NSensors})...'
                    print(processing_message + ' '*len(processing_message), end = '\r')

                s2_data     = FindS2(sens_id, nexusEvent, s2Table, TPC)
                s2_values   = s2_data.astype(np.float32) # [pes]

                self.SensorResponse[sens_id]    = s2_values # [pes]
            
            print(f'{processing_message} Signal created succesfully :)', end = '\r')

        else:
            print(f'This event produced no signal!')
            self.Time   = 0 *unit.ns# [ns]
            for jj, sens_id in enumerate(TPC.SensorsIDs[:]):
                self.SensorResponse[sens_id]    = 0



    def AddShapinAndSamplin(self, shapin_rise = 12 *unit.ns, shapin_tau = 28 *unit.ns, samplin_rate = 25 *unit.ns, t_binin = 0.1 *unit.ns):

        shapin_rise     = shapin_rise.to(unit.ns).magnitude
        shapin_tau      = shapin_tau.to(unit.ns).magnitude
        samplin_rate    = samplin_rate.to(unit.ns).magnitude
        t_binin         = t_binin.to(unit.ns).magnitude

        if t_binin > samplin_rate:
            raise ValueError("The time bining cannot be > the sampling rate!!")

        sensor_keys = list(self.SensorResponse.keys())
        n_sensors   = len(sensor_keys)

        time_data   = self.Time.to(unit.ns).magnitude # [ns]

        # for the s2 as deltas
        tail_in_ns  = shapin_tau*4 # [ns]
        bin_edges   = np.arange(time_data.min(), time_data.max() + tail_in_ns, t_binin)

        # for the Shaping
        bin_means               = (bin_edges[:-1] + bin_edges[1:])/2
        generic_sipm_response   = ResponseSiPM(1, bin_means, bin_means.mean(), shapin_rise, shapin_tau)

        # for the Samplin
        samplin_step = int(samplin_rate//t_binin)

        self.SignalShapedSampled    = {} 

        for jj, sensor in enumerate(sensor_keys[:]):

            if (((jj+1)%1 == 0) or jj == 0):
                processing_message = f'Shaping and sampling (Processing sensor {jj+1}/{n_sensors})...'
                print(processing_message + ' '*2*len(processing_message), end = '\r')

            s2_data     = self.SensorResponse[sensor] # [pes]

            # s2 as deltas
            s2_deltas, _    = np.histogram(time_data, bins=bin_edges, weights = s2_data)

            # Shaping
            s2_data_shaped          = np.convolve(s2_deltas, generic_sipm_response, mode='same')

            # Sample: Filter data
            s2_data_shaped_sampled  = s2_data_shaped[::samplin_step]
            s2_values               = s2_data_shaped_sampled

            self.SignalShapedSampled[sensor]    = s2_values # [pes]/[ns]
            self.SamplinRate                    = samplin_rate *unit.ns # [ns]
            self.ShapinDecayConstant            = shapin_tau *unit.ns # [ns]

        print(f'{processing_message} Shapin and samplin done succesfully :)' + ' '*10, end = '\r')


    def PrintWaveform(self, 
                      sensor, 
                      shaped_and_sampled = True, 
                      units = 'mV',  # (mV, pes/ns, mA, ...) only used if shaped_and_sampled = True
                      impedance_in_ohm = 50, # only used if units_in_mV = True
                      bin_width = 1 *unit.us, # only used if shaped_and_sampled = False
                      new_figure = True, 
                      comment = ''):

        if new_figure:
                fig, ax = plt.subplots(nrows = 1, ncols = 1, 
                                       figsize=(7, 7), 
                                       constrained_layout=True) # Create a new figure

        else:
            # Check if there's an existing figure and create it if there's none
            if plt.gcf().get_axes():
                ax = plt.gcf().get_axes()[0]
            else:
                fig, ax = plt.subplots(nrows = 1, ncols = 1, 
                                           figsize=(7, 7), 
                                           constrained_layout=True)

        font_size   = 22
        event       = self.EventID

        if shaped_and_sampled:

            nt  = len(self.SignalShapedSampled[sensor])
            dt  = self.SamplinRate.to(unit.us).magnitude
            t0  = self.Time.to(unit.us).magnitude.min()
            t   = np.linspace(t0, t0 + nt*dt, nt)*unit.us # [us]

            if units == 'mV':
                waveform    = ConvertTomV(self.SignalShapedSampled[sensor], impedance_in_ohm)

            elif units == 'mA':
                waveform    = ConvertTomA(self.SignalShapedSampled[sensor])

            elif units == 'pes/ns':
                waveform    = self.SignalShapedSampled[sensor]

            else:
                raise ValueError("Invalid units: choose among mV, pes/ns or mA")

            samplin_rate    = self.SamplinRate
            tau             = self.ShapinDecayConstant
            label           = ("SiPM response: \n"
                               f"Sampling rate = {samplin_rate.magnitude:.2f}[{samplin_rate.units:~}] \n"
                               fr"$\tau$ = {tau.magnitude:.2f}[{tau.units:~}]"
                               )
            ax.plot(t.magnitude, waveform, label = label)

        else:

            t           = self.Time.to(unit.us) # [us]
            bin_width   = bin_width.to(unit.us) # [us]
            binin       = np.arange(t.magnitude.min() - bin_width.magnitude, t.magnitude.max() + 2*bin_width.magnitude, bin_width.magnitude)

            waveform    = self.SensorResponse[sensor]

            _, _, _ = ax.hist(t.magnitude, 
                              bins      = binin, 
                              weights   = waveform,
                              label     = 's2 waveform un-sampled and un-shaped'
                              );

            units = 'pes'
            

        ax.set_title(f's2 of event {event} in {sensor}', fontsize = font_size);
        ax.set_xlabel(f'Time [{t.units:~}]', fontsize = font_size);
        ax.set_ylabel(f'Signal [{units}]', fontsize = font_size);

        ax.tick_params(axis='both', labelsize = font_size*2/3);

        return t, waveform, ax


    def CreateS2DataFrames(self, 
                           event_id, 
                           units = 'pes/ns',
                           impedance_in_ohm   = 50 # only used if units = mV
                          ):

        # Create a dictionary for this event's signal
        # Event info
        event_df                                = {}
        event_df['event_id']                    = [event_id]
        event_df['N_ie']                        = [self.TotalNie] # Number of ionization electrons
        event_df['event_type']                  = [self.EventType]
        event_df['radial_coordinate_in_mm']     = [self.PrimaryElectronsR]

        # Waveforms
        s2_df                   = {}
        s2_df['event_id']       = event_id

        
        # s2 waveform
        if units == 'mV':
            waveform    = ConvertTomV(list(self.SignalShapedSampled.values()), impedance_in_ohm)

        elif units == 'mA':
            waveform    = ConvertTomA(list(self.SignalShapedSampled.values()))

        elif units == 'pes/ns':
            waveform    = list(self.SignalShapedSampled.values())

        else:
            raise ValueError("Invalid units: choose among mV, pes/ns or mA")

        # Time
        nt                  = len(waveform[0])
        dt                  = self.SamplinRate.to(unit.ns).magnitude
        t0                  = self.Time.to(unit.ns).magnitude.min()
        t                   = np.linspace(t0, t0 + nt*dt, nt) # [ns]
        s2_df['time']       = t.astype(np.float32)
        
        # Dynamically add sensor data with the appropriate suffix
        for sensor_id, wvf in zip(self.SignalShapedSampled.keys(), waveform):
            column_name             = f"sensor{sensor_id}_s2"
            s2_df[column_name]      = wvf.astype(np.float32)
        
        # Create a DataFrame from the event dictionary
        event_df    = pd.DataFrame(event_df)
        s2_df       = pd.DataFrame(s2_df)
        
        return event_df, s2_df
    

def CreateSignalHDF5(hdf5output_path, 
                     s2table,
                     TPC,
                     Events_Dict,
                     fiducial_radio     = 490   *unit.mm,
                     shapin_rise        = 12    *unit.ns,
                     shapin_tau         = 28    *unit.ns,
                     samplin_rate       = 25    *unit.ns,
                     t_binin            = 1     *unit.ns,
                     units              = 'pes/ns',
                     impedance_in_ohm   = 50 # only used if units = mV
                    ):

        
    # Create an empty DataFrame
    empty_df = pd.DataFrame()

    # Save the empty DataFrame to a .h5 file
    empty_df.to_hdf(hdf5output_path, key='/s2simulation/configuration', mode='w', format = 'table')

    # Signal
    event_id    = 0
            
    for jj, event_type in enumerate(Events_Dict.keys()):
        for ii, event_path in enumerate(Events_Dict[event_type]):
                # Read only the event_id column
                event_ids = pd.read_hdf(event_path, "/MC/hits", columns=['event_id'])
                # Get the maximum value of event_id
                n_events = event_ids['event_id'].max() + 1

                for event in range(n_events):
                # for event in range(0, 1):

                    if (((event)%1 == 0) or event == 0):
                        processing_message = f'Storing {event_type}\'s signal in a DF (Processing event {event + 1}/{n_events} in file {ii+1}/{len(Events_Dict[event_type])}...)'
                        print(processing_message + ' '*2*len(processing_message), end = '\r\n', flush=True)

                    NexusEvent = nexusEvent(event_path, event, event_type)
                    if NexusEvent.Empty:
                        continue
                    
                    if NexusEvent.PrimaryElectronR > fiducial_radio:
                        print(f'Event {NexusEvent.EventID} discarded by fiducial cut')
                        continue

                    # NexusEvent.AddDriftAndDiffusion(TPC) # already included in s2Signal initialization

                    s2signal = s2Signal(s2table, TPC, NexusEvent)
                    s2signal.AddShapinAndSamplin(shapin_rise, shapin_tau, samplin_rate, t_binin)

                    event_df, signal_df    = s2signal.CreateS2DataFrames(event_id, units, impedance_in_ohm)

                    setup.safe_write_to_hdf(event_df, hdf5output_path, '/s2simulation/events')
                    setup.safe_write_to_hdf(signal_df, hdf5output_path, '/s2simulation/s2')
                    
                    event_id    = event_id + 1
                
    # Configuration dataframe
    config_df                   = {}
    config_df['ParameterName']  = ['NEvents',
                                   'NSensors',
                                   'S2LightYield',
                                   'S2TableStatistics',
                                   'FiducialRadio',
                                   'Time (see s2 table)', 
                                   'Signal (see s2 table)', 
                                   'SamplingRate', 
                                   'ShapingRiseTime',
                                   'ShapingDecayTime'
                                  ]
    config_df['ParameterValue'] = [event_id,
                                   TPC.NSensors,
                                   s2table.LightYield,
                                   s2table.Nie,
                                   float(fiducial_radio.magnitude),
                                   None,
                                   None,
                                   float(samplin_rate.magnitude),
                                   float(shapin_rise.magnitude),
                                   float(shapin_tau.magnitude)
                                  ]
    config_df['ParameterUnits'] = np.vectorize(str)([unit.dimensionless,
                                                     unit.dimensionless,
                                                     'photons/ie⁻',
                                                     'ie⁻',
                                                     fiducial_radio.units,
                                                     unit.ns,
                                                     units,
                                                     samplin_rate.units,
                                                     shapin_rise.units,
                                                     shapin_tau.units
                                                    ])
    config_df = pd.DataFrame(config_df)

    # Store the configuration DataFrame using pandas
    setup.safe_write_to_hdf(config_df, hdf5output_path, '/s2simulation/configuration')
    
    
    print(f"Data has been written to {hdf5output_path} :)")


# **************************************************************************************************************************************
# **************************************************************************************************************************************
# **************************************************************************************************************************************
# **************************************************************************************************************************************
# **************************************************************************************************************************************
# **************************************************************************************************************************************
# **************************************************************************************************************************************
# **************************************************************************************************************************************
# **************************************************************************************************************************************
# **************************************************************************************************************************************
# **************************************************************************************************************************************
# **************************************************************************************************************************************
    
def CreateHDF5(self, hdf5output_path):

    df_configuration = pd.DataFrame({'ParameterName': ['NEvents', 'NSensors', 'DynRangeUnits','FiducialRadio', 'FiducialRadioUnits'],
                                        'ParameterValue': np.vectorize(str)([self.NEvents, self.NSensors, str(self.Units),
                                                                    self.FiducialR.magnitude, str(self.FiducialR.units)])
                                    })

    # Convert s2MaxValues to a DataFrame
    df_s2MaxValues = pd.DataFrame(list(self.s2TotalEnergy.items()), columns=['event', 'energy'])

    # Store the s2MaxValues DataFrame using pandas
    df_s2MaxValues.to_hdf(hdf5output_path, 
                            key='/s2simulation/TotalEnergy', 
                            mode='w', 
                            format='table', 
                            index=False,
                            data_columns=True
                            )
    

    df_configuration.to_hdf(hdf5output_path, key='/s2simulation/configuration', 
                            mode='a',
                            format='table', 
                            index=False,
                            data_columns=True)

    print(f"Data has been written to {hdf5output_path} :)")



class DynamicRange:
    def __init__(self, 
                 s2table, 
                 TPC, 
                 List_Events_Paths,
                 fiducial_radio     = 490 *unit.mm,
                 shapin_rise        = 25 *unit.ns,
                 shapin_tau         = 25 *unit.ns,
                 samplin_rate       = 25*unit.ns,
                 t_binin            = 1*unit.ns,
                 units              = 'mV',
                 impedance_in_ohm   = 50 # only used if units = mV
                 ):

        self.NSensors       = TPC.NSensors
        self.FiducialR      = fiducial_radio
        self.s2MaxValues    = {}

        for ii, event_path in enumerate(List_Events_Paths):
            # Read only the event_id column
            event_ids = pd.read_hdf(event_path, "/MC/hits", columns=['event_id'])
            # Get the maximum value of event_id
            n_events = event_ids['event_id'].max() + 1

            for event in range(n_events):
            # for event in range(0, 1):

                if (((event)%1 == 0) or event == 0):
                    processing_message = f'Dynamic range (Processing event {event + 1}/{n_events} in file {ii+1}/{len(List_Events_Paths)}...)'
                    print(processing_message + ' '*2*len(processing_message), end = '\r\n', flush=True)

                NexusEvent = nexusEvent(event_path, event)
                if NexusEvent.Empty:
                    continue
                
                if NexusEvent.PrimaryElectronR > fiducial_radio:
                    print(f'Event {NexusEvent.EventID} discarded by fiducial cut')
                    continue

                # NexusEvent.AddDriftAndDiffusion(TPC)

                s2signal = s2Signal(s2table, TPC, NexusEvent)
                s2signal.AddShapinAndSamplin(shapin_rise, shapin_tau, samplin_rate, t_binin)

                if units == 'mV':
                    waveform    = ConvertTomV(list(s2signal.SignalShapedSampled.values()), impedance_in_ohm) *unit.mV

                elif units == 'mA':
                    waveform    = ConvertTomA(list(s2signal.SignalShapedSampled.values())) *unit.mA

                elif units == 'pes/ns':
                    waveform    = list(s2signal.SignalShapedSampled.values()) *(unit.ns**-1)

                else:
                    raise ValueError("Invalid units: choose among mV, pes/ns or mA")

                max_value                   = max(map(max, waveform.magnitude))
                event_id                    = ii*(10**len(f'{n_events}')) + event
                self.s2MaxValues[event_id]  = max_value 

        self.Units      = waveform.units
        self.NEvents    = len(self.s2MaxValues.keys())


    def PrintHist(self, bin_width = 250, logscale = False):

        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7,7), constrained_layout=True)
        font_size = 20

        s2 = np.array(list(self.s2MaxValues.values()))
        n_events = len(s2)
        binin = np.arange(s2.min() - bin_width, s2.max() + 2*bin_width, bin_width)

        events, bins, bars = ax.hist(s2, binin,
                                     density=False,
                                     label='s2 max value in each event distribution',
                                     histtype='step')


        ax.text(0.6, .85, 'max value =%.2f'%(s2.max()),
        transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))

        ax.text(0.6, .8, '$\mu$ = %.2f'%(s2.mean()),
        transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))

        ax.text(0.6, .75, '$\sigma$ = %.2f'%(s2.std()),
        transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))

        ax.text(0.6, .7, '$N_{entries}$ = %s'%(int(events.sum())),
        transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))

        ax.text(0.6, .65, f'Fiducial radio cut = {self.FiducialR.magnitude:.2f} [{self.FiducialR.units:~}]',
        transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))


        ax.set_title(f'Max s2 signal of all {self.NSensors} sensors for {self.NEvents} events', fontsize = font_size);
        ax.set_xlabel(f's2 signal max [{self.Units:~}]', fontsize = font_size);
        ax.set_ylabel('Counts', fontsize = font_size);
        
        if logscale:
            ax.set_yscale('log')

        ax.legend(fontsize=0.7*font_size, loc='best')

        ax.tick_params(axis='both', labelsize = font_size*2/3)

        return events, bins, ax

    def CreateHDF5(self, hdf5output_path):

        df_configuration = pd.DataFrame({'ParameterName': ['NEvents', 'NSensors', 'DynRangeUnits','FiducialRadio', 'FiducialRadioUnits'],
                                         'ParameterValue': np.vectorize(str)([self.NEvents, self.NSensors, str(self.Units),
                                                                       self.FiducialR.magnitude, str(self.FiducialR.units)])
                                        })

        # Convert s2MaxValues to a DataFrame
        df_s2MaxValues = pd.DataFrame(list(self.s2MaxValues.items()), columns=['event', 's2max'])

        # Store the s2MaxValues DataFrame using pandas
        df_s2MaxValues.to_hdf(hdf5output_path, 
                              key='/s2simulation/DynamicRange', 
                              mode='w', 
                              format='table', 
                              index=False,
                              data_columns=True
                              )
        

        df_configuration.to_hdf(hdf5output_path, key='/s2simulation/configuration', 
                                mode='a',
                                format='table', 
                                index=False,
                                data_columns=True)

        print(f"Data has been written to {hdf5output_path} :)")


class EnergyResolution:
    def __init__(self, 
                 s2table, 
                 TPC, 
                 List_Events_Paths,
                 fiducial_radio     = 490 *unit.mm,
                 shaped             = False,
                 shapin_rise        = 25 *unit.ns,
                 shapin_tau         = 25 *unit.ns,     # only used if shaped == True   
                 samplin_rate       = 25*unit.ns,       # only used if shaped == True
                 t_binin            = 1*unit.ns,        # only used if shaped == True
                 units              = 'mV',
                 impedance_in_ohm   = 50                # only used if units == mV
                 ):

        self.NSensors       = TPC.NSensors
        self.FiducialR      = fiducial_radio
        self.s2TotalEnergy  = {}

        for ii, event_path in enumerate(List_Events_Paths):
            # Read only the event_id column
            event_ids = pd.read_hdf(event_path, "/MC/hits", columns=['event_id'])
            # Get the maximum value of event_id
            n_events = event_ids['event_id'].max() + 1

            for event in range(n_events):
            # for event in range(0, 1):

                if (((event)%1 == 0) or event == 0):
                    processing_message = f'Energy resolution (Processing event {event + 1}/{n_events} in file {ii+1}/{len(List_Events_Paths)}...)'
                    print(processing_message + ' '*2*len(processing_message), end = '\r\n', flush=True)

                NexusEvent = nexusEvent(event_path, event)
                if NexusEvent.Empty:
                    continue
                
                if NexusEvent.PrimaryElectronR > fiducial_radio:
                    print(f'Event {NexusEvent.EventID} discarded by fiducial cut')
                    continue

                # NexusEvent.AddDriftAndDiffusion(TPC)

                s2signal = s2Signal(s2table, TPC, NexusEvent)

                if shaped:
                    s2signal.AddShapinAndSamplin(shapin_rise, shapin_tau, samplin_rate, t_binin)

                    if units == 'mV':
                        waveform    = ConvertTomV(list(s2signal.SignalShapedSampled.values()), impedance_in_ohm) *unit.mV

                    elif units == 'mA':
                        waveform    = ConvertTomA(list(s2signal.SignalShapedSampled.values())) *unit.mA

                    elif units == 'pes/ns':
                        waveform    = list(s2signal.SignalShapedSampled.values()) *(unit.ns**-1)

                    else:
                        raise ValueError("Invalid units: choose among mV, pes/ns or mA")
                    
                    nt      = len(waveform[0])
                    dt      = s2signal.SamplinRate.to(unit.us).magnitude
                    t0      = s2signal.Time.to(unit.us).magnitude.min()
                    t       = np.linspace(t0, t0 + nt*dt, nt) # [us]
                    wvf     = waveform.magnitude

                else:
                    orderin         = np.argsort(s2signal.Time.magnitude)
                    waveform        = list(s2signal.SensorResponse.values()) *unit.dimensionless # [pes]
                    t               = s2signal.Time.to(unit.us).magnitude[orderin]
                    wvf             = np.stack(waveform.magnitude)[:, orderin]
                

                total_energy                    = np.trapz(x = t, y = wvf).sum()
                event_id                        = ii*(10**len(f'{n_events}')) + event
                self.s2TotalEnergy[event_id]    = total_energy 

        self.Units      = waveform.units
        self.NEvents    = len(self.s2TotalEnergy.keys())

    def PrintSpectrum(self, bin_width = 10, logscale = False):

        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7,7), constrained_layout=True)
        font_size = 20

        energy = np.array(list(self.s2TotalEnergy.values()))
        n_events = len(energy)
        binin = np.arange(energy.min() - bin_width, energy.max() + 2*bin_width, bin_width)

        events, bins, bars = ax.hist(energy, binin,
                                    density=False,
                                    label='Total energy value in each event distribution',
                                    histtype='step')


        ax.text(0.6, .85, 'max value =%.2f'%(energy.max()),
        transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))

        ax.text(0.6, .8, '$\mu$=%.2f'%(energy.mean()),
        transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))

        ax.text(0.6, .75, '$\sigma$=%.2f'%(energy.std()),
        transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))

        ax.text(0.6, .7, '$N_{entries}$ = %s'%(int(events.sum())),
        transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))

        ax.text(0.6, .65, f'Fiducial radio cut = {self.FiducialR.magnitude:.2f} [{self.FiducialR.units:~}]',
        transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))


        ax.set_title(f'Total energy of all {self.NSensors} sensors for {self.NEvents} events', fontsize = font_size);
        ax.set_xlabel(f'Event energy [{self.Units:~}]', fontsize = font_size);
        ax.set_ylabel('Counts', fontsize = font_size);
        
        if logscale:
            ax.set_yscale('log')

        ax.legend(fontsize=0.7*font_size, loc='best')

        ax.tick_params(axis='both', labelsize = font_size*2/3)

        return events, bins, ax

    
    def CreateHDF5(self, hdf5output_path):

        df_configuration = pd.DataFrame({'ParameterName': ['NEvents', 'NSensors', 'DynRangeUnits','FiducialRadio', 'FiducialRadioUnits'],
                                         'ParameterValue': np.vectorize(str)([self.NEvents, self.NSensors, str(self.Units),
                                                                       self.FiducialR.magnitude, str(self.FiducialR.units)])
                                        })

        # Convert s2MaxValues to a DataFrame
        df_s2MaxValues = pd.DataFrame(list(self.s2TotalEnergy.items()), columns=['event', 'energy'])

        # Store the s2MaxValues DataFrame using pandas
        df_s2MaxValues.to_hdf(hdf5output_path, 
                              key='/s2simulation/TotalEnergy', 
                              mode='w', 
                              format='table', 
                              index=False,
                              data_columns=True
                              )
        

        df_configuration.to_hdf(hdf5output_path, key='/s2simulation/configuration', 
                                mode='a',
                                format='table', 
                                index=False,
                                data_columns=True)

        print(f"Data has been written to {hdf5output_path} :)")


