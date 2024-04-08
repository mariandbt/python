# S2 signal modules

from import_modules import *

import set_up as setup


# ________________________________________________________________________________________________________________
# Functions
# ________________________________________________________________________________________________________________


def set_global_parameters(global_vars, n_bb_files=None, n_bb_events_per_file=None, n_panels=None, n_sensors=None,
                          v_drift_EL=None, t_resolution_in_ns = 25):
    # If parameters are not provided, ask the user for input
    if n_bb_files is None:
        n_bb_files = int(input("Specify the number of bb data files: "))
    if n_bb_events_per_file is None:
        n_bb_events_per_file = int(input("Specify the number of bb events per file: "))
    if n_panels is None:
        n_panels = int(input("Specify the number of panels: "))
    if n_sensors is None:
        n_sensors = int(input("Specify the number of sensors: "))
    if v_drift_EL is None:
        v_drift_EL = input("Specify drift velocity of the EL region in [mm]/[ns]: ")
    # Add more parameters as needed

    dtheta = 2*np.pi/n_panels # rad
    dpos = int(n_sensors/n_panels) # number of sensors in 1 panel

    chunksize = int(2e5) # aprox length of an event to read the tables
    z_half_EL = -5.1 # [mm]

    t_res = t_resolution_in_ns # [ns]

    # Set global variables
    vars_names = ('n_bb_files', 'n_bb_events_per_file', 'n_panels', 'n_sensors', 'dtheta', 'dpos',
                  'chunksize','v_drift_EL', 'z_half_EL', 't_res')
    vars_values = (n_bb_files, n_bb_events_per_file, n_panels, n_sensors, dtheta, dpos,
                   chunksize, v_drift_EL, z_half_EL, t_res)

    setup.create_or_update_global_variable(global_vars, vars_names, vars_values, verbose = False)
    setup.create_or_update_global_variable(globals(), vars_names, vars_values, verbose = False)

    # Set more global variables as needed

    print("Global parameters set successfully :)")


def set_s2_table_specs(s2_table):

    s2_tab = s2_table[f'sens_200'] # all maps have the same specs, so we get whichever

    s2tab_x_nbins, s2tab_y_nbins = len(s2_tab.bin_initial_x.unique()), len(s2_tab.bin_initial_y.unique())

    s2tab_x_bin_width = (s2_tab.bin_final_x - s2_tab.bin_initial_x)[0]
    s2tab_y_bin_width = (s2_tab.bin_final_y - s2_tab.bin_initial_y)[0]

    s2tab_x_min = s2_tab.bin_initial_x.min()
    s2tab_y_min = s2_tab.bin_initial_y.min()

    s2tab_x_max = s2_tab.bin_final_x.max()
    s2tab_y_max = s2_tab.bin_final_y.max()


    # Set global variables
    vars_names = ('s2tab_x_nbins', 's2tab_y_nbins', 's2tab_x_bin_width',
                  's2tab_y_bin_width', 's2tab_x_min', 's2tab_y_min', 's2tab_x_max', 's2tab_y_max')
    vars_values = (s2tab_x_nbins, s2tab_y_nbins, s2tab_x_bin_width,
                   s2tab_y_bin_width, s2tab_x_min, s2tab_y_min, s2tab_x_max, s2tab_y_max)

    setup.create_or_update_global_variable(globals(), vars_names, vars_values, verbose = False)



def find_bb_ie(bb_file_path, start, event_id):

        stop = start + chunksize

        bb_particles = pd.read_hdf(bb_file_path, "/MC/particles",
                                   start = start, stop = stop,
                                   low_memory=True)

        bb_ie = bb_particles.query(f'(event_id == {event_id}) & (particle_name == "ie-")')
        prim_e = bb_particles.query(f'(event_id == {event_id}) & (primary == 1)')

        new_start = start + bb_particles.query(f'(event_id == {event_id})').index.max() - 1

        return new_start, bb_ie, prim_e

def find_rot(alpha_in_rad):
    rot = -10

    while not ((alpha_in_rad < (dtheta/2 + dtheta*rot)) &
               (alpha_in_rad > (-dtheta/2 + dtheta*rot))):
        rot += 1
        if rot > n_panels:
            error_coment = r'ERROR! Check that the angle is less than 2pi'
            print(error_coment)
            break

    return rot
find_rot = np.vectorize(find_rot) # Vectorize the function

def build_particle_dict(particles):
    # To build the dictionaries with the information of the ie⁻ of the bb event
    tt_dict = dict(zip(particles['particle_id'], particles['final_t']))
    xx_dict = dict(zip(particles['particle_id'], particles['final_x']))
    yy_dict = dict(zip(particles['particle_id'], particles['final_y']))
    zz_dict = dict(zip(particles['particle_id'], particles['final_z']))

    rr_dict = dict(zip(particles['particle_id'], np.sqrt(particles['final_x']**2 + particles['final_y']**2)))
    alpha_dict = dict(zip(particles['particle_id'], np.arctan2(particles['final_y'], particles['final_x'])))


    rr = np.array(list(rr_dict.values()))
    alpha = np.array(list(alpha_dict.values()))
    rotation = find_rot(alpha)

    new_alpha = alpha - rotation * dtheta
    new_xx = rr * np.cos(new_alpha)
    new_yy = rr * np.sin(new_alpha)


    rot_dict = dict(zip(particles['particle_id'], rotation))
    new_xx_dict = dict(zip(particles['particle_id'], new_xx))
    new_yy_dict = dict(zip(particles['particle_id'], new_yy))


    dict_names = ('tt_dict', 'xx_dict', 'yy_dict', 'zz_dict', 'rr_dict', 'alpha_dict', 'rot_dict', 'new_xx_dict', 'new_yy_dict')
    dicts = (tt_dict, xx_dict, yy_dict, zz_dict, rr_dict, alpha_dict, rot_dict, new_xx_dict, new_yy_dict)

    setup.create_or_update_global_variable(globals(), dict_names, dicts, verbose = False)


def build_s2_tab_dict(s2_table):

    s2_tab_dict = {}

    for sens_id in s2_table.keys():

        s2_tab = s2_table[sens_id]
        s2_matrix = s2_tab.s2.to_numpy().reshape(s2tab_x_nbins, s2tab_y_nbins)

        s2_tab_dict[sens_id] = s2_matrix

    setup.create_or_update_global_variable(globals(), 's2_tab_dict', s2_tab_dict, verbose = False)


def build_sensors_dict(sns_positions):
    # sensors
    x_dict = dict(zip(sns_positions['sensor_id'], sns_positions['x']))
    y_dict = dict(zip(sns_positions['sensor_id'], sns_positions['y']))

    theta_sens = np.round(np.arctan2(sns_positions['y'], sns_positions['x']).tolist(), 3)

    theta_dict = dict(zip(sns_positions['sensor_id'], theta_sens)) # theta vs sensor_id
    sens_dict = dict(zip(theta_sens, sns_positions['sensor_id'])) # sensor_id vs theta

    theta_to_pos_dict = dict(zip(sorted(theta_sens), np.arange(int(-n_sensors/2),
                                                               int(n_sensors/2)))) # theta vs position_id
    pos_to_theta_dict = dict(zip(np.arange(int(-n_sensors/2), int(n_sensors/2)),
                                 sorted(theta_sens))) # position_id vs theta

    dict_names = ('x_dict', 'y_dict', 'theta_dict', 'sens_dict', 'theta_to_pos_dict', 'pos_to_theta_dict')
    dicts = (x_dict, y_dict, theta_dict, sens_dict, theta_to_pos_dict, pos_to_theta_dict)

    setup.create_or_update_global_variable(globals(), dict_names, dicts, verbose = False)


def find_sensor(selected_sens, rot):

    theta = theta_dict[selected_sens]
    pos = theta_to_pos_dict[theta]

    new_pos = pos - rot*dpos

    if (new_pos == n_sensors/2):
        new_pos = -n_sensors/2

    if (new_pos > n_sensors/2):
        new_pos = new_pos%(n_sensors/2)

    if (new_pos < -n_sensors/2):
        new_pos = new_pos%(-n_sensors/2)

    new_theta = pos_to_theta_dict[new_pos]

    new_sens_id = sens_dict[new_theta]

    return new_sens_id


def find_s2(selected_sens, selected_particle):

    rr, alpha = rr_dict[selected_particle], alpha_dict[selected_particle]
    rot = rot_dict[selected_particle]
    new_xx, new_yy = new_xx_dict[selected_particle], new_yy_dict[selected_particle]

    new_sens_id = find_sensor(selected_sens, rot)
    s2_tab_matrix = s2_tab_dict[f'sens_{new_sens_id}']


    x_bin = int((new_xx - s2tab_x_min)//s2tab_x_bin_width)
    y_bin = int((new_yy - s2tab_y_min)//s2tab_y_bin_width)

    if ((new_xx > s2tab_x_max) or
        (new_xx < s2tab_x_min) or
        (new_yy > s2tab_y_max) or
        (new_yy < s2tab_y_min)
       ):
        s2_signal = 0.
    else:

        bin_content = s2_tab_matrix[x_bin][y_bin]
        s2_signal = np.random.poisson(bin_content, 1) # instead of just taking the contente of the bin, we add poisson fluctuations


    return s2_signal
find_s2 = np.vectorize(find_s2) # Vectorize the function





def create_s2_signal(s2_table, list_of_bb_file_paths, output_file_path, EL_ON = False):

    set_s2_table_specs(s2_table)
    build_s2_tab_dict(s2_table)

    file_index = 10**(int(math.log10(n_bb_events_per_file)) + 1)

    bin_width = t_res # [ns]

    # Open the HDF5 file in write mode
    with h5py.File(output_file_path, 'w') as file:

        for ii, bb_file_path in enumerate(list_of_bb_file_paths):

            bb_sns_pos, bb_sns_res = setup.read_fiber_sens(bb_file_path)

            build_sensors_dict(bb_sns_pos)

            start = 0

            for event in range(n_bb_events_per_file):

                start, bb_ie, prim_e = find_bb_ie(bb_file_path, start, event)

                if len(bb_ie) > 0:
                    group = file.create_group(str(file_index*ii + event))

                    build_particle_dict(bb_ie)

                    z_ie = np.array(list(zz_dict.values())) # [mm]
                    time_data = np.array(list(tt_dict.values())) # [ns]

                    t_delay = (z_ie - z_half_EL)/v_drift_EL # [ns]

                    time_data = time_data + t_delay # [ns]

                    bin_edges = np.arange(time_data.min() - bin_width, time_data.max() + 2*bin_width, bin_width)
                    t_values = (bin_edges[:-1] + bin_edges[1:])/2

                    prim_e_x = prim_e.initial_x.values[0] # [mm]
                    prim_e_y = prim_e.initial_y.values[0] # [mm]
                    prim_e_r = np.sqrt(prim_e_x**2 + prim_e_y**2) # [mm]

                    for jj, sens_id in enumerate(bb_sns_pos.sensor_id[:]):

                        if (((jj+1)%10 == 0) or jj == 0):
                            print(f'Sensor {jj+1}/{n_sensors}; Event {event+1}/{n_bb_events_per_file}; File {ii+1}/{n_bb_files}')

                        table_id = f'sens_{sens_id}'

                        s2_data = find_s2(sens_id, bb_ie['particle_id'])

                        # Create a histogram
                        s2_values, _ = np.histogram(time_data, bins=bin_edges,
                                                              weights = s2_data)


                        # Create the DataFrame after the loop using a dictionary
                        table_data = pd.DataFrame({'time': t_values, 's2': s2_values,
                                                #    'prim_e_x':prim_e_x, 'prim_e_y':prim_e_y,
                                                   'prim_e_r':prim_e_r
                                                   })

                        group.create_dataset(table_id, data=table_data)

    print('Done! s2 signal created :)')
