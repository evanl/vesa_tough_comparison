#Author - Evan Leister
import sys
import string
import t2_input_funcs as it2f
# TODO Nonobvious parameters in write_param
# TODO Capillary pressure fix? 
# TODO Collate all input parameters to be called from this driver
# TODO clean up processing routines into a returnable object that can be 
# called by other routines. 
# TODO comment post-processing routines

if __name__ == '__main__':
    """ t2_create_input is run in the following way from a command line 
        or shell script: 

        $ python t2_create_input.py <simulation title> 
    """
    if len(sys.argv) != 2:
        sys.exit( "Please specify a simulation title")
    sim_title = str(sys.argv[1])

    # if hydro_directory is False, it means that the initial condition will be 
    # generated using a hydrostatic gradient with a constant density.
    # If a directory is specified, the initial pressures and 
    # dissolved fractions will be taken from the 'hydro_directory' + '_dir/'.
    # in order to run a particular grid to hydrostatic equilibrium, 
    # specify hydro = True and use num_timesteps = 1 and days_per_timestep = 30
    hydro = False

    # if True, a uniform rectangular grid is generated.
    uniform = True

    # makes two-d case 
    two_d = False

    # Uses the IEAGHG data for Layer 9
    sleipner = False

    # Leaves the shale cells in for the sleipner case
    shale = True

    # If true, creates linear rel perms
    linear_rp = False

    # If True, removes capillary pressure (sets max cap of 10 Pa
    no_cap = False

    # If isothermal == True, the heat equation is not solved and 
    # the co2_enthalpy (specific enthalpy) is not called by TOUGH2
    isothermal = False
    co2_enthalpy = 293.e3 #J/kg

    # If True, ignores flow rate and fixes pressure and saturation
    # with an apparent CO2 saturation of sat_frac
    type1_source = False
    sat_frac = 0.50

    # performs a number of checks to process
    if hydro == True:
        hydro_directory = False
        # specifies no flow boundary conditions on the edges
        edge_bc_type = 2
    else:
        # directory names for different cases
        hydro_directory = 'u25_hydro'
        if sleipner == True:
            if shale == True:
                hydro_directory = 'sl_hydro'
            else:
                hydro_directory = 'sl_noshale_hydro_t32'
        if two_d == True:
            hydro_directory = 'sl_twod_hydro_32'

        # specifies fixed pressure boundary conditions on edge of the domain.
        edge_bc_type = 1
    if uniform == True:
        shale = True
        sleipner = False

    # specifies the number of output timesteps, and how many days are 
    # simulated in between each output.
    num_timesteps = 24
    days_per_timestep = 15

    t2inputgrid = it2f.write_input_files(sim_title, two_d = two_d, \
            uniform = uniform, \
            sleipner = sleipner, hydro = hydro, \
            hydro_directory = hydro_directory, \
            num_steps = num_timesteps, \
            days_per_step = days_per_timestep, \
            edge_bc_type = edge_bc_type, linear_rp = linear_rp,\
            no_cap = no_cap, shale = shale, tolerance = -5,\
            type1_source = type1_source, sat_frac = sat_frac,
            isothermal = isothermal, co2_enthalpy = co2_enthalpy)
