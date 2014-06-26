#Author - Evan Leister
import sys
import string
import t2_input_funcs as it2f

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

    # if hydro is False, it means that the initial condition will be 
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

    # these values are the five values required in the van-genuchten-mualem
    # model for capillary pressure in the tOUGH2 manual. 
    # If no capillary pressure is desired, set the 4th parameter to 1.e1
    cp_vals = [0.4, 0.0, 1.61e-3, 1.e7, 0.999]

    # initial dissolved CO2 fraction, if nonzero, minimizes dissolution
    solubility = 0.454104e-3
    #solubility = 0.0
    temp = 42

    # If isothermal == True, the heat equation is not solved and 
    # the co2_enthalpy (specific enthalpy) is not called by TOUGH2
    isothermal = False
    # TEMP [C] || Specific Enthalpy [J/kg]
    #    32    ||          587.e3
    #    37    ||          613.e3
    #    42    ||          652.e3
    #    47    ||          697.e3
    co2_enthalpy = 612.e3 #J/kg

    # If True, ignores flow rate and fixes pressure and saturation
    # with an apparent CO2 saturation of sat_frac
    type1_source = False
    sat_frac = 0.50

    # specifies the number of output timesteps, and how many days are 
    # simulated in between each output.
    num_timesteps = 5
    days_per_timestep = 60.

    # rock parameters
    porosity = 0.35
    permeability = 2.e-12

    phase = 'co2'
    # injection rate in kg/sec
    if hydro == True:
        mass_rate = 0.0
    else:
        if sleipner == True:
            # average mass influx rate for the sleipner injection 
            mass_rate = 4.475
        else:
            mass_rate = 1.00

    # performs a number of checks to process
    if hydro == True:
        hydro_directory = False
        # specifies no flow boundary conditions on the edges
        edge_bc_type = 2
    else:
        # directory names for different cases
        hydro_directory = 'u25_hydro'
        #hydro_directory = 'u25_t37_hydro'
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

    if uniform == True:
        injection_cell = 'yB212'
        nx = 25
        ny = 25
        nz = 25
        dx = 50
        dy = 50
        dz = 0.6
    else:
        e_cel, nx, ny, nz = re.read_eclipse()
        if self.two_d == True:
            injection_cell = 'aH732'
        else:
            injection_cell = 'JH732'


    # create an input file and input grid object
    t2input = it2f.T2Input(sim_title, cp_vals, two_d = two_d, \
              uniform = uniform, \
              sleipner = sleipner, hydro = hydro, \
              hydro_directory = hydro_directory, \
              num_steps = num_timesteps, \
              days_per_step = days_per_timestep, \
              edge_bc_type = edge_bc_type, linear_rp = linear_rp,\
              shale = shale, mass_rate = mass_rate, \
              tolerance = -5, type1_source = type1_source, sat_frac = sat_frac,
              isothermal = isothermal, co2_enthalpy = co2_enthalpy,
              porosity = porosity, permeability = permeability, 
              injection_cell = injection_cell, phase = phase)
    t2grid = it2f.T2InputGrid(nx, ny, nz)

    # create and write the mesh file or use an old one
    
    meshmaker = False
    if meshmaker == False:
        t2input.write_mesh_file(t2grid, dx = dx, dy = dy, dz = dz, temp = temp,\
                solubility = solubility)
        t2input.write_incon(t2grid)

    t2input.write_input_file(sleipner = sleipner, meshmaker = meshmaker, \
            nx = nx, ny = ny, nz = nz)

    if hydro == True:
        t2grid.plot_cells(show = False, two_d = two_d)
        t2grid.plot_cells_slice(direc = 2, ind = 0, show = False)

