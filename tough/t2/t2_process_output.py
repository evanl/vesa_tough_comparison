#Author - Evan Leister
import sys
import numpy as np
import matplotlib.pyplot as plt
import t2_output_funcs as ot2f

if __name__ == '__main__':
    print "PROCESSING TOUGH2 OUTPUT FOR " + str(sys.argv[1])
    if len(sys.argv) != 2:
        sys.exit( "Please specify a simulation title")
    sim_title = sys.argv[1]
    # for hydrostatic simulations, prints the pressure difference as a function
    # of elevation
    hydro = False

    # these three parameters should match those in t2_create_input.py
    two_d = False
    sleipner = False
    shale = False

    # if this is run to process two-file parallel output. 
    parallel = False

    # if the simulation is split into two runs, concatenate both of them into
    # one file in the following way (sim_title = test)
    # Files: 
    # OUTPUT_DATA_FIRST_8_STEPS
    # OUTPUT_DATA_LAST_3_STEPS
    #   Concatenate files by using
    #   $ cat OUTPUT_DATA_FIRST_8_STEPS OUTPUT_DATA_LAST_3_STEPS > \
    #         OUTPUT_DATA_test
    # OUTPUT_FIRST_8_STEPS
    # OUTPUT_LAST_3_STEPS
    #   $ cat OUTPUT_FIRST_8_STEPS OUTPUT_LAST_3_STEPS > \
    #         OUTPUT_test
    # 
    # After doing this, set split = 8 in order to correctly count steps
    split = 0

    # If there are two mass-balance blocks for each timestep, leave this True
    double_balance = True

    # number of timesteps from TOUGH2 output file.
    num_outputs = 5

    grid, time_steps = ot2f.process_t2_output(sim_title, parallel, split = split,\
            double_balance = double_balance, num_outputs = num_outputs)

    if two_d == True:
        eleme = 'aH732'
        k_layer = 0
        xind = 65/2
        yind = 119/2
    else:
        # sleipner
        if sleipner == True:
            eleme = 'JH732'
            if shale == True:
                k_layer = 3
            elif shale == False:
                k_layer = 3.
            xind = 32
            yind = 77
        else:
            eleme = 'yB212'
            k_layer = 0
            xind = 12
            yind = 12

    # choose format for plots.
    fmt = 'png'
    for basis in range(1,4):
        ot2f.plot_planar_contours(grid, time_steps, sim_title, fmt, \
                two_d = two_d, sleipner = sleipner, \
                shale = shale, axis = basis,\
                i_in = xind, j_in = yind, k_in = k_layer)

    ot2f.plot_wellhead_pressure(grid, time_steps, \
            well_cell_id = eleme)

    if hydro == True:
        ot2f.plot_incon_change_vs_index(grid, time_steps, vs_elev = True)
        ot2f.check_3d_hydro_pressure(grid, time_steps)

    ot2f.plot_mass_balance(grid, time_steps)

    #ot2f.write_viscosity(grid, time_steps)

    ot2f.move_files(sim_title, fmt, parallel = parallel)
