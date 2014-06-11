#Author - Evan Leister
import sys
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
import t2_output_funcs as ot2f

def process_t2_output(sim_title, parallel = False, split = 0,\
        double_balance = False):
    # this must be changed to reflect the input data
    # create an output grid object
    grid = ot2f.T2grid()
    # read the coordinate data from MESH
    grid.read_mesh()

    # number of timesteps from TOUGH2 output file.
    num_outputs = 24

    # read the file[s] for the timestep data
    if parallel == True:
        fname ='OUTPUT_DATA_' + sim_title
        gname = 'OUTPUT_' + sim_title
        print "fname " + fname
        print "gname " + gname
        f = open(fname)
        g = open(gname)
        grid.read_initial_mass(g)
    else:
        f = open(sim_title + '.out')
        g = open(sim_title) # JUST A PLACEHOLD, BAD CODING <<<<<<<
        grid.read_initial_mass(f)
    time_steps = []
    for i in range(num_outputs):
        if split !=0 and i > 7:
            year_add = 8
        else:
            year_add = 0
        timetest = ot2f.T2Timestep()
        double_read = double_balance
        timetest.readOutput(f, g, grid, parallel = parallel, \
                year_add = year_add, double_read = double_read)
        timetest.get_co2_density_viscosity()
        time_steps.append(timetest)
    f.close()
    g.close()

    return grid, time_steps

if __name__ == '__main__':
    print "PROCESSING TOUGH2 OUTPUT FOR " + str(sys.argv[1])
    if len(sys.argv) != 2:
        sys.exit( "Please specify a simulation title")
    sim_title = sys.argv[1]
    hydro = False
    two_d = False
    sleipner = False
    section = False
    shale = True
    parallel = False
    split = 0
    double_balance = True
    # creates a T2Grid object and a list of T2Timestep objects
    grid, time_steps = process_t2_output(sim_title, parallel, split = split,\
            double_balance = double_balance)
    # choose format for plots.
    fmt = 'png'
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
            if section == True:
                eleme = 'AE717'
                xind = 17
                yind = 47
                k_layer = 26
        else:
            eleme = 'yB212'
            k_layer = 0
            xind = 12
            yind = 12
        #modified for center
        #eleme = 'cA2 2'
        #k_layer = 2
        # column
        #eleme = 'dA0 0'
        #k_layer = 0
        # box

    for basis in range(1,4):
        ot2f.plot_planar_contours(grid, time_steps, sim_title, fmt, \
                two_d = two_d, sleipner = sleipner, section = section,\
                shale = shale, axis = basis,\
                i_in = xind, j_in = yind, k_in = k_layer)

    ot2f.plot_wellhead_pressure(grid, time_steps, \
            well_cell_id = eleme)

    if hydro == True:
        ot2f.plot_incon_change_vs_index(grid, time_steps, vs_elev = True)
        ot2f.check_3d_hydro_pressure(grid, time_steps)

    ot2f.plot_mass_balance(grid, time_steps)

    ot2f.write_viscosity(grid, time_steps)

    sim_title_dir = sim_title + '_dir'
    call(["mkdir",sim_title_dir])
    movestring = "mv " + "*" + fmt +" " + sim_title_dir
    call(movestring, shell=True)
    if parallel == False:
        call (["mv",sim_title,sim_title_dir])
        call (["mv",sim_title+'.out',sim_title_dir])
    else:
        fname ='OUTPUT_DATA_' + sim_title
        gname = 'OUTPUT_' + sim_title
        call (["mv",fname,sim_title_dir])
        call (["mv",gname,sim_title_dir])

    call (["mv",'INCON' ,sim_title_dir])
    call (["mv",'MESH' ,sim_title_dir])
    call (["mv",'SAVE' ,sim_title_dir])
    call (["mv",'rho_visc' ,sim_title_dir])
    call (["cp",'CO2TAB' ,sim_title_dir])
