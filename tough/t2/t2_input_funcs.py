#Author - Evan Leister
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import read_eclipse as re
import eclipse_cells as ec
from time import clock
import string
from subprocess import call

class T2Input(object):
    def __init__(self, sim_title, cp_vals, two_d = False, uniform = False,\
            sleipner = False, hydro = False, hydro_directory = False,\
            num_steps = 11, days_per_step = 365.25, \
            edge_bc_type = 1, linear_rp = False, 
            shale = True, mass_rate = 0.0, tolerance = -5,\
            type1_source = False, sat_frac = 0.8, \
            isothermal = True, co2_enthalpy = 50., \
            porosity = 0.35, permeability = 2.e-12,\
            injection_cell = 'yB212', phase = 'co2'):
        """ Creates a TOUGH2 injection simulation
            self.sim_title 
                Title of simulation

            self.two_d 
                If True, makes a vertical average of a fully discretized 3d domain
                with a single, thick grid cell in the vertical.

            self.uniform 
                If true, a uniform, rectangular grid will be generated using
                fill_uniform_grid
              
            self.hydro 
                If hydro is False, it means that the initial condition will be 
                generated using a hydrostatic gradient with a constant density.
                If a directory is specified, the initial pressures and 
                dissolved fractions will be taken from the 
                'hydro_directory' + '_dir/'.
                in order to run a particular grid to hydrostatic equilibrium, 
                specify hydro = True and use num_timesteps = 1 and 
                days_per_timestep = 30

            self.hydro_directory 
                Directory with grid that has been run to hydrostatic equilibrium

            self.num_steps 
                Number of timesteps that are output to 'sim_title.out'

            self.days_per_step 
                The number of days in between each timestep
                Example:
                    In order to run a 360 day simulation with output roughly every
                    month, use:
                    num_steps = 12
                    days_per_step = 30

            self.edge_bc_type 
                If edge_bc_type is 2: all domain boundaries are no flow
                if edge_bc_type is 1: the horizontal domain boundaries are 
                fixed-pressure boundaries. the top and bottom boundaries
                of the domain remain no-flow

            self.linear_rp 
                IF true, relative permeability functions will be linear

            self.cp_vals 
                List of capillary pressure values, currently only takes a 
                van-genuchten-mualem model
                cp_vals[0] = lambda
                cp_vals[1] = liquid residual saturation, 
                             left at 0.0 for numerical purposes as the manual
                             recommends specifying residual saturation less than
                             the relative permeability residual saturation
                cp_vals[2] = 1/ P_0
                cp_vals[3] = P_max, set at 10. for maximum capillary pressure = 0.
                cp_vals[4] = S_b_sat, 

            self.shale 
                If true, includes the shale cells in the Eclipse files
            self.mass_rate 
                Constant injection rate [kg/s]
            self.tolerance 
                Absolute error convergence criterion for numerical solutions, 
                accepts error of 10^tolerance.
                Must be given as a negative number

            self.type1_source 
                If True, will ignore injection rate and specify a type 1 boundary
                condition at the injection cell

            self.sat_frac 
                Specifies the fraction of available CO2 saturation used in type 1
                boundary condition simulations.
                For example, if the brine residual saturation is 0.2, then 
                the CO2 can have 0.8/1 of the total pore space. 
                Thus, if sat_frac is 0.5, the CO2 will be fixed at 0.4 (0.8 * 0.5)
                in the type1_source cell

            self.isothermal 
                If True, no specific enthalpy is specified
                If False, the CO2 specific enthalpy is used and the injection cell
                is specified to have an infinite density

            self.co2_enthalpy 
                Specifies the CO2 specific enthalpy in J/kg

            porosity
                specifies default rock porosity
            perm 
                absolute permeability in [m^2]
            injection_cell
                specifies a single injection cell in 5-character coordinates

        """
        self.sim_title = sim_title
        self.two_d = two_d
        self.uniform = uniform
        self.hydro = hydro
        self.hydro_directory = hydro_directory
        self.num_steps = num_steps
        self.days_per_step = days_per_step
        self.edge_bc_type = edge_bc_type
        self.linear_rp = linear_rp
        self.cp_vals = cp_vals
        self.shale = shale
        self.mass_rate = mass_rate
        self.tolerance = tolerance
        self.type1_source = type1_source
        self.sat_frac = sat_frac
        self.isothermal = isothermal
        self.co2_enthalpy = co2_enthalpy
        self.porosity = porosity
        self.perm = permeability
        self.injection_cell = injection_cell
        self.phase = phase
        if type1_source == True:
            self.type1_source_cell = injection_cell
        else:
            self.type1_source_cell = 'none'
    def write_input_file(self, sleipner = False):
        print 'CREATING TOUGH2 INPUT FILE'


        output_day_list = []
        for i in range(self.num_steps):
            output_day_list.append(self.days_per_step * (i+1))

        print 'Creating T2 input files for simulation: ' + self.sim_title
        f = open(self.sim_title,'w')
        f.write('*'+ self.sim_title +'*\n')

        #input parameters
        if sleipner == True:
            self.porosity = 0.35 
            xperm = self.perm
            yperm = self.perm
            zperm = self.perm /3.
        else:
            xperm = self.perm
            yperm = self.perm
            zperm = self.perm 

        if self.linear_rp == True:
            rel_perm = 'linear'
            rp_vals = [0.2, 0., 1., 1.]
        else:
            rel_perm = 'vanGenuchten'
            rp_vals = [0.8, 0.2, 1.0, 0.05]

        plot_relperm_cap(rp_vals, self.cp_vals, fmt = 'pdf', rp = rel_perm)

        write_separator(f, 'ROCKS')
        # SAND
        name = 'sands'
        sand_density = 2600.
        write_rocks(f, name, sand_density, self.porosity, xperm, yperm, zperm, \
            self.cp_vals, rp_vals, thermk = 2.51, specheat = 920., \
            rel_perm = rel_perm)

        # injector cell to make the heat equation solve correctly
        name = 'well '
        sand_density = 2600.e40
        write_rocks(f, name, sand_density, self.porosity, xperm, yperm, zperm, \
            self.cp_vals, rp_vals, thermk = 2.51, specheat = 920., \
            rel_perm = rel_perm)

        # SHALE 
        name = 'shale'
        shale_density = 2600.
        shale_porosity = 0.15
        xperm = 1.e-18
        yperm = 1.e-18
        zperm = 1.e-18
        write_rocks(f, name, shale_density, shale_porosity, xperm, yperm, zperm, \
            self.cp_vals, rp_vals, thermk = 2.51, specheat = 920., \
            rel_perm = rel_perm,\
            end = True)

        write_multi(f, isothermal = self.isothermal)
        write_selec(f)
        write_start(f)
        write_param(f, tmax = output_day_list[-1] * 24 * 3600,\
                tolexp = self.tolerance)
        linear_solver_integer = 5
        preprocess_integer = 1
        write_solvr(f, linear_solver_integer, preprocess_integer )

        # flow rate in megatons per year for the Sleipner injection
        #massinflow = [0.0198, 0.0405, 0.0437, 0.0540, 0.0740, 0.1030, \
                                    #0.1390, 0.1830, 0.2370, 0.2960, 0.370]
        #kg_per_sec_factor = 31.71 # kg/s per mt/yr
        #kgInflow = [ x * kg_per_sec_factor for x in massinflow]

        if self.type1_source == True:
            write_gener(f, self.injection_cell, phase = self.phase, \
                    mass_rate = 0.0,\
                    co2_enthalpy = self.co2_enthalpy, \
                    isothermal = self.isothermal)
        else:
            write_gener(f, self.injection_cell, phase = self.phase, \
                    mass_rate = self.mass_rate, \
                    co2_enthalpy = self.co2_enthalpy, \
                    isothermal = self.isothermal)
                    #, kg_inflow = kgInflow, times = output_day_list )

        write_times(f, output_day_list)
        write_foft(f)
        write_coft(f, sleipner)
        write_goft(f)
        write_separator(f, 'ENDCY')
        f.close()

        print "write_input_file COMPLETE"
        return 0

    def write_mesh_file(self, t2grid, dx = 50, dy = 50, dz = 0.6, \
            temp = 32, inj_cell = 'yB212', solubility = 0.0):
        if self.uniform == True:
            brine_density = 1016.4
            t2grid.fill_uniform_grid(self.porosity, dx, dy, dz, \
                    density = brine_density, solubility_limit = solubility, \
                    inj_cell = self.injection_cell, temperature = temp)
            e_cel = 'uniform'
            t2grid.write_mesh(e_cel, two_d = self.two_d, uniform = self.uniform,\
                    boundary_type = self.edge_bc_type, shale = self.shale,\
                    type1_source_cell = self.type1_source_cell)
        else:
            brine_density = 1019.35
            if self.hydro_directory == False:
                t2grid.fill_eclipse_grid(e_cel, temperature = temp, \
                        density = brine_density, two_d = self.two_d, \
                        solubility = solubility, shale = self.shale)
                t2grid.write_mesh(e_cel, two_d = self.two_d, uniform = self.uniform,\
                        boundary_type = self.edge_bc_type, shale = self.shale,\
                        type1_source_cell = self.type1_source_cell)
            else:
                mesh_string = self.hydro_directory + '_dir/MESH'
                call (["mv", mesh_string, './OLDMESH'])
                clean_oldmesh_file()
        return 0

    def write_incon(self, t2grid):
        if self.hydro_directory == False:
            t2grid.write_incon(shale = self.shale, porosity = self.porosity)
        else:
            t2grid.use_old_incon(self.hydro_directory, \
                    type1_source_cell = self.type1_source_cell,\
                    saturation_fraction = self.sat_frac)

        return 0
    # end T2Input class

# writing functions for different blocks of the input file
def write_separator(f, keyword):
    """ writes for the following keywords
    START ROCKS
    MULTI PARAM
    SOLVR GENER
    TIMES FOFT 
    GOFT  COFT 
    ELEME INCON
    ENDCY CONNE
    MESHMAKER
    SELEC MOP
    ENDFI
    one space is required after four letter keywords
    none required after MOP
    """
    if keyword == 'ENDFI':
        f.write(''.join([keyword,'----']))
        f = write_dashes(f)
    elif keyword == 'MOP':
        f.write('----*----1 MOP: 123456789*123456789*1234')
        f.write(' ---*----5----*----6----*----7----*----8')
        f.write('\r\n')
    elif keyword == 'SELEC':
        f.write('SELEC....2....3....4....5....6....7....8')
        f.write('....9...10...11...12...13...14...15...16')
        f.write('\r\n')
    elif keyword =='MESHMAKER':
        f.write('MESHMAKER')
        write_dashes(f)
        f.write('\r\n')
    elif (  
            keyword == 'START' or 
            keyword == 'ROCKS' or 
            keyword == 'MULTI' or 
            keyword == 'PARAM' or 
            keyword == 'SOLVR' or 
            keyword == 'GENER' or 
            keyword == 'TIMES' or 
            keyword == 'FOFT ' or 
            keyword == 'GOFT ' or 
            keyword == 'COFT ' or
            keyword == 'ELEME' or 
            keyword == 'CONNE' or 
            keyword == 'INCON' or  
            keyword == 'ENDCY'    ):
        f.write(''.join([keyword,'----']))
        f = write_dashes(f)
        f.write('\r\n')
    else :
        print "keyword = " + keyword
        raise NameError('Invalid Keyword')
    return f

def write_dashes(f):
    f.write('1----*----2----*----3----*----4----*----5----*----6----*----7----*----8')
    return f

def write_rocks(f, name, density, porosity, xperm, yperm, zperm, \
    cp_vals, rp_vals, thermk = 2.51, specheat = 920., \
    cap = 'vanGenuchten', rel_perm = 'vanGenuchten',\
    end = False):
    """
        writes a ROCKS entry with the following properties
        name: Choose a name of the rock i.e. sand, shale
        density: rock grain density [kg/m^3]
        porosity: default porosity if no porosity is specified in INCON
        xperm: x-permeability [m^2]
        yperm: y-permeability [m^2]
        zperm: z-permeability [m^2]
        thermk: formation thermal conductivity [W/(m degC)]
        specheat: rock specific heat [J/(kg degC)]
            NOTE: values greater than 10^4 will be ignored.
        cp_vals: list of values for capillary pressure functions.
            NOTE In this case, only van-genuchten-mualem models are 
            used for capillary pressure. 
        rp_vals: list of values for relative permeability functions
        cap: string that specifies the type of capillary pressure function
        rel_perm: value that specifies the type of relative permeability curves
        end: if this is the final rock type used, set end = True
    """

    # NOTE, the 2 that is hard-coded here means that two materials will be
    # read, if only one is used, change this value to 1
    f.write(name + '    2')
    if density < 10.e7:
        d = '{: 10.0f}'.format(density)
    else:
        d = '{: 10.1e}'.format(density)
    p = '{: 10.2f}'.format(porosity)
    xp = '{: 10.2e}'.format(xperm)
    yp = '{: 10.2e}'.format(yperm)
    zp = '{: 10.2e}'.format(zperm)
    l = '{: 10.2f}'.format(thermk)
    sph = '{: 10.2f}'.format(specheat)
    f.write(d + p + xp + yp + zp + l + sph + '\r\n')

    # compressibility = 0
    f.write('   0.0e-10\r\n')

    # Capillary pressure
    cp_type = 7

    # relative permeability
    if rel_perm == 'vanGenuchten':
        rp_type = 7
    else:
        rp_type = 1

    # write either form
    rp_str = 4 * ' ' + str(rp_type) + 5 * ' '
    s = ""
    for el in rp_vals:
        s += format_float_rocks(el)
    f.write(rp_str + s + '\r\n')
    cp_str = 4 * ' ' + str(cp_type) + 5 * ' '
    sc = ""
    for el in cp_vals:
        sc += format_float_rocks(el)
    f.write(cp_str + sc + '\r\n')

    # break the rocks line
    if end == True:
        f.write('\r\n')
    return f
def write_multi(f, isothermal = True):
    if isothermal == True:
        write_separator(f, 'MULTI')
        f.write('    3    3    3    6\r\n')
    else:
        write_separator(f, 'MULTI')
        f.write('    3    4    3    6\r\n')
    return f

def write_selec(f):
    write_separator(f, 'SELEC')
    f.write('    1     ' \
        + 10 * ' ' + 10 * ' ' + 10 * ' ' + \
          '         0' + '    0    0' + \
                '    0    0' + '    0    0\r\n')
    f.write('        .8        .8\r\n')
    return f

def write_start(f):
    write_separator(f,'START')
    write_separator(f,'MOP')
    return f

def write_param(f, pres = 110.5e5, salt = 3.2e-2, \
        co2 = 0.0, temp = 32., maxtstep = 9500, max_cpu_seconds = 9999, \
    tmax = 63.1152e6, tolexp = -5):
    """
        Writes the PARAM block of the input file
        Default parameters for an element (grid cell) if not specified in INCON
            pres: pressure 
            salt: dissolved NaCl concentration
            co2: CO2 dissolved fraction
            temp: temperature
        maxtstep: maximum number of timesteps
        max_cpu_seconds: maximum number of cup seconds
        tmax: end time of simulation
        tolexp: convergence criterion relative error exponent, default = 10^-5
            example: 10^-5 = 10^tolexp
    """
    write_separator(f, 'PARAM')
    # line 1

    tstep = '{:4d}'.format(maxtstep)
    f.write(4 * ' ' + tstep + 4 * ' ')
    # maximum duration in CPU seconds
    f.write('{:4d}'.format(max_cpu_seconds))
    # writes the MOP parameters, see the TOUGH2 manual for details
    # on what each of them are.
    #   MOP: 123456789*123456789*1234
    f.write('1000 00000000  4    3   \r\n')

    # line 2
    f.write(10 * ' ')
    
    b = '{: 10.5e}'.format(tmax)
    l = list(b)
    # l[2], l[1] = l[1], l[2]
    del l[-2]
    del l[-2]
    c = ''.join(l)
    f.write(c)

    # length of time steps in seconds
    # this seems to be an initial time step, and is reduced if the 
    # system fails to converge. 
    f.write('     0.8e6')

    # gravity
    f.write(20 * ' ' + '      9.81\r\n')

    # convergence criterion relative error
    tolstr = '{:03d}'.format(tolexp)
    # convergence criterion maximum error
    f.write('    1.E' + tolstr + '    1.E-01\r\n')

    # default initial conditions
    s1 = '{: 10.3e}'.format(pres)
    s2 = '{: 10.3e}'.format(salt)
    s3 = '{: 10.3e}'.format(co2 )
    s4 = '{: 10.1f}'.format(temp)
    f.write(10 * ' ' + s1)
    f.write(10 * ' ' + s2)
    f.write(10 * ' ' + s3)
    f.write(10 * ' ' + s4)
    f.write('\r\n')
    return f

def write_solvr(f, linsolve_int, preprocess_int, closur = -7):
    """
        writes the SOLVR block for solving equations
        linsolve_int: integer [2,6] that solves the 
            system. 
            2: biconjugate gradient
            3: Lanczos-type biconjugate gradient
            4: generalized minimum residual solver
            5: stabilized biconjugate gradient solver
            6: Direct solver 
        preprocess_int: preprocessing integer, 
            0: no preprocessing
            1: replace zeros on main diagonal with small constant
            2: make linear combinations to achieve nonzeros on main diagonal
            3: normalize equations, followed by 2
            4:  details in manual
        closur: closure criterion for CG iterations
    """
    write_separator(f,'SOLVR')
    if closur < -12 or closur > -6:
        print "closur must be between -12 and -6"
        return 1
    tolstr = '{:03d}'.format(closur)
    # the O0 in this line says no preconditioner will be used.
    # 1.0e-1 is the default value for RITMAX, details in manual
    # closur: convergence criterion for 
    f.write(str(linsolve_int) + \
            '  Z'+ str(preprocess_int) + \
            '   O0    1.0e-1   1.0e' + tolstr +'\r\n')
    return f

def write_gener(f, eleme, phase = 'brine', mass_rate = .1585, \
        kg_inflow = [], times = [], co2_enthalpy = 50.,\
        isothermal = True):
    """
        Writes GENER file
            So far, has only been used for the case when the 
            kg_inflow list is not given, thus implying a 
            constant injection rate
        mass_rate: injection rate of component in kg/s
        phase: either 'co2' or 'brine'
        co2_enthalpy: for nonisothermal cases
        isothermal: True or False
        kg_inflow: list of mass inflow rates ending at the corresponding
            index in times
        times: list of times that end variable injection periods.
    """
    write_separator(f,'GENER')
    if phase =='co2':
        p = 'COM3 '
    else:
        p = 'COM1 '

    if kg_inflow == [] :
        mr = '{: 10.4f}'.format(mass_rate)
        heat_rate = '{: 10.4f}'.format(mass_rate * 10.e5)
        if isothermal == True:
            se = ''
        else:
            se = '{: 10.3e}'.format(co2_enthalpy)
        f.write(eleme + 'inj 1'+ 19*' ' + '1' + 5 * ' ' + p + mr + se + '\r\n')
        #f.write(eleme + 'hot 1'+ 19*' ' + '1' + 5 * ' ' + \
                #'HEAT ' + heat_rate + '\r\n')
    else: 
        if len(kg_inflow) != len(times):
            print "mass rate and time arrays are not the same length"
            return 1

        # store the formatted strings in lists
        massprint = []
        timeprint = []
        for i in range(len(kg_inflow)):
            b = '{: 15.3e}'.format(times[i] * 24. * 3600.)
            l = list(b)
            del l[-2]
            if l[-2] == '+':
                l[-2] = '0'
            c = ''.join(l)
            timeprint.append(' ' + c)
            
            # quick and dirty formatting of mass inflow. 
            # only works if the mass flow rate is between -1 and +1. 
            b = '{: 15.3e}'.format(kg_inflow[i])
            l = list(b)
            del l[-2]
            if l[-2] == '+':
                l[-2] = '0'
            c = ''.join(l)
            massprint.append(' ' + c)

        # LTAB parameter
        ltab = '{:5d}'.format(len(kg_inflow))

        f.write(eleme + 'inj 1'+ 10 * ' '  + 5 * ' ' \
                + ltab + 5 * ' ' + p   + '\r\n')

        # write F1: Generation times
        for i in range(1, len(timeprint)+1):
            f.write(timeprint[i-1])
            if i % 5 == 0:
                f.write('\r\n')
        f.write('\r\n')

        #write F2: Generation rates
        for i in range(1, len(massprint)+1):
            f.write(massprint[i-1])
            if i % 5 == 0:
                f.write('\r\n')
        f.write('\r\n')

        # write F3: specific enthalpy
        specific_enthalpy = '{: 16.1e}'.format(co2_enthalpy)
        se = list(specific_enthalpy)
        del se[-3]
        specific_enthalpy = ''.join(se)
        for i in range(1, len(timeprint)+1):
            f.write(specific_enthalpy)
            if i % 5 == 0:
                f.write('\r\n')
        f.write('\r\n')
    f.write('\r\n')
    return f

def write_times(f, timelist):
    numtimes = len(timelist)

    write_separator(f,'TIMES')
    f.write('   '+ str(numtimes)+'\r\n')

    count = 0
    for i in range(numtimes):
        seconds = timelist[i] * 24. * 3600.
        b = '{: 10.4e}'.format(seconds)
        l = list(b)
        del l[-2]
        del l[-2]
        c = ''.join(l)
        f.write(' ' +  c)
        count +=1
        if count == 8: 
            f.write('\r\n')
            count = 0
    f.write('\r\n')
    return f

def write_foft(f):
    write_separator(f, 'FOFT ')
    # place element ID's here to get some output
    f.write('\r\n')
    return f

def write_coft(f, sleipner):
    """
        Gets output from a specific connection
    """
    write_separator(f, 'COFT ')
    if sleipner == True:
        #f.write('bH732cH732\r\n')
        dummy =1
    else:
        #f.write('aB212bB212\r\n')
        f.write('xB212yB212\r\n')
    f.write('\r\n')
    return f

def write_goft(f):
    """
        What is the generation rate output?
    """
    write_separator(f, 'GOFT ')
    f.write('\r\n')
    return f

def write_meshmaker(f , rect = True, flat = True, nx = 65, ny = 119, nz = 1):
    """
        uses internal MESHMAKER routine.
    """
    write_separator(f, 'MESHMAKER')
    f.write('XYZ\r\n')
    # degrees to rotate domain
    f.write('        0.\r\n')

    f.write('NX     '+ str(nx)+'\r\n')
    count = 0
    for i in range(nx):
        f.write('       50.')
        count += 1
        if count % 8 == 0:
            f.write('\r\n')
    f.write('\r\n')

    f.write('NY     '+ str(ny) + '\r\n')
    count = 0
    for i in range(ny):
        f.write('       50.')
        count += 1
        if count % 8 == 0:
            f.write('\r\n')
    f.write('\r\n')
    if flat == True:
        f.write('NZ       1      28.0\r\n') 
    else:
        f.write('NZ     '+ str(nz) + '\r\n')
        count = 0
        for i in range(nz):
            f.write('       50.')
            count +=1
            if count % 8 == 0:
                f.write('\r\n')

    f.write('\r\n')
    f.write('\r\n') 
    return f

def format_float_mesh(val):
    """ takes in a double, returns a string to be entered into the MESH file
            ( output 's' , val )
            ('-.5000E-00', -0.5)
            ('0.1000E+01', 1.0)
            ('0.3500E+03', 350.0)
            ('-.7500E-03', -0.00075)
            ('0.4100E-01', 0.041)
    """
    a = '{: .3E}'.format(val)
    l = list(a)
    l[1], l[2] = l[2], l[1]
    if l[-3] == '+':
        l[-1] = str( int(l[-1]) + 1 )
    elif l[-3] == '-':
        l[-1] = str( int(l[-1]) - 1 )
    if l[0] != '-':
        l[0] = '0'
    s = ''.join(l)
    return s 

def format_float_rocks(val):
    """
        formats floats to fit well into the ROCKS block.
    """
    return '{: .3e}'.format(val)

def format_float_incon(val):
    """ formats to fit in the INCON, includes blank chars before number
    only works for positive floating point values (I think)
    """
    a = '{:> 20.12E}'.format(val)
    l = list(a)
    l[2], l[3] = l[3], l[2]
    if l[-3] == '+':
        l[-1] = str( int(l[-1]) + 1 )
    elif l[-3] == '-':
        l[-1] = str( int(l[-1]) - 1 )
    else : 
        print "What you input was not a float or somthin"
        print "writingFunctions.format_float_mesh(val)"
        print "val = " + str(val)
    s = ''.join(l)
    return s

class T2InputGrid(object):
    """
        T2InputGrid contains the functions that handle the numerical grid
        while creating TOUGH2/ECO2N input files. 

         
        self.boundary = []
        - An ancillary list that keeps track of cells on the boundary 
          of the domain, used for making fixed-pressure boundary conditions
          during injection simulations

        self.el_array = []
        - 3d array containing all of the element keys in space. When called:
          self.el_array[i][j][k]
          returns the 5 character element key.
        self.x = {}
        - Dictionary that returns the x position of an element based on its key:

        self.y = {}
        - Dictionary that returns the y position of an element 
        
        self.z = {}
        - Dictionary that returns the z position of an element 

        self.corners = {}
        - Dictionary that contains a sorted list of corners from the eclipse
          grid file. A corner contains the x-y-z location of a polyhedron
          corner from an ECLIPSE file. 
          Each list has length of 8. The indices are labeled below:
          ----- t2_corner_sketch ------------
          High elevation    [0] *-----* [1]    ^ 
          plan view             |     |        | North
                                |     |
                            [2] *-----* [3]

          Low elevation     [4] *-----* [5]    ^ 
          plan view             |     |        | North
                                |     |
                            [6] *-----* [7]

        self.mat = {}     
        - Dictionary that returns the rock type of an element 

        self.vol = {}
        - Dictionary that returns the volume  of an element 

        self.ahtx = {}
        - Dictionary that returns the interface area for heat exchange 
          between semi infinite confining beds.
        - this parameter is auto-generated by the MESHMAKER routine internally
          in TOUGH2 and is generated in the same way here.
          The MESHMAKER routine, for a radial injection problem, uses the 
          relationship vol/ahtx = 50. 
          The same relationship is specified in the generation routine below.
        
        self.pres = {}
        - Dictionary that returns the fluid pressure of an element 

        self.na_cl = {}
        - Dictionary that returns the salt concentration of an element 

        self.x_co2 = {}
        - Dictionary that returns the dissolved CO2 concentration of an element

        self.temp = {}
        - Dictionary that returns the temperature of an element 

        self.num_elements = 0
         - Number of elements in domain. 
         To be used as a check after elements are created

    """
    def __init__(self, nx, ny, nz):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.elements = []
        self.boundary = []
        self.el_array = []
        self.x = {}
        self.y = {}
        self.z = {}
        self.corners = {}
        self.mat = {}
        self.vol = {}
        self.ahtx = {}
        self.pres = {}
        self.na_cl = {}
        self.x_co2 = {}
        self.temp = {}
        self.num_elements = 0

    class Corner(object):
        def __init__(self, x, y, z):
            self.x = x
            self.y = y
            self.z = z

        def get_x(self):
            return self.x

        def get_y(self):
            return self.y

        def get_z(self):
            return self.z


    def get_x_char(self, i):
        """ returns a 2-character string to be added to the element ID
            self.get_x_char(4)  -> ' 4'
            self.get_x_char(17) -> '17'
        """
        if i > 99 or i < 0 :
            print "only works for nx < 100"
            return 1
        else : 
            c1 = str(i)
            if len(c1) == 1:
                char = ' ' + c1
            elif len(c1) == 2:
                char = c1
            else:
                print "I'm pretty sure you probably put an index less than 1"
            return char

    def get_y_char(self, j):
        """ returns a 2-character string to be added to the element ID
            self.get_y_char(4)  -> 'A4'
            self.get_y_char(17) -> 'B7'
        """
        if j > 259 or j < 0:
            print "only works for ny < 260"
            return 1
        else : 
            letters = string.uppercase
            c1 = letters[ j/10]
            c2 = str( j - 10 * (j/10))
            char = c1 + c2
            return char

    def get_z_char(self, k):
        """ returns a 1-character string to be added to the element ID
            self.get_z_char(4)  -> 'e'
            self.get_z_char(17) -> 'r'
        """
        if k > 51 or k < 0:
            print 'Only works for nz < 52'
            return 1
        else:
            letters =  string.lowercase + string.uppercase
            char = letters[k]
            return char

    def get_element_chars(self, i, j, k):
        """ Returns a 5 character element ID for i,j,k
            get_element_chars(1, 2, 3)    -> 'dA2 1'
            get_element_chars(25, 65, 11) -> 'lG525'

            The scheme here was used to resemble the scheme used in the ECO2N
            manual but different schemes could be used to generate 
            element IDs.
            Here, if the element id is '12345', 
            '1' contains the z axis location
            '23' contains the y axis location and 
            '45' contains the x axis location.

        """
        c1 = self.get_x_char(i)
        c2 = self.get_y_char(j)
        c3 = self.get_z_char(k)
        el = c3 + c2 + c1
        if len(el) != 5:
            print "i, j, k"
            print i, j, k 
            print "c1, c2, c3"
            print c1, c2, c3
            print "failed to give correct character length"
            return 1
        else:
            return el

    def fill_uniform_grid(self, porosity, dx, dy, dz, density = 1000.,\
            solubility_limit = 0.454104e-3, inj_cell = 'none',\
            altered_cell = 'none', temperature = 32.):
        # populates grid 
        g = open('MESH', 'w')
        print "MESH created with:"
        for i in range(self.nx):
            temparray = []
            for j in range(self.ny):
                temparrayz = []
                for k in range(self.nz):
                    eleme = self.get_element_chars(i, j, k)
                    self.elements.append(eleme)

                    xlocal = (dx / 2 + dx * i )
                    ylocal = (dy / 2 + dy * j ) 
                    zlocal = -(dz / 2 + dz * k ) - 859.5
                    if altered_cell != 'none':
                        if eleme == altered_cell:
                            zlocal = zlocal + 5.
                    self.x[eleme] = xlocal
                    self.y[eleme] = ylocal
                    self.z[eleme] = zlocal
                    self.vol[eleme] = (dz * dy * dx)
                    # this thermal contact area was specified based on the 
                    # MESHMAKER routine in ECO2N
                    self.ahtx[eleme] = (dz * dy * dx) / 50.
                    self.pres[eleme] = -zlocal * density * 10.0
                    self.na_cl[eleme] = 3.2e-2
                    self.x_co2[eleme] = solubility_limit
                    self.temp[eleme] = temperature
                    temparrayz.append(eleme)

                    if eleme == inj_cell:
                        self.mat[eleme] = 'sands'
                        self.mat[eleme] = 'well '
                        print "INJECTION cell"
                        print eleme
                    else:
                        self.mat[eleme] = 'sands'

                    vw = format_float_mesh(self.vol[eleme])
                    aw = format_float_mesh(self.ahtx[eleme])
                    xw = format_float_mesh(self.x[eleme])
                    yw = format_float_mesh(self.y[eleme])
                    zw = format_float_mesh(self.z[eleme])
                temparray.append(temparrayz)
            self.el_array.append(temparray)

        return 0


    def e_cell_index(self,i,j,k):
        """ returns the index used to call the cell list"""
        nx = 65
        ny = 119
        ind = i + nx * j + nx * ny * k
        return ind 

    def fill_eclipse_grid(self, e_cells , temperature = 37., density = 1000.,\
            two_d = False,  gradient = 10 , solubility = 0.474e-1, \
            shale = True, inj_cell = 'JH732'):
        """Fills the 3d grid with x, y, and z
        'gradient' specifies the hydrostatic gradient. [MPa/km]

        'solubility' specifies the initial dissolved CO2 in aqueous phase.

        'two_d', 
        if true, ensures that the grid being created is two-dimensional.

        'thickness', if specified, returns a constant thickness, also requires 
        'depth' to be specified in order to get the correct initial pressure.
        'gradient' is the pressure gradient, given in units of MPa/km

        if thickness and depth are not specified, 
        the sleipner vertical geometry is considered.
        """
        print "FILLING 3D grid based on Sleipner input data........"

        if two_d == True:
            self.nz = 1

        count = 0
        for i in range(0, self.nx):
            temparray = []
            for j in range(self.ny-1, -1, -1):
                temparrayz = []
                sand_count = 0
                for k in range(0, self.nz):
                    ind = self.e_cell_index(i, j, k)
                    eleme = self.get_element_chars(i, j, k)
                    if eleme == inj_cell:
                        injector = True
                    else: 
                        injector = False
                    if shale == True or e_cells[ind].getXPermeability() > 1.:
                        if e_cells[ind].getXPermeability() > 1.:
                            sand_count +=1
                        count +=1
                        self.elements.append(eleme)
                        self.write_eclipse_cell(e_cells, i, j, k,\
                                temperature, density, injector ,\
                                two_d, gradient, solubility)
                        temparrayz.append(eleme)
                temparray.append(temparrayz)
            self.el_array.append(temparray)

        # quick debugging to make sure all elements were input
        self.num_elements = count
        check_element_length = len(self.el_array) * len(self.el_array[0]) * \
                len(self.el_array[0][0])
        if shale == True:
            input_length = self.nx * self.ny * self.nz
        else: 
            input_length = self.nx * self.ny * 34
        assert self.num_elements == check_element_length, \
                "num_elements created does not match array dimensions"
        assert self.num_elements == input_length, \
                "num_elements created does not match input nx*ny*nz"

        print "GRID FILLING COMPLETE"
        print str(count) + " ELEMENTS CREATED"
        print "__________________________"
        return 0

    def write_eclipse_cell(self, e_cells, i, j, k,\
            temperature = 37., density = 1000., injector = False,\
            two_d = False,  gradient = 10 , solubility = 0.474e-1):
        """ converts the corners from the eclipse_cell order into the TOUGH2 input
            order for convenience in doing volume calculations
        """

        eclipse_index = self.e_cell_index(i, j, k)
        eleme = self.get_element_chars(i, j, k)

        if e_cells[eclipse_index].getXPermeability() > 1 or two_d == True:
            self.mat[eleme] = 'sands'
        elif injector == True:
            self.mat[eleme] = 'well '
            #self.mat[eleme] = 'sands'
        else:
            self.mat[eleme] = 'shale'


        # takes original corners from eclipse_cell and converts them to 
        # local corners to be used in volume calcs
        # the goal is to get the x,y,z coordinates of the centroid and get 
        # the cell volume to 
        if two_d == True:
            # Finds uppermost grid cell and lowermost grid cell
            # from eclipse cells
            pc_count = []
            for k in range(self.nz):
                pc_ind = self.e_cell_index(i,j,k)
                if e_cells[pc_ind].getXPermeability() > 1.:
                    pc_ind_1 = self.e_cell_index(i,j,k+1)
                    pc_count.append(k)
            bot_ind = self.e_cell_index(i, j, pc_count[-1])
            top_ind = self.e_cell_index(i, j, pc_count[0])

            # takes high elevation corners from uppermost cell and low elevation
            # corners from lowermost cell and makes vertically averaged cell
            top_corners_in = e_cells[top_ind].getCorners()
            bot_corners_in = e_cells[bot_ind].getCorners()
            vert_avg_corners = []
            for corner_index in range(4):
                x, y = top_corners_in[corner_index].getXY()
                z = -top_corners_in[corner_index].getZ()
                nc = self.Corner(x, y, z)
                vert_avg_corners.append(nc)
            for corner_index in range(4,8):
                x,y = bot_corners_in[corner_index].getXY()
                z = -bot_corners_in[corner_index].getZ()
                nc = self.Corner(x,y,z)
                vert_avg_corners.append(nc)
            temp_new_corners = vert_avg_corners
        else:
            # takes the regular 3d, local scale corners instead
            original_corners = e_cells[eclipse_index].getCorners()
            temp_new_corners = []
            for c in original_corners:
                x, y = c.getXY()
                # Depth from eclipse file is converted to elevation
                z = - c.getZ()
                nc = self.Corner(x, y, z)
                temp_new_corners.append(nc)

        self.corners[eleme] = temp_new_corners
        x = self.get_x_centroid(temp_new_corners)
        y = self.get_y_centroid(temp_new_corners)
        z = self.get_z_centroid(temp_new_corners)
        self.x[eleme] = x
        self.y[eleme] = y
        self.z[eleme] = z
        self.vol[eleme] = self.get_volume(x, y, z, temp_new_corners)

        pressure = -z * gradient * density
        self.pres[eleme] = pressure

        # This thermal contact area is specified based on the internal
        # generation routines in the MESHMAKER capabilities of TOUGH2
        self.ahtx[eleme] = self.vol[eleme] / 50.
        self.na_cl[eleme] = 3.2e-2
        self.x_co2[eleme] = solubility
        self.temp[eleme] = temperature

    def get_x_centroid(self, corners):
        count = 0.
        sum_c = 0.
        for c in corners:
            count += 1.
            sum_c += c.get_x()
        return sum_c / count

    def get_y_centroid(self, corners):
        count = 0.
        sum_c = 0.
        for c in corners:
            count += 1.
            sum_c += c.get_y()
        return sum_c / count

    def get_z_centroid(self, corners):
        count = 0.
        sum_c = 0.
        for c in corners:
            count += 1.
            sum_c += c.get_z()
        return sum_c / count
  
    def get_volume(self, x, y, z, corners):
        """ uses the equation for volume of an orientable polyhedron
            V = 1/3 \sum_i x_i \dot n^hat_i A_i
            See 
            http://en.wikipedia.org/wiki/Polyhedron#Volume
            for details.
        """
        face_map = ['west', 'south', 'east', 'north', 'bot', 'top']

        v_sum = 0.0
        for face in face_map:
            a = self.get_area(corners, face)
            centroid = self.get_face_center(x, y, z, corners, face)
            cent = np.asarray(centroid)
            vec = self.get_normal_vector(x, y, z, corners, face)
            v_sum += np.dot(cent, vec) * a

        vol = 1./3. * v_sum
        return vol

    def get_area(self, corners, face):
        """ returns the area of a cell face, east, west, etc
          For the vertical faces (east, west, north, south), this diagram is 
          used to find the area by the trapezoid rule. 

          -------------------------------------------------------
          cross section     [x1,y1] *-----* [x2,y2]    ^ 
                                    |     |        | positive Z
                                    |     |
                            [x1,y3] *-----* [x2,y4]

          Note that the upper 'x' values are identical to lower ones since 
          the eclipse cells are essentially columns.

          For the horizontal face (top and bot), a planar fit is performed
          to find the function z = f(x,y) in order to do a surface integral
          in the manner seen here: 
          http://en.wikipedia.org/wiki/Surface_integral

          See the sketch below 't2_corner_sketch' in the T2InputGrid 
          description for a diagram of all corners in the grid cell.
        """ 
        if face == 'west':
            x1 = corners[2].get_y()
            x2 = corners[0].get_y()
            y1 = corners[2].get_z()
            y2 = corners[0].get_z()
            y3 = corners[6].get_z()
            y4 = corners[4].get_z()
            area = -self.get_area_side(x1, x2, y1, y2, y3, y4)
        elif face == 'south':
            x1 = corners[2].get_x()
            x2 = corners[3].get_x()
            y1 = corners[2].get_z()
            y2 = corners[3].get_z()
            y3 = corners[6].get_z()
            y4 = corners[7].get_z()
            area = -self.get_area_side(x1, x2, y1, y2, y3, y4)
        elif face == 'east':
            x1 = corners[3].get_y()
            x2 = corners[1].get_y()
            y1 = corners[3].get_z()
            y2 = corners[1].get_z()
            y3 = corners[7].get_z()
            y4 = corners[5].get_z()
            area = -self.get_area_side(x1, x2, y1, y2, y3, y4)
        elif face == 'north':
            x1 = corners[0].get_x()
            x2 = corners[1].get_x()
            y1 = corners[0].get_z()
            y2 = corners[1].get_z()
            y3 = corners[4].get_z()
            y4 = corners[5].get_z()
            area = -self.get_area_side(x1, x2, y1, y2, y3, y4)
        elif face == 'bot':
            nc = [corners[6], corners[7], corners[4], corners[5]]
            c, resid, rank, sigma = self.fit_plane(nc)
            mag = np.sqrt(pow(c[0],2.) + pow(c[1],2.) + 1)
            x1 = corners[2].get_x()
            x2 = corners[3].get_x()
            y1 = corners[2].get_y()
            y2 = corners[0].get_y()
            area = mag * ((x2 * y2 - x1 * y2) - (x2 * y1 - x1 * y1))
        elif face == 'top':
            nc = [corners[2], corners[3], corners[0], corners[1]]
            c, resid, rank, sigma = self.fit_plane(nc)
            mag = np.sqrt(pow(c[0],2.) + pow(c[1],2.) + 1)
            x1 = corners[6].get_x()
            x2 = corners[7].get_x()
            y1 = corners[6].get_y()
            y2 = corners[4].get_y()
            area = mag * ((x2 * y2 - x1 * y2) - (x2 * y1 - x1 * y1))
        else:
            raise Exception("Invalid Face, please specify" +  \
                    "one of the six faces in face_map \n\n")
        return area

    def get_area_side(self, x1, x2, y1, y2, y3, y4):
        """
            gets the area of a polygon using the trapezoid rule
        """
        h = x2 - x1
        b1 = y4 - y2
        b2 = y3 - y1
        return 0.5 * h * (b1 + b2)

    def get_face_center(self, xc, yc, zc, corners, face):
        """ center vector location relative to polyhedron center
        """
        if face == 'west':
            nc = [corners[0], corners[2], corners[4], corners[6]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
        elif face == 'south':
            nc = [corners[2], corners[3], corners[6], corners[7]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
            a = 2
        elif face == 'east':
            nc = [corners[3], corners[1], corners[7], corners[5]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
        elif face == 'north':
            nc = [corners[0], corners[1], corners[4], corners[5]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
        elif face == 'bot':
            nc = [corners[6], corners[7], corners[4], corners[5]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
        elif face == 'top':
            nc = [corners[2], corners[3], corners[0], corners[1]]
            xf = self.get_x_centroid(nc)
            yf = self.get_y_centroid(nc)
            zf = self.get_z_centroid(nc)
        else:
            raise Exception("Invalid Face, please specify" +  \
                    "one of the six faces in face_map \n\n")

        vec = [xf - xc, yf - yc, zf - zc]
        return vec

    def get_normal_vector(self, x, y, z, corners, face):
        """ gets normal vector of face
        """
        if face == 'west':
            vec = [-1., 0., 0.]
        elif face == 'south':
            vec = [0., -1., 0.]
        elif face == 'east':
            vec = [1., 0., 0.]
        elif face == 'north':
            vec = [0., 1., 0.]
        elif face == 'bot':
            nc = [corners[6], corners[7], corners[4], corners[5]]
            c, resid, rank, sigma = self.fit_plane(nc)
            mag = np.sqrt(pow(c[0], 2.) + pow(c[1],2.) + 1)
            vec = [c[0]/mag, c[1]/mag, -1./mag]
        elif face == 'top':
            nc = [corners[2], corners[3], corners[0], corners[1]]
            c, resid, rank, sigma = self.fit_plane(nc)
            mag = np.sqrt(pow(c[0], 2.) + pow(c[1],2.) + 1)
            vec = [-c[0]/mag, -c[1]/mag, 1./mag]
        else:
            raise Exception("Invalid Face, please specify" +  \
                    "one of the six faces in face_map \n\n")
        return vec

    def fit_plane(self, corners):
        """ takes four corner points and fits a plane least squares to them
            returns in form z = c[0] x + c[1] y + c[2]
        """
        x = []
        y = []
        z = []
        for c in corners:
            x.append(c.get_x())
            y.append(c.get_y())
            z.append(c.get_z())
        x = np.asarray(x)
        y = np.asarray(y)
        z = np.asarray(z)
        A = np.column_stack((x, y, np.ones(x.size)))
        c, resid, rank, sigma = np.linalg.lstsq(A, z)
        return c, resid, rank, sigma

    def write_mesh(self, e_cel, two_d = False, uniform = False,\
            boundary_type = 1, shale = True,\
            type1_source_cell = 'none'):
        """ populates grid and writes ELEME block of MESH
        """
        g = open('MESH', 'w')
        print "MESH created with:"
        print "writing ELEMENT block of input data"
        g.write('ELEME\r\n')
        if shale == False:
            nzi = 34
        else:
            nzi = self.nz
        count = 0 
        for i in range(0, self.nx):
            for j in range(0, self.ny):
                for k in range(0, nzi):
                    eleme = self.el_array[i][j][k]
                    # if the cell is interior, writes it.
                    # interior check is below
                    # if not, adds it to the boundary group to be 
                    # written at the end
                    if boundary_type == 1 and \
                            (i == 0 or i == (self.nx - 1 ) or \
                            j == 0 or j == (self.ny - 1 )):
                        self.boundary.append(eleme)
                    elif boundary_type == 1 and type1_source_cell == eleme:
                        print "YAY"
                        print "element ", eleme, " is on the boundary!~"
                        print "YAY"
                        self.boundary.append(eleme)
                    else: 
                        vw = format_float_mesh(self.vol[eleme])
                        aw = format_float_mesh(self.ahtx[eleme])
                        xw = format_float_mesh(self.x[eleme])
                        yw = format_float_mesh(self.y[eleme])
                        zw = format_float_mesh(self.z[eleme])
                        mat = self.mat[eleme]
                        g.write(eleme + 5 * ' ' + 5 * ' ' + mat + vw + aw + \
                            10 * ' ' + xw + yw + zw + '\r\n')
                        count +=1

        if boundary_type == 1:
            g.write('ina\r\n')

        for el in self.boundary:
            vw = format_float_mesh(self.vol[el])
            aw = format_float_mesh(self.ahtx[el])
            xw = format_float_mesh(self.x[el])
            yw = format_float_mesh(self.y[el])
            zw = format_float_mesh(self.z[el])
            mat = self.mat[el]
            g.write(el + 5 * ' ' + 5 * ' ' + mat + vw + aw + \
                10 * ' ' + xw + yw + zw + '\r\n')
            count += 1 
        g.write('\r\n')

        print "ELEME: " + str(count) + " elements"
        print "ELEMENTS COMPLETE --------------- "
        print "Writing Connections.........."

        g.write('CONNE\r\n')
        count = 0 
        for i in range(0, self.nx ):
            for j in range(0, self.ny ):
                for k in range(0, nzi):
                    # z connection
                    if k != nzi -1:
                        g, count = self.write_z_connection(g, count, \
                                i, j, k, uniform)
                    # y connection
                    if j != self.ny -1:
                        g, count = self.write_y_connection(g, count, \
                                i, j, k, two_d, uniform, e_cel )
                    # x connection
                    if i != self.nx -1: 
                        g, count = self.write_x_connection(g, count, \
                                i, j, k, two_d, uniform, e_cel )

        g.write('\r\n')
        g.write('\r\n')
        
        g.close()
        print "CONNE: " + str(count) +  " connections"
        print "CONNECTIONS COMPLETE -----------------"
        return 0

    def write_z_connection(self, g, count, i, j, k, uniform = False):
        """ should be legit now. 
        """
        if uniform == False:
            direc = 3  
            el1 = self.el_array[i][j][k]
            el2 = self.el_array[i][j][k+1]
            z1 = self.z[el1]
            z2 = self.z[el2]
            area = self.get_area(self.corners[el1], 'bot')
            mid = self.get_face_center(self.x[el1], self.y[el1], self.z[el1],\
                    self.corners[el1], 'bot')
            zmid = z1 + mid[2]
            dz1 = z1 - zmid
            dz2 = zmid - z2
            dz1w = format_float_mesh(dz1)
            dz2w = format_float_mesh(dz2)
            aw = format_float_mesh(area)
            beta = 1.0
            betaw = format_float_mesh(beta) 
            z_string = el1 + el2 + 19 * ' ' + str(direc) + \
                dz1w + dz2w + aw + betaw + '\r\n'
        else:
            direc = 3  
            el1 =  self.el_array[i][j][k]
            el2 =  self.el_array[i][j][k+1]
            z1 = self.z[el1]
            z2 = self.z[el2]
            dzn = (z1 - z2)/2.
            dz1 = format_float_mesh((z1 - 0.5 * ( z1 + z2)))
            dz2 = format_float_mesh(-(z2 - 0.5 * (z1 + z2)))
            a1 = self.vol[el1] / dzn
            a2 = self.vol[el2] / dzn
            a = min(a1, a2)
            dzw = format_float_mesh(dzn) 
            aw = format_float_mesh(a)
            beta = 1.0
            betaw = format_float_mesh(beta) 
            z_string = el1 + el2 + 19 * ' ' + str(direc) + \
                            dz1 + dz2 + aw + betaw + '\r\n'
        g.write(z_string)
        count +=1
        return g, count

    def write_y_connection(self, g, count, i, j, k, \
            two_d = False, uniform = False, e_cells = 'uniform' ):
        """ need to incorporate z, inside this loop
        """ 
        direc = 2  
        el1 =  self.el_array[i][j][k]
        el2 = self.el_array[i][j+1][k]
        if uniform == True or two_d == True:
            y1 = self.y[el1]
            y2 = self.y[el2]
            dyn = (y2 - y1)/2.
            a1 = self.vol[el1] / dyn
            a2 = self.vol[el2] / dyn
            a = min(a1, a2)
            dyw = format_float_mesh(dyn) 
            aw = format_float_mesh(a)
            z1 = self.z[el1]
            z2 = self.z[el2]
            dz = z1 - z2
            beta = dz / np.sqrt(pow(dz,2) + pow(dyn,2))
            betaw = format_float_mesh(beta) 
            y_string = el1 + el2 + 19 * ' ' + str(direc) + \
                dyw + dyw + aw + betaw + '\r\n'
            g.write(y_string)
            count +=1
        else:
            area = self.get_area(self.corners[el1], 'south')
            dy1 = self.y[el1] - self.corners[el1][3].get_y()
            dy2 = self.corners[el2][1].get_y() - self.y[el2] 
            dy1w = format_float_mesh(dy1)
            dy2w = format_float_mesh(dy2)
            aw = format_float_mesh(area)
            dz = self.z[el1] - self.z[el2]
            beta = dz / np.sqrt(pow(dz,2) + pow(dy1 + dy2,2))
            betaw = format_float_mesh(beta)
            y_string = el1 + el2 + 19 * ' ' + str(direc) + \
                dy1w + dy2w + aw + betaw + '\r\n'
            g.write(y_string)
            count +=1
            
        return g, count

    def write_x_connection(self, g, count, i, j, k, \
            two_d = True, uniform = False, e_cells = 'uniform' ):
        """ need to incorporate z, note, inside this loop
        """
        direc = 1
        el1 =  self.el_array[i][j][k]
        el2 =  self.el_array[i+1][j][k]
        if uniform == True or two_d == True:
            x1 = self.x[el1]
            z1 = self.z[el1]
            x2 = self.x[el2]
            dxn = (x2 - x1)/2.
            a1 = self.vol[el1] / dxn
            a2 = self.vol[el2] / dxn
            a = min(a1, a2)
            dxw = format_float_mesh(dxn) 
            aw = format_float_mesh(a)
            z2 = self.z[el2]
            dz = z1 - z2
            beta = dz / np.sqrt(pow(dz,2) + pow(dxn,2))
            betaw = format_float_mesh(beta) 
            x_string = el1 + el2 + 19 * ' ' + str(direc) + \
                dxw + dxw + aw + betaw + '\r\n'
            g.write(x_string)
            count +=1
        else:
            area = self.get_area(self.corners[el1], 'east')
            dx1 = self.corners[el1][3].get_x() - self.x[el1]
            dx2 = self.x[el2] - self.corners[el2][2].get_x()
            dx1w = format_float_mesh(dx1)
            dx2w = format_float_mesh(dx2)
            aw = format_float_mesh(area)
            dz = self.z[el1] - self.z[el2]
            beta = dz / np.sqrt(pow(dz,2) + pow(dx1 + dx2,2))
            betaw = format_float_mesh(beta)
            x_string = el1 + el2 + 19 * ' ' + str(direc) + \
                dx1w + dx2w + aw + betaw + '\r\n'
            g.write(x_string)
            count +=1
             
        return g, count

    def write_incon(self, shale = True, porosity = []):
        print "Writing INCON FILE"
        h = open('INCON','w')
        h.write('INCON -- INITIAL CONDITIONS FOR ' + str(len(self.elements)) + \
            ' ELEMENTS AT TIME  .000000E+00\r\n')
        if shale == False:
            nzi = 34
        else:
            nzi = self.nz
        for i in range(0, self.nx):
            for j in range(0, self.ny):
                for k in range(nzi):
                    el = self.el_array[i][j][k]
                    if porosity != []:
                        poro = '{: 15.8E}'.format(porosity)
                    else:
                        poro = '{: 15.8E}'.format(self.poro[el])
                    h.write(el + 10 * ' '+ poro + '\r\n')
                    x1 =  format_float_incon(self.pres[el])
                    x2 =  format_float_incon(self.na_cl[el])
                    x3 =  format_float_incon(self.x_co2[el])
                    x4 =  format_float_incon(self.temp[el])
                    h.write( x1 + x2 + x3 + x4 + '\r\n')
        h.write('\r\n')
        h.write('\r\n')
        h.close()
        print "INCON COMPLETE"
        return 0 

    def use_old_incon(self, hydro_directory, type1_source_cell = 'none',\
            saturation_fraction= 0.8):
        """ note that the hydro_directory does not require the '_dir/' suffix
        """
        print "Writing INCON with SAVE file from: " + hydro_directory
        with open(hydro_directory + '_dir/' + 'SAVE','r') as f:
            content = f.readlines()
        content[-2] = '\n'
        content[-1] = '\n'
        del content[0]
        f = open('INCON','w')
        f.write('INCON -- INITIAL CONDITIONS FOR ' + str(len(self.elements)) + \
            ' ELEMENTS AT TIME  .000000E+00\n')
        saturate = False
        for line in content:
            s = line.split()
            if saturate == True:
                print "saturating"
                print line
                gas_sat = 10 +  saturation_fraction
                brine_sat = 0.1 *(1. - saturation_fraction)
                brine_res = format_float_incon(brine_sat)
                gas_sat = format_float_incon(gas_sat)
                pres = format_float_incon(float(s[0]))
                temp = format_float_incon(float(s[3]))
                line = pres + brine_res + gas_sat + temp + '\r\n'
                print "changed line"
                print line
                saturate = False
            if s != [] and len(s[0]) == 3:
                s[0] = s[0] + '_' + s[1]
            if s != [] and type1_source_cell == s[0]:
                print "found one!"
                print line
                print type1_source_cell, s[0]
                saturate = True
            f.write(line)
        f.close()
        print "INCON from SAVE complete"
        return 0
    
    # grid Plotting routines ------------------------------------------------------
    def plot_cells(self, show = True, two_d = False):
        print " Creating scatter plot of all cells."
        x = []
        y = []
        z = []
        pres = []
        for el in self.elements:
            x.append(self.x[el])
            y.append(self.y[el])
            z.append(self.z[el])
            pres.append(self.pres[el])

        xs = np.asarray(x)
        ys = np.asarray(y)
        zs = np.asarray(z)
        ps = np.asarray(pres)
        if two_d == True:
            xlist = []
            ylist = []
            zlist = []
            tempzlist = []
            plist = []
            tempplist = []
            count = 0
            for n in range(len(z)):

                tempzlist.append(z[n])
                tempplist.append(pres[n])
                if n < 119:
                    ylist.append(y[n])
                if n%119 == 0:
                    xlist.append(x[n])
                count +=1
                if count%119 == 0:
                    zlist.append(tempzlist)
                    plist.append(tempplist)
                    tempzlist = []
                    tempplist = []

            xg, yg = np.meshgrid(xlist, ylist)
            zg = np.asarray(zlist)
            zg = np.transpose(zg)
            pg = np.asarray(plist)
            pg = np.transpose(pg)


        if two_d == True:
            fig = plt.figure(figsize=(7.5,10), dpi = 480, facecolor = 'w',\
                edgecolor ='k')
            ax = fig.add_subplot(111)
            CS = ax.contourf(xg, yg, pg)
            CB = plt.colorbar(CS, shrink = 0.8, extend = 'both')
            title = 'contour_cells'
        else:
            fig = plt.figure(figsize=(11.,8.5), dpi = 960, facecolor = 'w',\
                edgecolor ='k')
            ax = fig.add_subplot(111, projection = '3d')
            p1 = ax.scatter(xs,ys,zs, s=10, facecolor = (.7,.7,.7), c=zs)
            plt.gray()
            ax.set_zlabel('Z')
            title = 'scatter_cells'

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        if show == True:
            plt.show()
        plt.savefig(title + '.pdf', format = 'pdf')
        plt.clf()
        plt.close()
        print "scatter plot complete"
        return 0

    def plot_cells_slice(self, direc = 1, ind = 0, show = True):
        """direc is the direction, ind is the index of the other direction
        for example, if direc == 1 and ind == 4,
        this routine will plot the cells in the x-z plane with y index of 4
        """
        print "plotting slice"
        x_mod = []
        z_mod = []
        for el in self.elements:
            if direc == 1:
                chars = self.get_y_char(ind)
                if el[1:3] == chars:
                    x_mod.append(self.x[el])
                    z_mod.append(self.z[el])
            elif direc == 2:
                chars = self.get_x_char(ind)
                if el[3:5] == chars:
                    x_mod.append(self.y[el])
                    z_mod.append(self.z[el])
            else:
                print "get the direction right, either 1 or two"

        xs = np.asarray(x_mod)
        zs = np.asarray(z_mod)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        #ax.set_aspect('equal')
        p1 = ax.scatter(xs,zs, s=20)
        if direc == 1:
            ax.set_xlabel('x')
        else:
            ax.set_xlabel('y')

        ax.set_ylabel('z')
        if show == True:
            plt.show()
        plt.savefig( str(direc) + '_' + str(ind) + 'scatter_cells_slice.png')
        plt.clf()
        plt.close()
        print "Scatter cells slice complete"
        return 0
    # end T2InputGrid class

def plot_relperm_cap(rp_vals, cp_vals, fmt = 'png', rp = 'linear'):

    print "PLOTTING RELATIVE PERMEABILITY AND CAPILLARY PRESSURE CURVES"
    nvals = 100
    if rp == 'linear':
        sbres = rp_vals[0]
    else:
        sbres = rp_vals[1]
    sat = np.linspace(sbres, 1., nvals)
    pcap = np.zeros(nvals)
    krl = np.zeros(nvals)
    krg = np.zeros(nvals)

    for i in range(len(sat)):
        pcap[i] = -cap_vangenuchten(sat[i], cp_vals)

    if rp == 'linear':
        for i in range(len(sat)):
            krl[i], krg[i] = rel_perms_linear(sat[i], rp_vals)
    else:
        for i in range(len(sat)):
            krl[i], krg[i] = rel_perms_vangenuchten(sat[i], rp_vals)

    font = { 'size'   : 14}
    matplotlib.rc('font', **font)
    f = plt.figure(num=None , dpi = 480, \
        facecolor = 'w', edgecolor = 'k')
    ax = f.add_subplot(111)
    ax.set_xlabel('Sw []')
    ax.set_ylabel('Capillary Pressure [Pa]')
    p = plt.plot(sat, pcap)
    f.savefig('capillary_pressure' + '.' + fmt)
    plt.clf()
    plt.close()

    g = plt.figure(num=None , dpi = 480, \
        facecolor = 'w', edgecolor = 'k')
    #g.suptitle('Relative Permeability Curves')
    ax = g.add_subplot(111)
    ax.set_xlabel('Sw []')
    ax.set_ylabel('Relative Permeability []')
    p = plt.plot(sat, krl, label = 'liquid')
    p = plt.plot(sat, krg, label = 'gas')
    plt.legend()
    g.savefig('rel_perm_curves' + '.' + fmt)
    plt.clf()
    plt.close()
    return 0

def cap_linear(s, cp_vals):
    if cp_vals[2] < cp_vals[1]:
        print 'cp_vals[2] must be greater than cp_vals[1]'
        return 1
    if s <= cp_vals[1]:
        pcap = -cp_vals[0]
    elif s >= cp_vals[2]:
        pcap = 0
    else: 
        pcap = -cp_vals[0] *(cp_vals[2] - s) / (cp_vals[2] - cp_vals[1])
    return pcap

def cap_vangenuchten(sl, cp_vals):
    lamb = cp_vals[0]
    slr = cp_vals[1]
    p_0 = 1. / cp_vals[2]
    pmax = cp_vals[3]
    sls = cp_vals[4]

    ss = (sl - slr) / (sls - slr)

    pcap = -p_0 * pow(pow(ss,-1./lamb) - 1, 1. - lamb)
    if pcap < -pmax:
        pcap = -pmax
    elif pcap > 0.:
        pcap = 0.

    return pcap

def rel_perms_linear(sl, rp_vals):
    sg = 1. - sl
    l_slope = 1 / (rp_vals[2] - rp_vals[0])
    g_slope = 1 / (rp_vals[3] - rp_vals[1])

    if sl <= rp_vals[0]:
        krl = 0.
    elif sl >= rp_vals[2]:
        krl = 1.
    else:
        krl = (sl - rp_vals[0]) * l_slope

    if sg <= rp_vals[1]:
        krg = 0.
    elif sg >= rp_vals[3]:
        krg = 1.
    else:
        krg = (sg - rp_vals[1]) * g_slope

    return krl, krg

def rel_perms_vangenuchten(sl, rp_vals):
    lamb = rp_vals[0]
    s_lr = rp_vals[1]
    s_ls = rp_vals[2]
    s_gr = rp_vals[3]

    ss = (sl - s_lr) / (s_ls - s_lr)
    sh = (sl - s_lr) / ( 1. - s_lr - s_gr)

    if sl < s_ls:
        krl = np.sqrt(ss) * pow(1 - pow(1 - pow(ss, 1./lamb), lamb), 2.)
    else:
        krl = 1.

    if s_gr == 0.:
        krg = 1. - krl
    elif s_gr > 0.:
        krg = pow(1. - sh, 2.) * (1 - pow(sh, 2.))
    else: 
        print "s_gr < 0 in rel_perms_vangenuchten"
        return 1

    return krl, krg

def clean_oldmesh_file():
    """ takes the mesh file from a previous simulation and removes the excess
    text
        The old meshfile must be called OLDMESH
    """
    old = open('OLDMESH', 'r')
    new = open('MESH', 'w')
    line = old.readline()
    s = line.split()
    while s != []: 
        new.write(line)
        line = old.readline()
        s = line.split()
    new.write(line)
    line = old.readline()
    s = line.split()
    while s[0] != '+++':
        new.write(line)
        line = old.readline()
        s = line.split()
    new.write('\r\n')
    old.close()
    new.close()

if __name__ == '__main__':
    # tests the initialization of an object
    nx = 65
    ny = 119
    nz = 43
    grid = T2InputGrid(nx, ny, nz)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                char = grid.get_element_chars(i, j, k)
                if char == 'JH732':
                    print char, 'JH732'
                    print i, j, k
