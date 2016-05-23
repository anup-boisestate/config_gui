import gi
from gi.repository import Gtk

gi.require_version('Gtk', '3.0')
from utilities import *


class GIN3DConfigWriter:
    def __init__(self, builder):
        self.cnfg = {}
        self.win = builder

    def to_file(self, filename):
        print("Filename: ", filename)
        self.get_data()
        self.writeToFile(filename)

    def text(self, obj_id, evaluate=False):
        obj = self.win.get_object(obj_id)
        if obj is not None:
            if isinstance(obj, Gtk.ComboBoxText):
                return obj.get_active_id()
            elif isinstance(obj, Gtk.FileChooserButton):
                return obj.get_filename()

            else:
                val = obj.get_text()
                if val != '' and evaluate:
                    return str(eval(val))
                return val
        return ''

    def active(self, obj_id):
        return self.win.get_object(obj_id).get_active()

    def get_obj(self, obj_id):
        return self.win.get_object(obj_id)

    # Read Data from the GUI
    def get_data(self):
        self.cnfg['Main'] = self.get_main_parameters()
        self.cnfg['BC'] = self.get_bc_parameters()
        self.cnfg['IC'] = self.get_ic_parameters()
        self.cnfg['FS'] = self.get_fs_parameters()
        self.cnfg['SO'] = self.get_so_parameters()
        self.cnfg['TP'] = self.get_tp_parameters()
        self.cnfg['TM'] = self.get_tm_parameters()
        self.cnfg['GE'] = self.get_ge_parameters()
        print(self.cnfg)

    def get_main_parameters(self):
        cnfg_main = {'Grid': "%s %s %s" % (self.text("domain_lx"), self.text("domain_ly"), self.text("domain_lz")),
                     'Dimensions': "%s %s %s" % (self.text("mesh_nx"), self.text("mesh_ny"), self.text("mesh_nz")),
                     'ReferenceLength': self.text("refLength", True), 'Nu': self.text("laminar_viscosity", True),
                     'StopTime': self.text("sim_endTime"), 'Temperature': self.active("simParam_switch_temperature"),
                     'Turbulence': self.active("simParam_switch_turbulence"),
                     'SolidGeometry': self.active("simParam_switch_geometry")}

        self.forcing = "off"  # variable set to check when writing to file
        forcing_off = self.active("simParam_rb_forcing_off")

        if not forcing_off:
            forcing_cpg = self.active("simParam_rb_forcing_constPresGrad")
            if forcing_cpg:
                self.forcing = "constPresGrad"
                cnfg_main['Forcing'] = "%s %s %s" % (
                    self.text("simParam_txtbox_fx", True), self.text("simParam_txtbox_fy", True),
                    self.text("simParam_txtbox_fz", True))
            else:
                self.forcing = "constMassFlowRate"
                cnfg_main['NumberForcingRegions'] = self.text("no_of_forcing_regions")
                numforcings = int(cnfg_main['NumberForcingRegions'])
                heights = [0.0]
                mdot = []
                for i in range(numforcings):
                    h = first_child("simParam_forcing_cmfr_box", "height_", i, self.win)
                    if h is not None:
                        heights.append(eval(h.get_text()))
                        tmp = []
                        x = first_child("simParam_forcing_cmfr_box", "frfr_x_", i, self.win)
                        if x is not None:
                            tmp.append(str(eval(x.get_text())))

                        y = first_child("simParam_forcing_cmfr_box", "frfr_y_", i, self.win)
                        if y is not None:
                            tmp.append(str(eval(y.get_text())))

                        z = first_child("simParam_forcing_cmfr_box", "frfr_z_", i, self.win)
                        if z is not None:
                            tmp.append(str(eval(z.get_text())))

                        mdot.append(tmp)

                cnfg_main['ForcingRegionBounds'] = heights
                cnfg_main['ForcingRegionFlowRate'] = mdot

        return cnfg_main

    def get_bc_parameters(self):
        cnfg_bc = {'Face_West': self.text("west_face"), 'Face_East': self.text("east_face"),
                   'Face_South': self.text("south_face"), 'Face_North': self.text("north_face"),
                   'Face_Bottom': self.text("bottom_face"), 'Face_Top': self.text("top_face"),
                   'InletVelocities': ("%s %s %s") % (self.text("inletProfile_u"),
                                                      self.text("inletProfile_v"), self.text("inletProfile_w")),
                   'VariableProfile': self.active("inletProfile_switch_variable"),
                   'PerturbTurbInflow': self.active("inletProfile_switch_turbulent"),
                   'PerturbZoneOffset': self.text('inletProfile_offset'),
                   'PerturbZoneRichardson': self.text('inletProfile_richardson'),
                   'PerturbZoneDepth': self.text('inletProfile_depth'),
                   'PerturbBoxLength': self.text('inletProfile_length'),
                   'PerturbBoxWidth': self.text('inletProfile_width'),
                   'PerturbBoxHeight': self.text('inletProfile_height')}

        if cnfg_bc['VariableProfile'] is True:
            val = "Inlet"
            if self.active("rb_inletProfile_roughWall"):
                val = "RoughLogLawInlet"
            elif self.active("rb_inletProfile_smoothWall"):
                val = "SmoothLogLawInlet"
            elif self.active("rb_inletProfile_powerLaw"):
                val = "OneSeventhsInlet"
            elif self.active("rb_inletProfile_parabolic"):
                val = "Parabolic"
            elif self.active("rb_inletProfile_dns"):
                val = "DNS2000"

            if cnfg_bc['Face_South'] == "Inlet":
                cnfg_bc['Face_South'] = val

            if cnfg_bc['Face_West'] == "Inlet":
                cnfg_bc['Face_West'] = val

        return cnfg_bc

    def get_ic_parameters(self):
        cnfg_ic = {}

        if self.active('ic_rb_sameAsInletProfile'):
            # removed 'Parabolic' since it has not been implemented
            inlet_values = ['Inlet', 'RoughLogLawInlet', 'SmoothLogLawInlet', 'OneSeventhsInlet', 'DNS2000']

            cnfg_bc = self.cnfg['BC']
            faceW = cnfg_bc['Face_West']
            faceS = cnfg_bc['Face_South']

            if faceW in inlet_values or faceS in inlet_values:
                if faceW == "Inlet" or faceS == "Inlet":
                    cnfg_ic['InitialCondition'] = 'Uniform ' + self.cnfg['BC']['InletVelocities']
                elif faceW == "RoughLogLawInlet" or faceS == "RoughLogLawInlet":
                    cnfg_ic['InitialCondition'] = 'ComplexTerrain'
                elif faceW == "SmoothLogLawInlet" or faceS == "SmoothLogLawInlet":
                    cnfg_ic['InitialCondition'] = 'TurbulentChannelXMean'
                    # Don't have a power initial condition implemented yet, so we stick with smooth wall log law
                elif faceW == "OneSeventhsInlet" or faceS == "OneSeventhsInlet":
                    cnfg_ic['InitialCondition'] = 'TurbulentChannelXMean'
                    # Same with DNS2000
                elif faceW == "DNS2000" or faceS == "DNS2000":
                    cnfg_ic['InitialCondition'] = 'TurbulentChannelXMean'
            else:
                self.get_obj('ic_rb_sameAsInletProfile').set_active('False')
        elif self.active('ic_rb_zeroField'):
            cnfg_ic['InitialCondition'] = 'ZeroField'
        elif self.active('ic_rb_uniform'):
            cnfg_ic['InitialCondition'] = 'Uniform ' + self.text('ic_uniform_u') + " " + self.text(
                'ic_uniform_v') + " " + self.text('ic_uniform_w')

        elif self.active('ic_rb_specialCases'):
            special_case = self.text('ic_list_specialCases')
            # Smooth Wall Log Law w/Sinusoidal Perturbations
            if special_case == '1':
                cnfg_ic['InitialCondition'] = 'ComplexTerrainSinPerturb'
            # Rough Wall Log Law w/Sinusoidal Perturbations
            elif special_case == '2':
                cnfg_ic['InitialCondition'] = 'TurbulentChannelX'
            # ABC Flow
            elif special_case == '3':
                cnfg_ic['InitialCondition'] = 'ABCFlow'

        return cnfg_ic

    def get_fs_parameters(self):

        cnfg_fs = {'AdvectionScheme': '0.00',
                   'TemperatureAdvectionScheme': '0.00'}

        # Advection Scheme

        # Momentum
        if self.active('fsp_rb_momentum_cds'):
            cnfg_fs['AdvectionScheme'] = '0.00'
        elif self.active('fsp_rb_momentum_fou'):
            cnfg_fs['AdvectionScheme'] = '1.00'
        elif self.active('fsp_rb_momentum_quick'):
            cnfg_fs['AdvectionScheme'] = '-1.00'
        elif self.active('fsp_rb_momentum_kappa'):
            cnfg_fs['AdvectionScheme'] = ''
        elif self.active('fsp_rb_momentum_hybrid'):
            cnfg_fs['AdvectionScheme'] = self.text('fsp_text_momentum_hybrid')

        # Temperature
        if self.active('fsp_rb_temperature_cds'):
            cnfg_fs['TemperatureAdvectionScheme'] = '0.00'
        elif self.active('fsp_rb_temperature_fou'):
            cnfg_fs['TemperatureAdvectionScheme'] = '1.00'
        elif self.active('fsp_rb_temperature_quick'):
            cnfg_fs['TemperatureAdvectionScheme'] = '-1.00'
        elif self.active('fsp_rb_temperature_kappa'):
            cnfg_fs['TemperatureAdvectionScheme'] = ''
        elif self.active('fsp_rb_temperature_hybrid'):
            cnfg_fs['TemperatureAdvectionScheme'] = self.text('fsp_text_temperature_hybrid')

        # Time Marching
        cnfg_fs['TimeMethod'] = self.text('fsp_timeMarching')

        # Time Step Size
        cnfg_fs['TimeStepSize'] = self.text('fsp_timeStepSize')

        if cnfg_fs['TimeStepSize'] == 'Variable':
            cnfg_fs['ConstTimeStep'] = '-0.1'
            cnfg_fs['CFL'] = self.text('fsp_cfl')
            cnfg_fs['DTStability'] = self.text('fsp_viscousLimit')
        else:
            cnfg_fs['ConstTimeStep'] = self.text('fsp_dt')
            cnfg_fs['CFL'] = '0.4'
            cnfg_fs['DTStability'] = '0.4'

        # Poisson Solver
        cnfg_fs['SolverMethod'] = self.text('fsp_poissonSolver')

        if cnfg_fs['SolverMethod'] == 'Multigrid':

            cnfg_fs['MG_Loops'] = self.text('fsp_text_poissonSolver_outerLoops')
            cnfg_fs['MG_Levels'] = self.text('fsp_text_poissonSolver_mgLevels')
            cnfg_fs['MG_Smoother_Iterations'] = [self.text('fsp_text_poissonSolver_innerLoops_down'),
                                                 self.text('fsp_text_poissonSolver_innerLoops_up')]
            cnfg_fs['MG_Jacobi_Weight'] = self.text('fsp_text_poissonSolver_mg_jacobiWeight')

        else:
            cnfg_fs['Iterative_Loops'] = self.text('fsp_text_poissonSolver_iterations')
            cnfg_fs['Iterative_Weight'] = self.text('fsp_text_poissonSolver_pj_jacobiWeight')
        return cnfg_fs

    def get_so_parameters(self):
        cnfg_so = {'OutputPath': self.text('so_folder_outputPath'),
                   'OutputPrefix': self.text('so_text_outputPrefix')
                   }

        # Restart Solution
        if self.active('so_switch_restartSolution'):
            cnfg_so['RestartDirectory'] = self.text('so_folder_restartDirectory')

        outputFormat = []

        # Formatted Output
        if self.active('so_cb_fileFormat_hdf5'):
            outputFormat.append('HDF5_Multiple')

        # need to confirm this
        if self.active('so_cb_fileFormat_vti'):
            outputFormat.append('VTK')

        if self.active('so_cb_fileFormat_hdf5') or self.active('so_cb_fileFormat_vti'):

            if self.active('so_cb_formatOutput_timeStep'):
                cnfg_so['OutputFrequencyVTK1'] = self.text('so_text_formatOutput_timeStep_interval1')
                if self.active('so_cb_formatOutput_timeStep_switch'):
                    cnfg_so['OutputFrequencyVTK2'] = self.text('so_text_formatOutput_timeStep_interval2')
                    cnfg_so['OutputFrequencyVTKSwitch'] = self.text('so_text_formatOutput_timeStep_switch')

            if self.active('so_cb_formatOutput_physicalTime'):
                cnfg_so['OutputTimeIntervalVTK1'] = self.text('so_text_formatOutput_physicalTime_interval1')
                if self.active('so_cb_formatOutput_physicalTime_switch'):
                    cnfg_so['OutputTimeIntervalVTK2'] = self.text('so_text_formatOutput_physicalTime_interval2')
                    cnfg_so['OutputTimeIntervalVTKSwitch'] = self.text('so_text_formatOutput_physicalTime_switch')

            outputVariablesVTK = []

            if self.active('so_cb_outputVTK_u'):
                outputVariablesVTK.append('u')
            if self.active('so_cb_outputVTK_p'):
                outputVariablesVTK.append('p')
            if self.active('so_cb_outputVTK_t'):
                outputVariablesVTK.append('t')
            if self.active('so_cb_outputVTK_n'):
                outputVariablesVTK.append('n')
            if self.active('so_cb_outputVTK_z'):
                outputVariablesVTK.append('z')
            if self.active('so_cb_outputVTK_q'):
                outputVariablesVTK.append('q')
            if self.active('so_cb_outputVTK_g'):
                outputVariablesVTK.append('g')
            if self.active('so_cb_outputVTK_c'):
                outputVariablesVTK.append('c')
            if self.active('so_cb_outputVTK_i'):
                outputVariablesVTK.append('i')

            cnfg_so['OutputVariablesVTK'] = outputVariablesVTK

        # Raw Data
        if self.active('so_cb_rawData_binary'):
            outputFormat.append('Binary')
        elif self.active('so_cb_rawData_parallelBinary'):
            outputFormat.append('PBinary')

        if self.active('so_cb_rawData_ASCII'):
            outputFormat.append('VTK')
        elif self.active('so_cb_rawData_parallelASCII'):
            outputFormat.append('PVTK')

        if len(set(outputFormat).intersection(set(['Binary', 'PBinary', 'VTK', 'PVTK']))) > 0:
            if self.active('so_cb_rawData_timeStep'):
                cnfg_so['OutputFrequencyText1'] = self.text('so_text_rawData_timeStep_interval1')
                if self.active('so_cb_rawData_timeStep_switch'):
                    cnfg_so['OutputFrequencyText2'] = self.text('so_text_rawData_timeStep_interval2')
                    cnfg_so['OutputFrequencyTextSwitch'] = self.text('so_text_rawData_timeStep_switch')

            if self.active('so_cb_rawData_physicalTime'):
                cnfg_so['OutputTimeIntervalText1'] = self.text('so_text_rawData_physicalTime_interval1')
                if self.active('so_cb_formatOutput_physicalTime_switch'):
                    cnfg_so['OutputTimeIntervalText2'] = self.text('so_text_rawData_physicalTime_interval2')
                    cnfg_so['OutputTimeIntervalTextSwitch'] = self.text('so_text_rawData_physicalTime_switch')

            outputVariablesTXT = []

            if self.active('so_cb_outputTXT_u'):
                outputVariablesTXT.append('u')
            if self.active('so_cb_outputTXT_p'):
                outputVariablesTXT.append('p')
            if self.active('so_cb_outputTXT_t'):
                outputVariablesTXT.append('t')
            if self.active('so_cb_outputTXT_n'):
                outputVariablesTXT.append('n')
            if self.active('so_cb_outputTXT_z'):
                outputVariablesTXT.append('z')
            if self.active('so_cb_outputTXT_q'):
                outputVariablesTXT.append('q')
            if self.active('so_cb_outputTXT_g'):
                outputVariablesTXT.append('g')
            if self.active('so_cb_outputTXT_c'):
                outputVariablesTXT.append('c')
            if self.active('so_cb_outputTXT_i'):
                outputVariablesTXT.append('i')

            cnfg_so['OutputVariablesTXT'] = outputVariablesTXT

        # Screen Dump
        if self.active('so_cb_screenDump'):
            outputFormat.append('Matrix')

        # Output Plane
        if self.active('so_rb_screenDump_xy'):
            cnfg_so['OutputPlane'] = 'XY'
            cnfg_so['OutputSlice'] = self.text('so_text_screenDump_xy')
        elif self.active('so_rb_screenDump_xz'):
            cnfg_so['OutputPlane'] = 'XZ'
            cnfg_so['OutputSlice'] = self.text('so_text_screenDump_xz')
        elif self.active('so_rb_screenDump_yz'):
            cnfg_so['OutputPlane'] = 'YZ'
            cnfg_so['OutputSlice'] = self.text('so_text_screenDump_yz')

        if self.active('so_rb_screenDump_u'):
            cnfg_so['MatrixData'] = 'u'
        elif self.active('so_rb_screenDump_v'):
            cnfg_so['MatrixData'] = 'v'
        elif self.active('so_rb_screenDump_w'):
            cnfg_so['MatrixData'] = 'w'
        elif self.active('so_rb_screenDump_p'):
            cnfg_so['MatrixData'] = 'p'
        elif self.active('so_rb_screenDump_temperature'):
            cnfg_so['MatrixData'] = 'phi'

        cnfg_so['OutputFormat'] = outputFormat

        # Sampling
        # -- need to do this section

        # Time Series OPtions
        # Collect Statistics
        cnfg_so['TimeAvg'] = str(self.active('so_switch_collectStatistics'))

        # Sampling Time Begin
        cnfg_so['StartingTime'] = self.text('so_text_samplingTimeBegin')

        # Spatial Averaging
        spatialAvg = self.active('so_switch_spatialAveraging')
        # -- need to review this
        turbStats = []
        if self.active('so_rb_firstOrderStatistics'):
            turbStats.append('m')
        else:
            if self.active('so_cb_stats_s'):
                turbStats.append('s')
            if self.active('so_cb_stats_d'):
                turbStats.append('d')
            if self.active('so_cb_stats_i'):
                turbStats.append('i')

        if spatialAvg:
            turbStats = [e.upper() for e in turbStats]

        cnfg_so['TurbStats'] = turbStats

        # Mean quantity files for fluctuation calculations
        if self.active('so_cb_mean_velocity'):
            cnfg_so['MeanProfile_U'] = self.text('so_file_meanprofile_u')
            cnfg_so['MeanProfile_V'] = self.text('so_file_meanprofile_v')
            cnfg_so['MeanProfile_W'] = self.text('so_file_meanprofile_w')

        if self.active('so_cb_mean_temperature'):
            cnfg_so['MeanProfile_Phi'] = self.text('so_file_meanprofile_phi')

        if self.active('so_cb_mean_pressure'):
            cnfg_so['MeanProfile_P'] = self.text('so_file_meanprofile_p')

        return cnfg_so

    def get_tp_parameters(self):
        cnfg_tp = {}
        if self.cnfg['Main']['Temperature']:
            cnfg_tp = {'Temp_West' : [self.text('temp_west_face'), self.text('temp_text_westFace')],
                       'Temp_East': [self.text('temp_east_face'), self.text('temp_text_eastFace')],
                       'Temp_North': [self.text('temp_north_face'), self.text('temp_text_northFace')],
                       'Temp_South': [self.text('temp_south_face'), self.text('temp_text_southFace')],
                       'Temp_Top': [self.text('temp_top_face'), self.text('temp_text_topFace')],
                       'Temp_Bottom': [self.text('temp_bottom_face'), self.text('temp_text_bottomFace')],
                       'Beta' : self.text('temp_text_thermalExpansionCoefficient'),
                       'Gamma' : self.text('temp_text_thermalDiffusivity'),
                       'Source_PHI' : self.text('temp_text_soruceMagnitude')}

        return cnfg_tp

    def get_tm_parameter(self):
        cnfg_tm = {}
        if self.cnfg['Main']['Turbulence']:
            cnfg_tm = {}
        return cnfg_tm

    def get_ge_parameter(self):
        cnfg_ge = {}
        if self.cnfg['Main']['SolidGeometry']:
            cnfg_ge = {}
        return cnfg_ge

    # Write Data to File
    def writeToFile(self, filename):

        try:
            fp = open(filename.replace('.cfg', '') + ".cfg", "w")

            # Main Parameters
            cnfg_main = self.cnfg["Main"]

            # Item 1
            fp.write("# Physical domain size, requires three floating point arguments after\n")
            fp.write("# Dimensions, separated by spaces.\n")
            fp.write("Dimensions " + cnfg_main["Dimensions"])
            fp.write("\n\n")

            # Item 2
            fp.write("# Simple mesh size, requires three integer arguments after Grid, separated by\n")
            fp.write("# spaces.\n")
            fp.write("Grid " + cnfg_main["Grid"])
            fp.write("\n\n")

            # Item 3
            fp.write("# Reference Length\n")
            fp.write("# This is used as the characteristic length 'L' in calculations of the\n")
            fp.write("# Reynolds and Rayleigh numbers.  Leave unset if you do not know what it is.\n")
            fp.write("#\n")
            fp.write("# Enter either 'LX', 'LY', 'LZ', or a number.\n")
            fp.write("ReferenceLength " + cnfg_main["ReferenceLength"])
            fp.write("\n\n")

            # Item 4
            fp.write("# Kinematic viscosity\n")
            fp.write("Nu " + cnfg_main["Nu"])
            fp.write("\n\n")

            # Item 5
            fp.write("# Specify a specific physical time (in seconds) for stopping the simulation.\n")
            fp.write("StopTime  " + cnfg_main['StopTime'])
            fp.write("\n\n")

            # Item 6
            # do nothing, will use TurbulenceModel tag later to choose model

            # Item 7
            fp.write("# Whether or not to solve for and output temperature\n")
            fp.write("Temperature   " + str(cnfg_main['Temperature']))
            fp.write("\n\n")

            # Item 8
            # Don't do anything for this one yet
            # for testing
            fp.write("Turbulence   " + str(cnfg_main['Turbulence']))

            # Item 9
            # Do nothing
            # for testing
            fp.write("SolidGeometry   " + str(cnfg_main['SolidGeometry']))

            # Item 10
            if self.forcing == "constPresGrad":
                fp.write("# Constant forcing, requires three floating point arguments after Forcing,\n")
                fp.write("# separated by spaces.\n")
                fp.write("Forcing   " + cnfg_main["Forcing"])
                fp.write("\n\n")

            # Item 11
            if self.forcing == "constMassFlowRate":
                fp.write("# To drive a periodic flow with a constant mass flow rate, we can define \n")
                fp.write("# different regions that will be given independent forcing values. The regions\n")
                fp.write("# are defined by the ForcingRegionBounds tag which requires two values, the\n")
                fp.write("# first for the lower bound and the second for the upper bound. The bounds are\n")
                fp.write("# determined by the distance field so ensure that is loaded. The forcing will \n")
                fp.write("# be adjusted by the difference between a prescribed mass flow rate and the \n")
                fp.write("# simulated mass flow rate. The prescribed mass flow rate is initially set by\n")
                fp.write("# the three values given to the ForcingRegionFlowRate tag. The three values\n")
                fp.write("# correspond to the x-, y- and z-directions, respectively.  The number of \n")
                fp.write("# ForcingRegionBounds and ForcingRegionFlowRate tag sets need to be equal to\n")
                fp.write("# the number of regions given to the NumberForcingRegions tag.  The first\n")
                fp.write("# ForcingRegionFlowRate tag applies to the first ForcingRegionBounds tag.\n")
                fp.write("ConstantMassFlowRate  True\n")

                numforcings = int(cnfg_main["NumberForcingRegions"])
                fp.write("NumberForcingRegions   " + str(numforcings) + "\n")

                heights = cnfg_main["ForcingRegionBounds"]
                mdot = cnfg_main["ForcingRegionFlowRate"]
                for i in range(numforcings):
                    fp.write("ForcingRegionBounds    %s %s\n" % (str(heights[i]), str(heights[i + 1])))
                    fp.write("ForcingRegionFlowRate  %s %s %s \n" % (str(mdot[i][0]), str(mdot[i][1]), str(mdot[i][2])))
                fp.write("\n")

            # Boundary Conditions
            cnfg_bc = self.cnfg["BC"]

            # Items
            fp.write("#-------------------------------------------------------------------------------\n")
            fp.write("#                          Domain Boundary Conditions\n")
            fp.write("#-------------------------------------------------------------------------------\n")
            fp.write("\n")
            fp.write("#    Face type           Boundary setting\n")
            fp.write("#    ------------------  -------------------------------------\n")
            fp.write("#    NoSlip              velocity 0 at boundary\n")
            fp.write("#    FreeSlip            velocity unchanged at boundary\n")
            fp.write("#    Inlet               BC = inlet\n")
            fp.write("#    Outlet              BC = interior velocity\n")
            fp.write("#    ConvectiveOutlet    See Ferziger (2001)\n")
            fp.write("#    Driven              BC = inlet in normal direction\n")
            fp.write("#    Periodic            BC = opposite side\n")
            fp.write("#    SmoothLogLawInlet   an inlet with a constant smooth-wall logarithmic profile\n")
            fp.write("#                          in the z-direction.\n")
            fp.write("#    RoughLogLawInlet    an inlet with a constant rough-wall logarithmic profile\n")
            fp.write("#                          in the z-direction.\n")
            fp.write("#    OneSeventhsInlet    an inlet with a constant 1/7 power law profile\n")
            fp.write("#                          in the z-direction.\n")
            fp.write("#    DataBase            inflow database BC, usage: DataBase path_to_database\n")
            fp.write("#    WallModel           Set boundary condition to Schumann wall model \n")
            fp.write("\n")
            fp.write("# Only one each of an Inlet and Driven are allowed.\n")
            fp.write("\n")
            fp.write("# Bottom is Z=0, Top   is in the Z+ direction (w component of velocity).\n")
            fp.write("# South  is Y=0, North is in the Y+ direction (v component of velocity).\n")
            fp.write("# West   is X=0, East  is in the X+ direction (u component of velocity).\n")
            fp.write("\n")
            fp.write("Face_West     " + cnfg_bc['Face_West'] + "\n")
            fp.write("Face_East     " + cnfg_bc['Face_East'] + "\n")
            fp.write("Face_South    " + cnfg_bc['Face_South'] + "\n")
            fp.write("Face_North    " + cnfg_bc['Face_North'] + "\n")
            fp.write("Face_Bottom   " + cnfg_bc['Face_Bottom'] + "\n")
            fp.write("Face_Top      " + cnfg_bc['Face_Top'] + "\n")
            fp.write("\n")

            # Items
            fp.write("# InletVelocities requires three floating point values separated by spaces.\n")
            fp.write("# These values apply to Inlet and LogarithmicInlet boundary conditions. With \n")
            fp.write("# the LogarithmicInlet, it is only used to determine the angle of the flow so\n")
            fp.write("# magnitude does not matter.  The magnitude is found from the FrictionVelocity.\n")
            fp.write("# Note: These are signed vectors, negative is west/south\n")
            fp.write("\n")
            fp.write("InletVelocities  " + cnfg_bc["InletVelocities"] + "\n")
            fp.write("\n")

            # Items
            fp.write("# The perturbation cell turbulent inflow conditions (first proposed by \n")
            fp.write("# Munoz-Esparza et al. 2014) applies random perturbations to groups of \n")
            fp.write("# temperature grid cells referred to as perturbation cells. These\n")
            fp.write("# perturbation cells need their length, width and height defined with\n")
            fp.write("# PerturbBoxLength, PerturbBoxWidth and PerturbBoxHeight, respectively.\n")
            fp.write("# The number of perturbation cells from the inflow is defined by the\n")
            fp.write("# PerturbZoneDepth tag as an integer value.\n")
            fp.write("# The perturbation amplitude is defined by PerturbZoneRichardson.\n")
            fp.write("# The perturbation cell method is activated by setting a boundary to an\n")
            fp.write("# inlet and setting PerturbTurbInflow to True.\n")
            fp.write("# Note: PerturbBoxLength, PerturbBoxWidth and PerturbBoxHeight can\n")
            fp.write("# be set to the strings dx, dy and dz, respectively as a shortcut to\n")
            fp.write("# set the perturbation dimensions equal to the grid spacing.\n")
            fp.write("\n")
            fp.write("PerturbTurbInflow      " + str(cnfg_bc["PerturbTurbInflow"]) + "\n")
            fp.write("PerturbZoneOffset      " + cnfg_bc["PerturbZoneOffset"] + "\n")
            fp.write("PerturbZoneRichardson  " + cnfg_bc["PerturbZoneRichardson"] + "\n")
            fp.write("PerturbZoneDepth       " + cnfg_bc["PerturbZoneDepth"] + "\n")
            fp.write("PerturbBoxLength  " + cnfg_bc["PerturbBoxLength"] + "\n")
            fp.write("PerturbBoxWidth   " + cnfg_bc["PerturbBoxWidth"] + "\n")
            fp.write("PerturbBoxHeight  " + cnfg_bc["PerturbBoxHeight"] + "\n")
            fp.write("\n")

            # Initial Conditions
            cnfg_ic = self.cnfg["IC"]

            fp.write("#-------------------------------------------------------------------------------\n")
            fp.write("#                             Initial Conditions\n")
            fp.write("#-------------------------------------------------------------------------------\n")
            fp.write("\n")
            fp.write("# Set the initial condition. Choose from the following:\n")
            fp.write("#    File                      Load from file. Fill in Initial_* tags below.\n")
            fp.write("#    TaylorGreenVortex         u =  sin(2*PI*x)*cos(2*PI*y)\n")
            fp.write("#                              v = -cos(2*PI*x)*sin(2*PI*y)\n")
            fp.write("#                              w = 0.0\n")
            fp.write("#    TurbulentChannelX         Note: Give FrictionVelocity and ReferenceVelocity tags\n")
            fp.write(
                "#                              u = u_tau ( 1/kappa * log(z+) + 5.2) + 2.0 * C/wx * cos(wx*x)*sin(wy*y)sin(wz*z)\n")
            fp.write("#                              v = -C/wy * sin(wx*x)*cos(wy*y)sin(wz*z)\n")
            fp.write("#                              w = -C/wz * sin(wx*x)*sin(wy*y)cos(wz*z)\n")
            fp.write("#                              wx = round(0.5*LX) * Pi / LX\n")
            fp.write("#                              wy = round(0.5*LY) * Pi / LY\n")
            fp.write("#                              wz = round(3.0*LZ) * Pi / LZ\n")
            fp.write("#                              C = 0.15 * Vref\n")
            fp.write("#    TurbulentChannelXMean     u = u_tau ( 1/kappa * log(z+) + 5.2)\n")
            fp.write("#                              v = 0.0\n")
            fp.write("#                              w = 0.0\n")
            fp.write("#    TurbulentChannelXPoise    u = 1.5 * Uinlet * ( 1 - y^2 / h^2 ) + perturbation\n")
            fp.write("#                              v = perturbation\n")
            fp.write("#                              w = perturbation\n")
            fp.write(
                "#    TurbulentChannelZ         u =  u_tau( 1/kappa * log(x+) + 5.5) + sin(Pi*x)*cos(2*z)*sin(2*y)\n")
            fp.write("#                              v =  -( 1 + cos(Pi*x) ) * sin(z) * sin(4.1*y)\n")
            fp.write("#                              w = -0.5*Pi*sin(z) * sin(Pi*x) * cos(1.25*Pi*y)\n")
            fp.write("#    ComplexTerrain            Note: theta determined from InletVelocity u and v components\n")
            fp.write("#                              u = u_tau/kappa * log(z/z_not) * cos(theta)\n")
            fp.write("#                              v = u_tau/kappa * log(z/z_not) * sin(theta)\n")
            fp.write("#                              w = 0.0\n")
            fp.write("#    ComplexTerrainSinPerturb  Note: theta determined from InletVelocity u and v components\n")
            fp.write(
                "#                              u = u_tau/kappa * log(z/z_not) * cos(theta) + 2.0*cos(x)*sin(y)*sin(z)\n")
            fp.write(
                "#                              v = u_tau/kappa * log(z/z_not) * sin(theta) - sin(x)*cos(y)*sin(z)\n")
            fp.write("#                              w = -sin(x)*sin(y)*( 1.0 + cos(z) )\n")
            fp.write("#    ABCFlow                   u = cos(y) + sin(z)\n")
            fp.write("#                              v = sin(x) + cos(z)\n")
            fp.write("#                              w = cos(x) + sin(y)\n")
            fp.write("# An option not listed or no option specified at all results in an\n")
            fp.write("# initial field of zero velocity.\n")
            if self.active('ic_rb_sameAsInletProfile'):
                fp.write("#<Same as Inlet Profile>\n")
            fp.write("InitialCondition  " + cnfg_ic['InitialCondition'] + "\n")

            # Flow Solver Parameters
            cnfg_fs = self.cnfg['FS']

            fp.write("#-------------------------------------------------------------------------------\n")
            fp.write("#                             Flow Solver Parameters\n")
            fp.write("#-------------------------------------------------------------------------------\n")
            fp.write("\n")
            fp.write("#Choose the finite difference method for spatial derivatives\n")
            fp.write("#0.00 means CDS, 1.00 means FOU, values between are allowed\n")
            fp.write("AdvectionScheme " + cnfg_fs['AdvectionScheme'] + '\n')
            fp.write("\n\n")

            fp.write("TemperatureAdvectionScheme " + cnfg_fs['TemperatureAdvectionScheme'] + '\n')
            fp.write("\n\n")

            fp.write("# The time derivative method\n")
            fp.write("#   Euler        (first order forward Euler)\n")
            fp.write("#   AdamsBash    (second order Adams-Bashforth, aka 'AB2')\n")
            fp.write("TimeMethod  " + cnfg_fs['TimeMethod'] + '\n')
            fp.write("\n\n")

            fp.write("# Specify a computational time step. Set a to negative number if dynamic time\n")
            fp.write("# step adjustment is desired.\n")
            fp.write("ConstTimeStep   " + cnfg_fs["ConstTimeStep"] + '\n')
            fp.write("\n\n")

            fp.write("#CFL*dz/velmax is the convective dt limit\n")
            fp.write("CFL      " + cnfg_fs["CFL"] + '\n')
            fp.write("\n\n")

            fp.write("#DTStability*dz*dz/Nu is the viscous dt limit\n")
            fp.write("#DTStability   " + cnfg_fs['DTStability'] + '\n')
            fp.write("\n\n")

            fp.write("# Choose either Iterative or Multigrid\n")
            fp.write("SolverMethod	          " + cnfg_fs['SolverMethod'] + '\n')
            fp.write("\n")
            fp.write("# Reasonable values:\n")
            if cnfg_fs['SolverMethod'] == 'Multigrid':
                fp.write("#   Cycle         V, W, or W1, recommend V\n")
                fp.write("#   Smoother      Jacobi or GS\n")
                fp.write("#   Loops         1-10\n")
                fp.write("#   Levels        2-N, recommend about 2-3 less than the maximum\n")
                fp.write("#   SmoothIters   1-10 1-10. Good choices include '2 1' and '4 2'\n")
                fp.write("#   JacobiWeight  0.55 - 0.95, recommend 0.86 (from Trottenberg says 0.857)\n")
                fp.write("#   SORWeight     1.0 - 1.8\n")
                fp.write("MG_Cycle                V\n")
                fp.write("MG_Smoother             Jacobi\n")
                fp.write("MG_Loops                " + cnfg_fs['MG_Loops'] + "\n")
                fp.write("MG_Levels               " + cnfg_fs['MG_Levels'] + "\n")  # default should be 50
                fp.write("MG_Smoother_Iterations  " + ' '.join(cnfg_fs['MG_Smoother_Iterations']) + "\n")
                fp.write("MG_Jacobi_Weight        " + cnfg_fs['MG_Jacobi_Weight'] + "\n")
                fp.write("MG_SOR_Weight           1.00\n")
            else:
                fp.write("#            Solver   Jacobi     GS\n")
                fp.write("#            Loops    10-40      5-30\n")
                fp.write("#            Weight   1.0        1.0-1.8\n")
                fp.write("Iterative_Solver        Jacobi \n")
                fp.write("Iterative_Loops         " + cnfg_fs['Iterative_Loops'] + "\n")
                fp.write("Iterative_Weight        " + cnfg_fs['Iterative_Weight'] + "\n")
                fp.write("\n")

            # Solution Output Parameters
            cnfg_so = self.cnfg['SO']
            fp.write("\n")
            fp.write("#-------------------------------------------------------------------------------\n")
            fp.write("#                           Solution Output Parameters\n")
            fp.write("#-------------------------------------------------------------------------------\n")

            op = cnfg_so['OutputPath']
            if op is not None:
                fp.write("# An optional path for the results\n")
                fp.write("OutputPath     " + op + "\n")
            fp.write("# Optional prefix (default is 'gin3d_soln')\n")
            fp.write("OutputPrefix     " + cnfg_so['OutputPrefix'] + "\n")
            fp.write("\n")

            fp.write("# Which 2D plane to output (XY, XZ, YZ) or Mesh for the entire 3D domain.\n")
            fp.write("# Can select multiple planes.\n")
            fp.write("OutputPlane Mesh\n")
            fp.write("OutputPlane " + cnfg_so['OutputPlane'] + "\n")
            fp.write("\n")
            fp.write("# If a plane is selected, choose the location of the slice along the\n")
            fp.write("# direction perpendicular to the plane as a decimal (0.0-1.0)\n")
            fp.write("OutputSlice " + cnfg_so['OutputSlice'] + "\n")
            fp.write("\n")
            fp.write("# Data output for the Matrix output: u, v, w, p, phi\n")
            fp.write("MatrixData " + cnfg_so['MatrixData'] + "\n")
            fp.write("\n")

            try:
                rd = cnfg_so['RestartDirectory']
            except KeyError:
                rd = None
            if rd is not None:
                fp.write("# To restart a simulation, simply give the directory to the ouptut from the previous\n")
                fp.write("# simulation. The most recent timestep present in the directory will be chosen.\n")
                fp.write("# Full path starting from / is recommended. Ensure the last / is also included. This\n")
                fp.write("# does not automatically place in the / if missing.\n")
                fp.write("RestartDirectory   " + rd + "\n")

            fp.write("\n")
            fp.write("# Output formats.  Output formats can be combined by giving one line per type.\n")
            fp.write("# Supported format types are:\n")
            fp.write("#    Matrix          display values on the screen\n")
            fp.write("#    DataFile        writes data to an ASCII file with .dat suffix\n")
            fp.write("#    Binary          writes a single binary file using MPI-IO with a .bin\n")
            fp.write("#                       suffix\n")
            fp.write("#    PBinary         writes parallel binary files using MPI-IO with a .bin\n")
            fp.write("#                       suffix\n")
            fp.write("#    VTK             Visualization Toolkit Image file (.vti)\n")
            fp.write("#    PVTK            Parallel VTK Image file (.pvti and .vti). Number of vti\n")
            fp.write("#                       files generated is the same as the number of GPUs.\n")
            fp.write("#    HDF5_Multiple   Write out heavy data in HDF5 format with XDMF files for\n")
            fp.write("#                       visualization software. Produces multiple HDF5 files,\n")
            fp.write("#                       one for each time step.\n")
            fp.write("#\n")
            fp.write("# Multiple output formats can be requested (e.g. both Matrix and VTK). The\n")
            fp.write("# PVTK option cannot be combined with its sequential counterpart.  Doing so\n")
            fp.write("# defaults to the VTK option. Also, VTK cannot be chosen when using more than\n")
            fp.write("# eight GPUs and will default to PVTK.\n")
            [fp.write("OutputFormat  " + e + "\n") for e in cnfg_so['OutputFormat']]

            fp.write("\n")
            fp.write("# Two methods of controlling output frequency: iterations and physical time.\n")
            fp.write("# Setting an interval negative will disable that method of control.\n")
            fp.write("\n")
            fp.write("# Output by iteration:\n")
            fp.write("# Write files every N timesteps. OutputFrequencyVTK is for VTK and PVTK file\n")
            fp.write("# types. OutputFrequencyText is for DataFile, Binary, and PBinary file types.\n")
            fp.write("# The output frequency can be changed part way through a simulation by\n")
            fp.write("# setting OutputFrequencyVTKSwitch or OutputFrequencyTextSwitch. WTK output\n")
            fp.write("# frequency will be every OutputFrequencyVTK1 timesteps until\n")
            fp.write("# OutputFrequencyVTKSwitch timesteps is reached, then the frequency will be\n")
            fp.write("# every OutputFrequencyVTK2 timesteps.\n")
            fp.write("\n")
            try:
                fp.write("OutputFrequencyVTK1          " + cnfg_so['OutputFrequencyVTK1'] + "\n")
            except KeyError:
                pass
            try:
                fp.write("OutputFrequencyVTK2          " + cnfg_so['OutputFrequencyVTK2'] + "\n")
                fp.write("OutputFrequencyVTKSwitch     " + cnfg_so['OutputFrequencyVTKSwitch'] + "\n")
            except KeyError:
                pass
            try:
                fp.write("OutputFrequencyText1          " + cnfg_so['OutputFrequencyText1'] + "\n")
            except KeyError:
                pass
            try:
                fp.write("OutputFrequencyText2          " + cnfg_so['OutputFrequencyText2'] + "\n")
                fp.write("OutputFrequencyTextSwitch     " + cnfg_so['OutputFrequencyTextSwitch'] + "\n")
            except KeyError:
                pass
            fp.write("\n")
            fp.write("# Output by physical time:\n")
            fp.write("# Writes files every time interval (in seconds) rather than every N\n")
            fp.write("# See avobe for description of interval switching and corresponding file\n")
            fp.write("# types.\n")
            try:
                fp.write("OutputTimeIntervalVTK1          " + cnfg_so['OutputTimeIntervalVTK1'] + "\n")
            except KeyError:
                pass
            try:
                fp.write("OutputTimeIntervalVTK2          " + cnfg_so['OutputTimeIntervalVTK2'] + "\n")
                fp.write("OutputTimeIntervalVTKSwitch     " + cnfg_so['OutputTimeIntervalVTKSwitch'] + "\n")
            except KeyError:
                pass
            try:
                fp.write("OutputTimeIntervalText1          " + cnfg_so['OutputTimeIntervalText1'] + "\n")
            except KeyError:
                pass
            try:
                fp.write("OutputTimeIntervalText2          " + cnfg_so['OutputTimeIntervalText2'] + "\n")
                fp.write("OutputTimeIntervalTextSwitch     " + cnfg_so['OutputTimeIntervalTextSwitch'] + "\n")
            except KeyError:
                pass
            fp.write("\n")

            fp.write("# Output options for OutputVariablesVTK and OutputVariablesTXT tags:\n")
            fp.write("#   *   All variables\n")
            fp.write("#   @   Only variables needed for restarting simulation\n")
            fp.write("#   u   Velocity\n")
            fp.write("#   p   Pressure\n")
            fp.write("#   t   Temperature\n")
            fp.write("#   q   Q-criterion\n")
            fp.write("#   d   Density (PGM Solid Representation)\n")
            fp.write("#   g   IB Flags(Immersed Boundary Method)\n")
            fp.write("#   z   Distance field\n")
            fp.write("#   n   Eddy Viscosity\n")
            fp.write("#   c   Smagorinsky model coefficient (Dynamic C_s)\n")
            fp.write("#   i   Inflow Perturbation\n")
            fp.write("# Choose any combination of these, each character separated by spaces.  The\n")
            fp.write("# time-averaged values are output in the final files, if applicable.\n")
            fp.write("# For example, OutputVariablesVTK u p t will output velocity, pressure,\n")
            fp.write("# instantaneous temperature.\n")
            fp.write("# If no or incorrect options are given, * is the default.\n")
            fp.write("\n")
            fp.write("# Choose output variables for VTK\n")
            fp.write("OutputVariablesVTK  ")
            try:
                fp.write(" ".join(cnfg_so['OutputVariablesVTK']))
            except:
                pass
            fp.write("\n")
            fp.write("# Choose output variables for DataFile and/or Binary\n")
            fp.write("OutputVariablesTXT  ")
            try:
                fp.write(" ".join(cnfg_so['OutputVariablesTXT']))
            except:
                pass
            fp.write("\n")
            fp.write("\n")

            fp.write("# -------------------------------------------------------------------------------\n")
            fp.write("#                               Time Series Options\n")
            fp.write("# -------------------------------------------------------------------------------\n")
            fp.write("# Set averaging to true or false.  Set the starting time for averaging to start.\n")
            fp.write("# Starting time does not matter if PrecursorSim is set to True.\n")
            fp.write("TimeAvg      " + cnfg_so['TimeAvg'] + "\n")
            fp.write("StartingTime " + cnfg_so['StartingTime'] + "\n")
            fp.write("\n")
            fp.write("# Set statistical quantity to calculate with TurbStats tags:\n")
            fp.write("#   m   time averaging of primitive quantities\n")
            fp.write("#   s   symmetric components of Reynolds stress tensor (shear)\n")
            fp.write("#   d   diagonal components of Reynolds stress tensor (normal)\n")
            fp.write("#   i   two-point covariance of u component in x- and y-directions\n")
            fp.write("TurbStats " + " ".join(cnfg_so['TurbStats']) + "\n")
            fp.write("\n")
            fp.write("# Pass in mean velocity profile to compute turbulent fluctuations in the\n")
            fp.write("# z-direction.\n")
            try:
                mp_u = cnfg_so['MeanProfile_U']
                mp_v = cnfg_so['MeanProfile_V']
                mp_w = cnfg_so['MeanProfile_W']
                if mp_u is not None and mp_v is not None and mp_w is not None:
                    fp.write("MeanProfile_U   " + mp_u + "\n")
                    fp.write("MeanProfile_V   " + mp_v + "\n")
                    fp.write("MeanProfile_W   " + mp_w + "\n")
            except KeyError:
                pass
            try:
                fp.write("MeanProfile_P   " + cnfg_so['MeanProfile_P'] + "\n")
            except (KeyError, TypeError):
                pass
            try:
                fp.write("MeanProfile_Phi " + cnfg_so['MeanProfile_Phi'] + "\n")
            except (KeyError, TypeError):
                pass
            fp.write("\n")

            # Temperature Equation Parameters
            if cnfg_main['Temperature']:
                try:
                    cnfg_tp = self.cnfg['TP']

                    fp.write("# -------------------------------------------------------------------------------\n")
                    fp.write("#                              Temperature\n")
                    fp.write("# -------------------------------------------------------------------------------\n")
                    fp.write("\n")
                    fp.write("# Notes:\n")
                    fp.write("#    Re = (ReferenceVelocity * L / Nu)\n")
                    fp.write("#    Ra = (g * Beta * (Tmax - Tmin) * L^3) / (Gamma * Nu)\n")
                    fp.write("#    Pr = (Nu * RhoInf) / Gamma\n")
                    fp.write("#    Gr = Ra / Pr\n")
                    fp.write("\n")
                    #fp.write("# Gravity magnitude\n")
                    #fp.write("Gravity 9.801\n")
                    fp.write("# Thermal expansion coefficient (1/T for ideal gas)\n")
                    fp.write("Beta " +  cnfg_tp['Beta'] + "\n")
                    #fp.write("# Typically 1.0.\n")
                    #fp.write("# Setting this to 0 will mean temperature will not drive momentum.\n")
                    #fp.write("Rho_Infinity 1\n")
                    fp.write("# Thermal Diffusivity  (Nu / Prandtl number)\n")
                    fp.write("Gamma " + cnfg_tp['Gamma'] + "\n")
                    fp.write("\n")
                    fp.write("# Scalar transport for temperature\n")
                    fp.write("# Note that gravity is assumed to act in the negative Z direction.\n")
                    fp.write("# The Boussinesq approximation is included in the w-momentum equation.\n")
                    fp.write("Temp_West   " + " ".join(cnfg_tp['Temp_West']) + "\n")
                    fp.write("Temp_East   " + " ".join(cnfg_tp['Temp_East']) + "\n")
                    fp.write("Temp_South  " + " ".join(cnfg_tp['Temp_South']) + "\n")
                    fp.write("Temp_North  " + " ".join(cnfg_tp['Temp_North']) + "\n")
                    fp.write("Temp_Bottom " + " ".join(cnfg_tp['Temp_Bottom']) + "\n")
                    fp.write("Temp_Top    " + " ".join(cnfg_tp['Temp_Top']) + "\n")
                    fp.write("\n")
                    fp.write("# Source term for temperature\n")
                    fp.write("Source_PHI " + cnfg_tp['Source_PHI']+ "\n")
                    fp.write("\n")

                except Exception as err:
                    log(err, 'e')

            # Turbulence Model Parameters
            if cnfg_main['Turbulence']:
                try:
                    pass
                except Exception as err:
                    log(err, 'e')

            # Geometry Parameters
            if cnfg_main['SolidGeometry']:
                try:
                    pass
                except Exception as err:
                    log(err, 'e')

            fp.close()
        except IndexError as err:
            log(err, 'e')
            # display(err, 'e',)
