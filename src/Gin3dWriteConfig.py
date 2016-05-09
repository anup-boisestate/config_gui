#!/bin/python3
import sys

try:
    from gi.repository import Gtk
    from gi.repository import Gdk
except:
    print('Gtk not available')
    sys.exit(1)


class Gin3dWriteConfig:
    def __init__(self, builder, gtk):
        self.cnfg = {}
        self.win = builder
        self.gtk = gtk

    def on_filemenu_saveas_activate(self, obj, data=None):
        print("Gather all data")
        self.get_data()

        print("Write to file")
        self.writeToFile("template-test")

    def get_data(self):
        self.cnfg['Main'] = self.get_main_parameters()
        self.cnfg['BC'] = self.get_bc_parameters()
        self.cnfg['IC'] = self.get_ic_parameters()
        print(self.cnfg)

    def get_main_parameters(self):
        cnfg_main = {'Grid': "%s %s %s" % (self.text("domain_lx"), self.text("domain_ly"), self.text("domain_lz")),
                     'Dimensions': "%s %s %s" % (self.text("mesh_nx"), self.text("mesh_ny"), self.text("mesh_nz")),
                     'ReferenceLength': self.text("refLength", True), 'Nu': self.text("laminar_viscosity", True),
                     'StopTime': self.text("sim_endTime"), 'Temperature': self.active("simParam_switch_temperature")}
        # cnfg_main['Turbulence'] = self.active("simParam_switch_turbulence")
        # cnfg_main['SolidGeometry'] = self.active("simParam_switch_geometry")

        self.forcing = "off"  # variable set to check when writing to file
        forcing_off = self.active("simParam_rb_forcing_off")

        if not forcing_off:
            forcing_cpg = self.active("simParam_rb_forcing_constPresGrad")
            if forcing_cpg:
                self.forcing = "constPresGrad"
                cnfg_main['Forcing'] = "%s %s %s" % (
                    self.text("simParam_txtbox_fx", True), self.text("simParam_txtbox_fy", True), self.text("simParam_txtbox_fz", True))
            else:
                self.forcing = "constMassFlowRate"
                cnfg_main['NumberForcingRegions'] = self.text("no_of_forcing_regions")
                numforcings = int(cnfg_main['NumberForcingRegions'])
                heights = [0.0]
                mdot = []
                for i in range(numforcings):
                    h = self.fchild("simParam_forcing_cmfr_box", "height_", i)
                    if h is not None:
                        heights.append(eval(h.get_text()))
                        tmp = []
                        x = self.fchild("simParam_forcing_cmfr_box", "frfr_x_", i)
                        if x is not None:
                            tmp.append(str(eval(x.get_text())))

                        y = self.fchild("simParam_forcing_cmfr_box", "frfr_y_", i)
                        if y is not None:
                            tmp.append(str(eval(y.get_text())))

                        z = self.fchild("simParam_forcing_cmfr_box", "frfr_z_", i)
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
                val = "RoughWall"
            elif self.active("rb_inletProfile_smoothWall"):
                val = "SmoothWall"
            elif self.active("rb_inletProfile_powerLaw"):
                val = "PowerLaw"
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
        return cnfg_ic

    def text(self, obj_id, evaluate=False):
        obj = self.win.get_object(obj_id)
        if obj is not None:
            if isinstance(obj, Gtk.ComboBoxText):
                return obj.get_active_id()
            else:
                val = obj.get_text()
                if val != '':
                    return str(eval(val))
                return val
        return ''

    def active(self, obj_id):
        return self.win.get_object(obj_id).get_active()

    def obj(self, obj_id):
        return self.win.get_object(obj_id)

    def fchild(self, p_id, c_id, index):
        p = self.obj(p_id)
        gridList = p.get_children()
        for c in gridList[index + 1].get_children():
            if self.gtk.Buildable.get_name(c) is not None:
                if self.gtk.Buildable.get_name(c) == c_id + str(index + 1):
                    return c

    def writeToFile(self, filename):
        fp = open(filename + ".cfg", "w")

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

        # Item 9
        # Do nothing

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

        fp.close()
