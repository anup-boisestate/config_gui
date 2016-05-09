#!/bin/python3
import sys

try:
    from gi.repository import Gtk
    from gi.repository import Gdk
except:
    print('Gtk not available')
    sys.exit(1)

from Gin3dWriteConfig import Gin3dWriteConfig
from Gin3dReadConfig import Gin3dReadConfig


# noinspection PyPep8Naming
class Gin3dGUI:
    def make_custom_gui_changes(self):

        WHITE = Gdk.Color(65535, 65535, 65535)

        # changing background color
        # self.window.modify_bg(0, WHITE) # Main Window

    def __init__(self, gladefile):

        self.builder = Gtk.Builder()
        self.builder.add_from_file(gladefile)
        self.window = self.builder.get_object("MainWindow")

        # making few custom changes to GUI
        self.make_custom_gui_changes()

        # Listen to callback signals
        self.builder.connect_signals(self)

        # Connect other signals to appropriate listeners
        self.specify_listeners()

        # Hide the Turbulence Model, Temperature Equation
        # and Geometry tabs
        notebook = self.builder.get_object("Notebook")
        notebook.get_nth_page(5).hide()
        notebook.get_nth_page(6).hide()
        notebook.get_nth_page(7).hide()

        self.builder.get_object("simParam_forcing_cmfr_box").hide()

    def show_gui(self):
        self.window.show()

    def get_gtk(self):
        return Gtk

    def close_window(self, obj, data, msg):
        print(msg)
        Gtk.main_quit()

    def specify_listeners(self):
        # GIN3D Configuration File Writer
        writer = Gin3dWriteConfig(self.builder, Gtk)
        reader = Gin3dReadConfig(self.window)

        # filemenu -> saveas
        self.builder.get_object("filemenu_saveas").connect("activate", writer.on_filemenu_saveas_activate)
        self.builder.get_object("filemenu_open").connect("activate", reader.on_filemenu_open_activate)

    ## MainWindow signal callback functions

    def on_MainWindow_destroy(self, obj, data=None):
        self.close_window(obj, data, "Closing application!")

    #####################################

    ## FileMenu signal callback functions

    def on_filemenu_quit_activate(self, filemenu_quit, data=None):
        self.close_window(filemenu_quit, data, "Closing application from file menu!")

    #####################################

    ## Text Entry signal callback and helper functions

    def on_entry_out_isInt(self, obj, event):
        data = obj.get_text()
        if len(data) < 1:
            return
        if not self.isInt(data):
            obj.set_text('')

    def on_entry_out_isInt_isPositive(self, obj, event):
        data = obj.get_text()
        if len(data) < 1:
            return
        if not self.isInt(data) or not self.isPositive(data):
            obj.set_text('')

    def on_entry_out_isFloat(self, obj, event):
        data = obj.get_text()
        if len(data) < 1:
            return
        if not self.isFloat(data):
            obj.set_text('')

    def on_entry_out_isFloat_isPositive(self, obj, event):
        data = obj.get_text()
        if len(data) < 1:
            return
        if not self.isFloat(data) or not self.isPositive(data):
            obj.set_text('')

    def isInt(self, data):
        try:
            value = int(data)
            return True
        except ValueError:
            pass
        return False

    def isFloat(self, data):
        try:
            value = float(eval(str(data)))
            return True
        except:
            pass
        return False

    def isPositive(self, data):
        try:
            value = eval(str(data))
            if value < 0:
                raise ValueError('Negative value')
            return True
        except:
            pass
        return False

    #####################################

    ## Simulation Parameters: Switch callback functions

    def on_simParam_switch_temperature_toggle(self, obj, data):
        self.toggle_tab(5, not obj.get_active())

    def on_simParam_switch_turbulence_toggle(self, obj, data):
        self.toggle_tab(6, not obj.get_active())

    def on_simParam_switch_geometry_toggle(self, obj, data):
        self.toggle_tab(7, not obj.get_active())

    def toggle_tab(self, tabNo, state):
        notebook = self.builder.get_object("Notebook")
        if state:
            notebook.get_nth_page(tabNo).show()
        else:
            notebook.get_nth_page(tabNo).hide()

    #####################################

    ## Simulation Parameters: Forcing callback functions

    def toggle_entry_constPresGrad(self, state):
        fx = self.builder.get_object("simParam_txtbox_fx")
        fy = self.builder.get_object("simParam_txtbox_fy")
        fz = self.builder.get_object("simParam_txtbox_fz")

        fx.set_sensitive(state)
        fy.set_sensitive(state)
        fz.set_sensitive(state)

    def delete_forcing_region_entries(self, n):

        parentBox = self.builder.get_object("simParam_forcing_cmfr_box")

        childList = parentBox.get_children()
        noOfChildren = len(childList)
        for i in range(1, n + 1):
            childList[noOfChildren - i].destroy()

    def create_forcing_region_entries(self, n, start):

        parentBox = self.builder.get_object("simParam_forcing_cmfr_box")

        for i in range(start, n + 1):
            hgrid = Gtk.Grid()
            hgrid.set_column_spacing(2)
            hgrid.set_row_spacing(2)

            label1 = Gtk.Label()
            label1.set_text("Zone #" + str(i) + " Upper Height:")
            label1.set_halign(Gtk.Align.END)
            label1.set_valign(Gtk.Align.CENTER)

            entry1 = Gtk.Entry()
            entry1.set_width_chars(7)
            entry1.connect("focus-out-event", self.on_entry_out_isFloat_isPositive)
            Gtk.Buildable.set_name(entry1, "height_" + str(i))

            hgrid.add(label1)
            hgrid.attach(entry1, 1, 0, 1, 1)

            label2 = Gtk.Label()
            label2.set_text("Zone #" + str(i) + " Mass Flow Rates:")
            label2.set_halign(Gtk.Align.END)
            label2.set_valign(Gtk.Align.CENTER)

            entry2 = Gtk.Entry()
            entry3 = Gtk.Entry()
            entry4 = Gtk.Entry()

            entry2.set_width_chars(7)
            entry3.set_width_chars(7)
            entry4.set_width_chars(7)

            Gtk.Buildable.set_name(entry3, "frfr_y_" + str(i))
            Gtk.Buildable.set_name(entry2, "frfr_x_" + str(i))
            Gtk.Buildable.set_name(entry4, "frfr_z_" + str(i))

            entry2.set_placeholder_text("x")
            entry3.set_placeholder_text("y")
            entry4.set_placeholder_text("z")

            entry4.set_text("0.0")
            entry4.set_sensitive(False)

            entry2.connect("focus-out-event", self.on_entry_out_isFloat_isPositive)
            entry3.connect("focus-out-event", self.on_entry_out_isFloat_isPositive)
            entry4.connect("focus-out-event", self.on_entry_out_isFloat_isPositive)

            hgrid.attach_next_to(label2, label1, Gtk.PositionType.BOTTOM, 1, 1)
            hgrid.attach_next_to(entry2, label2, Gtk.PositionType.RIGHT, 1, 1)
            hgrid.attach_next_to(entry3, entry2, Gtk.PositionType.RIGHT, 1, 1)
            hgrid.attach_next_to(entry4, entry3, Gtk.PositionType.RIGHT, 1, 1)

            parentBox.pack_start(hgrid, True, True, 0)

        parentBox.show_all()

        return

    def on_forcingregion_changed(self, obj, data=None):
        noOfForcingRegions = 1
        try:
            noOfForcingRegions = int(obj.get_text())

            if not self.isPositive(noOfForcingRegions) or noOfForcingRegions == 0:
                noOfForcingRegions = 1
                raise ValueError
        except ValueError:
            obj.set_text(str(noOfForcingRegions))

        parentBox = self.builder.get_object("simParam_forcing_cmfr_box")

        # Subtracting 1 because one of the children is there by default
        # which is the line to input "Number of Forcing Regions" in GUI
        prev_noOfForcingRegions = len(parentBox.get_children()) - 1

        if noOfForcingRegions < prev_noOfForcingRegions:
            noOfRegionsToDelete = prev_noOfForcingRegions - noOfForcingRegions
            self.delete_forcing_region_entries(noOfRegionsToDelete)

        elif noOfForcingRegions > prev_noOfForcingRegions:
            start = prev_noOfForcingRegions + 1
            self.create_forcing_region_entries(noOfForcingRegions, start)

        return

    def activate_forcing_region(self):
        obj_fr = self.builder.get_object("forcing_region")
        if obj_fr is None:
            obj_entry_fr = self.builder.get_object("no_of_forcing_regions")
            if obj_entry_fr is not None:
                self.on_forcingregion_changed(obj_entry_fr)
            else:
                print('Error: Object Not Found - "no_of_forcing_regions"')

    def toggle_entry_constMassFlowRate(self, state):
        if state:
            self.builder.get_object("simParam_forcing_cmfr_box").show()
            self.activate_forcing_region()
        else:
            self.builder.get_object("simParam_forcing_cmfr_box").hide()

    def on_toggled_forcing(self, obj):
        state = obj.get_active()
        radioButton = obj.get_label()

        if radioButton == "Constant Pressure Gradient":
            self.toggle_entry_constPresGrad(state)
        elif radioButton == "Constant Mass Flow Rate":
            self.toggle_entry_constMassFlowRate(state)

    #####################################


    #########################################
    # signals form Boundary Conditions Window
    #########################################

    def on_face_value_change(self, obj, data=None):
        if obj is not None:
            if obj.get_active_text() == "inlet":
                self.show_inlet_profile()
            else:
                self.check_face_value_for_inlet()

    def show_inlet_profile(self):
        inletBox = self.builder.get_object("inletProfile_box")
        if inletBox is not None:
            if not inletBox.is_visible():
                inletBox.show()

    def check_face_value_for_inlet(self):
        inletBox = self.builder.get_object("inletProfile_box")
        if inletBox is not None:
            if inletBox.is_visible():
                faceTypeList = ["west", "east", "north", "south", "top", "bottom"]
                for each in faceTypeList:
                    obj = self.builder.get_object(each + "_face")
                    if obj is not None:
                        if obj.get_active_text() == "inlet":
                            return

                inletBox.hide()

    def on_toggle_variable_profile_for_inlet(self, obj, data=None):
        vprofile_grid = self.builder.get_object("variable_profile_inlet_grid")
        if vprofile_grid is not None and obj is not None:
            state = not obj.get_active()
            if state:
                vprofile_grid.show()
            else:
                vprofile_grid.hide()

    def on_toggle_rb_variable_profile_for_inlet(self, obj):
        if obj is None:
            return

        state = obj.get_active()
        radioButton = obj.get_label()

        if radioButton == "Rough Wall":
            # Always disabled - uncomment to allow users to change Kappa
            # rough_k = self.builder.get_object("inletProfile_rough_k");
            # if(rough_k != None):
            #    rough_k.set_sensitive(state)

            rough_z_o = self.builder.get_object("inletProfile_rough_z_o")
            if rough_z_o is not None:
                rough_z_o.set_sensitive(state)
                # Uncomment to enable options for specific type
                # elif radioButton == "Smooth Wall":
                # Always disabled - uncomment to allow users to change Kappa
                # smooth_k = self.builder.get_object("inletProfile_smooth_k");
                # if(smooth_k != None):
                #    smooth_k.set_sensitive(state)
                # Always disabled - uncomment to allow users to change b
                # smooth_b = self.builder.get_object("inletProfile_smooth_b");
                # if(smooth_b != None):
                #    smooth_b.set_sensitive(state)

                # elif radioButton == "Power Law":
                # Always disabled - uncomment to allow users to change n
                # power_n = self.builder.get_object("inletProfile_power_n");
                # if(power_n != None):
                #    power_n.set_sensitive(state)

                # elif radioButton == "Parabolic":
                # Always disabled - uncomment to allow users to change a,b,c
                # parabolic_a = self.builder.get_object("inletProfile_parabolic_a");
                # if(parabolic_a != None):
                #    parabolic_a.set_sensitive(state)

                # parabolic_b = self.builder.get_object("inletProfile_parabolic_b");
                # if(parabolic_b != None):
                #    parabolic_b.set_sensitive(state)

                # parabolic_c = self.builder.get_object("inletProfile_parabolic_c");
                # if(parabolic_c != None):
                #    parabolic_c.set_sensitive(state)

    def on_toggle_turbulent_for_inlet(self, obj, data=None):
        turbulent_grid = self.builder.get_object("turbulent_inlet_grid")
        if turbulent_grid is not None and obj is not None:
            state = not obj.get_active()
            if state:
                turbulent_grid.show()
            else:
                turbulent_grid.hide()

    # checks whether the value for offset from inlet # of
    # mesh points is <= min(NX,NY) or if NX,NY is not
    # given caps the value at 50
    def on_entry_out_check_value(self, obj, data=None):
        val = obj.get_text()
        if self.isInt(val) and self.isPositive(val):
            val = int(val)
        else:
            return

        nx_text = self.builder.get_object("mesh_nx").get_text()
        ny_text = self.builder.get_object("mesh_ny").get_text()
        nx = 50 if nx_text == "" else int(nx_text)
        ny = 50 if ny_text == "" else int(ny_text)
        obj.set_text(str(val) if val < min(nx,ny) else "")
    #########################################


    #########################################
    # signals form Initial Conditions Window
    #########################################

    def on_toggle_rb_ic(self, obj):
        if obj is None:
            return

        state = obj.get_active()
        radioButton = obj.get_label()

        if radioButton == "Uniform":
            u = self.builder.get_object("ic_uniform_u")
            if u is not None:
                u.set_sensitive(state)
            v = self.builder.get_object("ic_uniform_v")
            if v is not None:
                v.set_sensitive(state)
            w = self.builder.get_object("ic_uniform_w")
            if w is not None:
                w.set_sensitive(state)

        elif radioButton == "Special Cases":
            sc = self.builder.get_object("ic_list_specialCases")
            if sc is not None:
                sc.set_sensitive(state)
                detailBox = self.builder.get_object("ic_text_specialCases_detail")
                # detailBox.set_editable(False)
                if detailBox is not None:
                    if not state:
                        detailBox.hide()
                    else:
                        detailBox.show()

    def on_changed_ic_specialCases(self, obj, data=None):
        sc = self.builder.get_object("ic_list_specialCases")
        detailBox = self.builder.get_object("ic_text_specialCases_detail")
        if sc is not None and detailBox is not None:
            a_id = sc.get_active_id()
            if a_id == "1":
                detailBox.get_buffer().set_text('u = u_tau ( 1/kappa * log(z+) + 5.2)\nv = 0.0\nw = 0.0')
            elif a_id == "2":
                detailBox.get_buffer().set_text(
                    "Note: Give FrictionVelocity and ReferenceVelocity tags\nu = u_tau ( 1/kappa * log(z+) + 5.2) + 2.0 * C/wx * cos(wx*x)*sin(wy*y)sin(wz*z)\nv = -C/wy * sin(wx*x)*cos(wy*y)sin(wz*z)\nw = -C/wz * sin(wx*x)*sin(wy*y)cos(wz*z)\nwx = round(0.5*LX) * Pi / LX\nwy = round(0.5*LY) * Pi / LY\nwz = round(3.0*LZ) * Pi / LZ\nC = 0.15 * Vref")
            elif a_id == "3":
                detailBox.get_buffer().set_text(
                    "Note: theta determined from InletVelocity u and v components\nu = u_tau/kappa * log(z/z_not) * cos(theta)\nv = u_tau/kappa * log(z/z_not) * sin(theta)\nw = 0.0")
            elif a_id == "4":
                detailBox.get_buffer().set_text(
                    "Note: theta determined from InletVelocity u and v components\nu = u_tau/kappa * log(z/z_not) * cos(theta) + 2.0*cos(x)*sin(y)*sin(z)\nv = u_tau/kappa * log(z/z_not) * sin(theta) - sin(x)*cos(y)*sin(z)\nw = -sin(x)*sin(y)*( 1.0 + cos(z) )")
            elif a_id == "5":
                detailBox.get_buffer().set_text("u = cos(y) + sin(z)\nv = sin(x) + cos(z)\nw = cos(x) + sin(y)")

    #########################################

    #########################################
    # signals from Flow Solver Parameters Window
    #########################################

    def on_toggled_rb_momentum(self, obj):

        id = Gtk.Buildable.get_name(obj)
        state = obj.get_active()

        if id == "fsp_rb_momentum_kappa":
            kappa_text = self.builder.get_object("fsp_text_momentum_kappa")
            if kappa_text is not None:
                kappa_text.set_sensitive(state)

        elif id == "fsp_rb_momentum_hybrid":
            hybrid_text = self.builder.get_object("fsp_text_momentum_hybrid")
            if hybrid_text is not None:
                hybrid_text.set_sensitive(state)

    def on_toggled_rb_temperature(self, obj):

        id = Gtk.Buildable.get_name(obj)
        state = obj.get_active()

        if id == "fsp_rb_temperature_kappa":
            kappa_text = self.builder.get_object("fsp_text_temperature_kappa")
            if kappa_text is not None:
                kappa_text.set_sensitive(state)

        elif id == "fsp_rb_temperature_hybrid":
            hybrid_text = self.builder.get_object("fsp_text_temperature_hybrid")
            if hybrid_text is not None:
                hybrid_text.set_sensitive(state)

    def on_changed_timeStepSize(self, obj, data=None):

        timeStepSize = obj.get_active_text()
        variable_grid = self.builder.get_object("fsp_grid_timeStepSize_variable")
        constant_grid = self.builder.get_object("fsp_grid_timeStepSize_constant")

        if variable_grid is not None and constant_grid is not None:

            if timeStepSize == "Variable":
                variable_grid.show()
                constant_grid.hide()

            elif timeStepSize == "Constant":
                variable_grid.hide()
                constant_grid.show()

    def on_changed_poissonSolver(self, obj, data=None):

        poissonSolver = obj.get_active_text()
        multigrid_grid = self.builder.get_object("fsp_grid_poissonSolver_multigrid")
        pointJacobi_grid = self.builder.get_object("fsp_grid_poissonSolver_pointJacobi")

        if multigrid_grid is not None and pointJacobi_grid is not None:

            if poissonSolver == "Multigrid":
                multigrid_grid.show()
                pointJacobi_grid.hide()

            elif poissonSolver == "Point Jacobi":
                multigrid_grid.hide()
                pointJacobi_grid.show()

    def on_entry_out_cflLimit(self, obj, event):
        data = obj.get_text()
        cfl = -1
        try:
            cfl = float(data)
        except ValueError:
            pass

        if cfl <= 0 or cfl >= 4:
            obj.set_text("")

    def on_entry_out_viscousLimit(self, obj, event):
        data = obj.get_text()
        viscousLimit = -1
        try:
            viscousLimit = float(data)
        except ValueError:
            pass

        if viscousLimit <= 0 or viscousLimit >= 1:
            obj.set_text("")

    def on_entry_out_iterationsLimit(self, obj, event):
        data = obj.get_text()
        iterationsLimit = -1
        try:
            iterationsLimit = int(data)
        except ValueError:
            pass

        if iterationsLimit < 1:
            obj.set_text("")

    def on_entry_out_hybridLimit(self, obj, event):
        data = obj.get_text()
        hybridLimit = -1
        try:
            hybridLimit = float(data)
        except ValueError:
            pass

        if hybridLimit < 0 or hybridLimit > 1:
            obj.set_text("")

    #########################################

    #########################################
    # Signals form Solution Output Window
    #########################################

    def on_toggle_switch_restartSolution(self, obj, data=None):
        state = not obj.get_active()

        folderChooser = self.builder.get_object("so_folderChooser_restartSolution")
        if folderChooser is not None:
            folderChooser.set_sensitive(state)

    def on_toggle_cb_writeRestartFiles(self, obj, data=None):
        velocity_cb = self.builder.get_object("so_cb_rawData_velocity")
        pressure_cb = self.builder.get_object("so_cb_rawData_pressure")
        temperature_cb = self.builder.get_object("so_cb_rawData_temperature")
        state = obj.get_active()
        if state:
            velocity_cb.set_active(state)
            pressure_cb.set_active(state)
            temperature_cb.set_active(state)

        velocity_cb.set_sensitive(not state)
        pressure_cb.set_sensitive(not state)

    def on_button_clicked_rawData_selectAll(self, obj, data=None):

        grid = self.builder.get_object("so_grid_rawData_outputVariables")
        if grid is not None:
            gList = grid.get_children()
            for var in gList:
                var.set_active(True)

    def on_button_clicked_rawData_deselectAll(self, obj, data=None):

        grid = self.builder.get_object("so_grid_rawData_outputVariables")
        if grid is not None:
            gList = grid.get_children()
            for var in gList:
                if var.get_sensitive():
                    var.set_active(False)

    def on_button_clicked_formatOutput_selectAll(self, obj, data=None):

        grid = self.builder.get_object("so_grid_formatOutput_outputVariables")
        if grid is not None:
            gList = grid.get_children()
            for var in gList:
                var.set_active(True)

    def on_button_clicked_formatOutput_deselectAll(self, obj, data=None):
        grid = self.builder.get_object("so_grid_formatOutput_outputVariables")
        if grid is not None:
            gList = grid.get_children()
            for var in gList:
                if var.get_sensitive():
                    var.set_active(False)

    def on_toggle_cb_formatOutput_timeStep(self, obj, data=None):
        state = obj.get_active()
        self.builder.get_object("so_text_formatOutput_timeStep_start").set_sensitive(state)
        self.builder.get_object("so_text_formatOutput_timeStep_interval1").set_sensitive(state)
        self.builder.get_object("so_text_formatOutput_timeStep_end").set_sensitive(state)
        self.builder.get_object("so_cb_formatOutput_timeStep_switch").set_sensitive(state)
        if not state:
            self.builder.get_object("so_cb_formatOutput_timeStep_switch").set_active(state)

    def on_toggle_cb_formatOutput_timeStep_switch(self, obj, data=None):
        state = obj.get_active()
        self.builder.get_object("so_text_formatOutput_timeStep_switch").set_sensitive(state)
        self.builder.get_object("so_text_formatOutput_timeStep_interval2").set_sensitive(state)

    def on_toggle_cb_formatOutput_physicalTime(self, obj, data=None):
        state = obj.get_active()
        self.builder.get_object("so_text_formatOutput_physicalTime_start").set_sensitive(state)
        self.builder.get_object("so_text_formatOutput_physicalTime_interval1").set_sensitive(state)
        self.builder.get_object("so_text_formatOutput_physicalTime_end").set_sensitive(state)
        self.builder.get_object("so_cb_formatOutput_physicalTime_switch").set_sensitive(state)
        if not state:
            self.builder.get_object("so_cb_formatOutput_physicalTime_switch").set_active(state)

    def on_toggle_cb_formatOutput_physicalTime_switch(self, obj, data=None):
        state = obj.get_active()
        self.builder.get_object("so_text_formatOutput_physicalTime_switch").set_sensitive(state)
        self.builder.get_object("so_text_formatOutput_physicalTime_interval2").set_sensitive(state)

    def on_toggle_cb_rawData_timeStep(self, obj, data=None):
        state = obj.get_active()
        self.builder.get_object("so_text_rawData_timeStep_start").set_sensitive(state)
        self.builder.get_object("so_text_rawData_timeStep_interval1").set_sensitive(state)
        self.builder.get_object("so_text_rawData_timeStep_end").set_sensitive(state)
        self.builder.get_object("so_cb_rawData_timeStep_switch").set_sensitive(state)
        if not state:
            self.builder.get_object("so_cb_rawData_timeStep_switch").set_active(state)

    def on_toggle_cb_rawData_timeStep_switch(self, obj, data=None):
        state = obj.get_active()
        self.builder.get_object("so_text_rawData_timeStep_switch").set_sensitive(state)
        self.builder.get_object("so_text_rawData_timeStep_interval2").set_sensitive(state)

    def on_toggle_cb_rawData_physicalTime(self, obj, data=None):
        state = obj.get_active()
        self.builder.get_object("so_text_rawData_physicalTime_start").set_sensitive(state)
        self.builder.get_object("so_text_rawData_physicalTime_interval1").set_sensitive(state)
        self.builder.get_object("so_text_rawData_physicalTime_end").set_sensitive(state)
        self.builder.get_object("so_cb_rawData_physicalTime_switch").set_sensitive(state)
        if not state:
            self.builder.get_object("so_cb_rawData_physicalTime_switch").set_active(state)

    def on_toggle_cb_rawData_physicalTime_switch(self, obj, data=None):
        state = obj.get_active()
        self.builder.get_object("so_text_rawData_physicalTime_switch").set_sensitive(state)
        self.builder.get_object("so_text_rawData_physicalTime_interval2").set_sensitive(state)

    def on_toggle_cb_rawData_binary(self, obj, data=None):
        state = obj.get_active()
        if state and obj.get_label() == "Binary":
            self.builder.get_object("so_cb_rawData_parallelBinary").set_active(False)
        elif state and obj.get_label() == "Parallel Binary":
            self.builder.get_object("so_cb_rawData_binary").set_active(False)

    def on_toggle_cb_rawData_ascii(self, obj, data=None):
        state = obj.get_active()
        if state and obj.get_label() == "ASCII Dump":
            self.builder.get_object("so_cb_rawData_parallelASCII").set_active(False)
        elif state and obj.get_label() == "Parallel ASCII Dump":
            self.builder.get_object("so_cb_rawData_ASCII").set_active(False)

    def on_out_text_rawData_timeStep_switch(self, obj, data=None):
        try:
            switchAt = int(obj.get_text())
            start = int(self.builder.get_object("so_text_rawData_timeStep_start").get_text())
            if switchAt <= start:
                obj.set_text("")
        except:
            obj.set_text("")

    def on_out_text_rawData_timeStep_end(self, obj, data=None):
        try:
            end = int(obj.get_text())
            start = int(self.builder.get_object("so_text_rawData_timeStep_start").get_text())
            if end <= start:
                obj.set_text("")
        except:
            obj.set_text("")

    def on_out_text_rawData_physicalTime_switch(self, obj, data=None):
        try:
            switchAt = int(obj.get_text())
            start = int(self.builder.get_object("so_text_rawData_timeStep_start").get_text())
            if switchAt <= start:
                obj.set_text("")
        except:
            obj.set_text("")

    def on_out_text_rawData_physicalTime_end(self, obj, data=None):
        try:
            end = int(obj.get_text())
            start = int(self.builder.get_object("so_text_rawData_physicalTime_start").get_text())
            if end <= start:
                obj.set_text("")
        except:
            obj.set_text("")

    def on_out_text_formatOutput_timeStep_switch(self, obj, data=None):
        try:
            switchAt = int(obj.get_text())
            start = int(self.builder.get_object("so_text_formatOutput_timeStep_start").get_text())

            if switchAt <= start:
                obj.set_text("")
        except:
            obj.set_text("")

    def on_out_text_formatOutput_timeStep_end(self, obj, data=None):
        try:
            end = int(obj.get_text())
            start = int(self.builder.get_object("so_text_formatOutput_timeStep_start").get_text())
            if end <= start:
                obj.set_text("")
        except:
            obj.set_text("")

    def on_out_text_formatOutput_physicalTime_switch(self, obj, data=None):
        try:
            switchAt = int(obj.get_text())
            start = int(self.builder.get_object("so_text_formatOutput_timeStep_start").get_text())
            if switchAt <= start:
                obj.set_text("")
        except:
            obj.set_text("")

    def on_out_text_formatOutput_physicalTime_end(self, obj, data=None):
        try:
            end = int(obj.get_text())
            start = int(self.builder.get_object("so_text_formatOutput_physicalTime_start").get_text())
            if end <= start:
                obj.set_text("")
        except:
            obj.set_text("")

    def on_toggle_rb_screenDump_outputPlane(self, obj, data=None):
        state = obj.get_active()
        label = obj.get_label()
        if label == "X-Y":
            self.builder.get_object("so_text_screenDump_xy").set_sensitive(state)
        elif label == "X-Z":
            self.builder.get_object("so_text_screenDump_xz").set_sensitive(state)
        elif label == "Y-Z":
            self.builder.get_object("so_text_screenDump_yz").set_sensitive(state)

    def on_out_text_screenDump_outputPlane(self, obj, data=None):
        try:
            val = float(obj.get_text())
            if val < 0 or val > 1:
                obj.set_text("")
        except:
            obj.set_text("")
