from utilities import *
from writer import *
import re
import gi
from gi.repository import Gtk
from gi.repository import Gdk

gi.require_version('Gtk', '3.0')


# ------------ Helper Functions ---------------

def isInt(data):
    try:
        value = int(data)
        return True
    except ValueError:
        pass
    return False


def isFloat(data):
    try:
        value = float(eval(str(data)))
        return True
    except:
        pass
    return False


def isPositive(data):
    try:
        value = eval(str(data))
        if value < 0:
            raise ValueError('Negative value')
        return True
    except:
        pass
    return False


def text(obj_id, evaluate=False):
    obj = Handler.gtkBuilder.get_object(obj_id)
    if obj is not None:
        if isinstance(obj, Gtk.ComboBoxText):
            return obj.get_active_id()
        else:
            val = obj.get_text()
            if val != '' and evaluate:
                return str(eval(val))
            return val
    return ''


def get_obj(obj_id):
    return Handler.gtkBuilder.get_object(obj_id)


# ---------------------------------------------

class WindowHandler:
    @staticmethod
    def on_MainWindow_destroy(*args):
        log('GIN3D config window closing...', 'd')
        Gtk.main_quit(*args)


class FileMenuHandler:
    @staticmethod
    def on_filemenu_quit_activate(*args):
        log('GIN3D config window closing...', 'd')
        Gtk.main_quit(*args)

    @staticmethod
    def on_filemenu_open_activate(obj, data=None):
        filters = {'Gin3D Config': '*.cfg'}
        filename = file_chooser(Handler.gtkConfigWindow, filters)

        if filename is not None:
            from GIN3DConfigWin import GIN3DConfigWin
            GIN3DConfigWin(filename)

    @staticmethod
    def on_filemenu_saveas_activate(obj, data=None):
        filters = {'Gin3D Config': '*.cfg'}
        filename = file_saveas(Handler.gtkConfigWindow, filters)
        GIN3DConfigWriter(Handler.gtkBuilder).to_file(filename)


class MainParametersHandler:
    # Switch callback functions
    def on_simParam_switch_temperature_toggle(self, obj, data):
        self.toggle_tab(5, obj.get_active(), obj.get_ancestor(Gtk.Notebook))

    def on_simParam_switch_turbulence_toggle(self, obj, data):
        self.toggle_tab(6, obj.get_active(), obj.get_ancestor(Gtk.Notebook))

    def on_simParam_switch_geometry_toggle(self, obj, data):
        self.toggle_tab(7, obj.get_active(), obj.get_ancestor(Gtk.Notebook))

    @staticmethod
    def toggle_tab(tabNo, state, notebook):

        if not state:
            notebook.get_nth_page(tabNo).show()
        else:
            notebook.get_nth_page(tabNo).hide()

    # Forcing callback functions

    def toggle_entry_constPresGrad(self, state):
        fx = get_obj("simParam_txtbox_fx")
        fy = get_obj("simParam_txtbox_fy")
        fz = get_obj("simParam_txtbox_fz")

        fx.set_sensitive(state)
        fy.set_sensitive(state)
        fz.set_sensitive(state)

    def delete_forcing_region_entries(self, n):

        parentBox = get_obj("simParam_forcing_cmfr_box")

        childList = parentBox.get_children()
        noOfChildren = len(childList)
        for i in range(1, n + 1):
            childList[noOfChildren - i].destroy()

    def create_forcing_region_entries(self, n, start):

        parentBox = get_obj("simParam_forcing_cmfr_box")

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
            entry1.connect("focus-out-event", Handler.is_float_positive)
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

            entry2.connect("focus-out-event", Handler.is_float_positive)
            entry3.connect("focus-out-event", Handler.is_float_positive)
            entry4.connect("focus-out-event", Handler.is_float_positive)

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

            if not isPositive(noOfForcingRegions) or noOfForcingRegions == 0:
                noOfForcingRegions = 1
                raise ValueError
        except ValueError:
            obj.set_text(str(noOfForcingRegions))

        parentBox = get_obj("simParam_forcing_cmfr_box")

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
        obj_fr = get_obj("forcing_region")
        if obj_fr is None:
            obj_entry_fr = get_obj("no_of_forcing_regions")
            if obj_entry_fr is not None:
                self.on_forcingregion_changed(obj_entry_fr)
            else:
                print('Error: Object Not Found - "no_of_forcing_regions"')

    def toggle_entry_constMassFlowRate(self, state):
        if state:

            get_obj("simParam_forcing_cmfr_box").show()
            self.activate_forcing_region()
        else:
            get_obj("simParam_forcing_cmfr_box").hide()

    def on_toggled_forcing(self, obj):
        state = obj.get_active()
        radioButton = obj.get_label()

        if radioButton == "Constant Pressure Gradient":
            self.toggle_entry_constPresGrad(state)
        elif radioButton == "Constant Mass Flow Rate":
            self.toggle_entry_constMassFlowRate(state)


class BoundaryConditionsHandler:
    def on_face_value_change(self, obj, data=None):
        if obj is not None:
            if obj.get_active_text() == "inlet":
                self.show_inlet_profile()
            else:
                self.check_face_value_for_inlet()

    @staticmethod
    def show_inlet_profile():
        inletBox = get_obj("inletProfile_box")
        if inletBox is not None:
            if not inletBox.is_visible():
                inletBox.show()
                # enable 'Same as Inlet Profile' radio button in Initial Conditions
                get_obj("ic_rb_sameAsInletProfile").set_sensitive(True)

    @staticmethod
    def check_face_value_for_inlet():
        inletBox = get_obj("inletProfile_box")
        if inletBox is not None:
            if inletBox.is_visible():
                faceTypeList = ["west", "east", "north", "south", "top", "bottom"]
                for each in faceTypeList:
                    obj = get_obj(each + "_face")
                    if obj is not None:
                        if obj.get_active_text() == "inlet":
                            return

                inletBox.hide()
                # enable 'Same as Inlet Profile' radio button in Initial Conditions
                get_obj("ic_rb_sameAsInletProfile").set_sensitive(False)
                get_obj("ic_rb_zeroField").set_active(True)

    @staticmethod
    def on_toggle_variable_profile_for_inlet(obj, data=None):
        vprofile_grid = get_obj("variable_profile_inlet_grid")
        if vprofile_grid is not None and obj is not None:
            state = not obj.get_active()
            if state:
                vprofile_grid.show()
            else:
                vprofile_grid.hide()

    @staticmethod
    def on_toggle_rb_variable_profile_for_inlet(obj):
        if obj is None:
            return

        state = obj.get_active()
        radioButton = obj.get_label()

        if radioButton == "Rough Wall":
            # Always disabled - uncomment to allow users to change Kappa
            # rough_k = get_obj("inletProfile_rough_k");
            # if(rough_k != None):
            #    rough_k.set_sensitive(state)

            rough_z_o = get_obj("inletProfile_rough_z_o")
            if rough_z_o is not None:
                rough_z_o.set_sensitive(state)
                # Uncomment to enable options for specific type
                # elif radioButton == "Smooth Wall":
                # Always disabled - uncomment to allow users to change Kappa
                # smooth_k = get_obj("inletProfile_smooth_k");
                # if(smooth_k != None):
                #    smooth_k.set_sensitive(state)
                # Always disabled - uncomment to allow users to change b
                # smooth_b = get_obj("inletProfile_smooth_b");
                # if(smooth_b != None):
                #    smooth_b.set_sensitive(state)

                # elif radioButton == "Power Law":
                # Always disabled - uncomment to allow users to change n
                # power_n = get_obj("inletProfile_power_n");
                # if(power_n != None):
                #    power_n.set_sensitive(state)

                # elif radioButton == "Parabolic":
                # Always disabled - uncomment to allow users to change a,b,c
                # parabolic_a = get_obj("inletProfile_parabolic_a");
                # if(parabolic_a != None):
                #    parabolic_a.set_sensitive(state)

                # parabolic_b = get_obj("inletProfile_parabolic_b");
                # if(parabolic_b != None):
                #    parabolic_b.set_sensitive(state)

                # parabolic_c = get_obj("inletProfile_parabolic_c");
                # if(parabolic_c != None):
                #    parabolic_c.set_sensitive(state)

    @staticmethod
    def on_toggle_turbulent_for_inlet(obj, data=None):
        turbulent_grid = get_obj("turbulent_inlet_grid")
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
        if isInt(val) and isPositive(val):
            val = int(val)
        else:
            return

        nx_text = text("mesh_nx")
        ny_text = text("mesh_ny")
        nx = 50 if nx_text == "" else int(nx_text)
        ny = 50 if ny_text == "" else int(ny_text)
        obj.set_text(str(val) if val < min(nx, ny) else "")


class InitialConditionsHandler:
    @staticmethod
    def on_toggle_rb_ic(obj):
        if obj is None:
            return

        state = obj.get_active()
        radioButton = obj.get_label()

        if radioButton == "Uniform":
            u = get_obj("ic_uniform_u")
            if u is not None:
                u.set_sensitive(state)
            v = get_obj("ic_uniform_v")
            if v is not None:
                v.set_sensitive(state)
            w = get_obj("ic_uniform_w")
            if w is not None:
                w.set_sensitive(state)

        elif radioButton == "Special Cases":
            sc = get_obj("ic_list_specialCases")
            if sc is not None:
                sc.set_sensitive(state)
                detailBox = get_obj("ic_text_specialCases_detail")
                # detailBox.set_editable(False)
                if detailBox is not None:
                    if not state:
                        detailBox.hide()
                    else:
                        detailBox.show()

    @staticmethod
    def on_changed_ic_specialCases(obj, data=None):
        sc = get_obj("ic_list_specialCases")
        detailBox = get_obj("ic_text_specialCases_detail")
        if sc is not None and detailBox is not None:
            a_id = sc.get_active_id()

            # if a_id == "":
            #    detailBox.get_buffer().set_text(re.sub(' +', ' ',
            #                                           """ u = u_tau ( 1/kappa * log(z+) + 5.2)
            #                                           v = 0.0
            #                                           w = 0.0"""))
            if a_id == "1":
                detailBox.get_buffer().set_text(re.sub(' +', ' ',
                                                       """ Note: Give FrictionVelocity and ReferenceVelocity tags
                                                       u = u_tau ( 1/kappa * log(z+) + 5.2) + 2.0 * C/wx * cos(wx*x)*sin(wy*y)*sin(wz*z)
                                                       v = -C/wy * sin(wx*x)*cos(wy*y)*sin(wz*z)
                                                       w = -C/wz * sin(wx*x)*sin(wy*y)*cos(wz*z)
                                                       wx = round(0.5*LX) * Pi / LX
                                                       wy = round(0.5*LY) * Pi / LY
                                                       wz = round(3.0*LZ) * Pi / LZ
                                                       C = 0.15 * Vref"""))
            # elif a_id == "":
            #    detailBox.get_buffer().set_text(re.sub(' +', ' ',
            #                                           """ Note: theta determined from InletVelocity u and v components
            #                                           u = u_tau/kappa * log(z/z_not) * cos(theta)
            #                                           v = u_tau/kappa * log(z/z_not) * sin(theta)
            #                                           w = 0.0"""))
            elif a_id == "2":
                detailBox.get_buffer().set_text(re.sub(' +', ' ',
                                                       """ Note: theta determined from InletVelocity u and v components
                                                       u = u_tau/kappa * log(z/z_not) * cos(theta) + 2.0*cos(x)*sin(y)*sin(z)
                                                       v = u_tau/kappa * log(z/z_not) * sin(theta) - sin(x)*cos(y)*sin(z)
                                                       w = -sin(x)*sin(y)*( 1.0 + cos(z) )"""))
            elif a_id == "3":
                detailBox.get_buffer().set_text(re.sub(' +', ' ',
                                                       """ u = cos(y) + sin(z)
                                                       v = sin(x) + cos(z)
                                                       w = cos(x) + sin(y)"""))


class FlowSolverParametersHandler:
    @staticmethod
    def on_toggled_rb_momentum(obj):

        id = Gtk.Buildable.get_name(obj)
        state = obj.get_active()

        if id == "fsp_rb_momentum_kappa":
            kappa_text = get_obj("fsp_text_momentum_kappa")
            if kappa_text is not None:
                kappa_text.set_sensitive(state)

        elif id == "fsp_rb_momentum_hybrid":
            hybrid_text = get_obj("fsp_text_momentum_hybrid")
            if hybrid_text is not None:
                hybrid_text.set_sensitive(state)

    @staticmethod
    def on_toggled_rb_temperature(obj):

        id = Gtk.Buildable.get_name(obj)
        state = obj.get_active()

        if id == "fsp_rb_temperature_kappa":
            kappa_text = get_obj("fsp_text_temperature_kappa")
            if kappa_text is not None:
                kappa_text.set_sensitive(state)

        elif id == "fsp_rb_temperature_hybrid":
            hybrid_text = get_obj("fsp_text_temperature_hybrid")
            if hybrid_text is not None:
                hybrid_text.set_sensitive(state)

    @staticmethod
    def on_changed_timeStepSize(obj, data=None):

        timeStepSize = obj.get_active_text()
        variable_grid = get_obj("fsp_grid_timeStepSize_variable")
        constant_grid = get_obj("fsp_grid_timeStepSize_constant")

        if variable_grid is not None and constant_grid is not None:

            if timeStepSize == "Variable":
                variable_grid.show()
                constant_grid.hide()

            elif timeStepSize == "Constant":
                variable_grid.hide()
                constant_grid.show()

    @staticmethod
    def on_changed_poissonSolver(obj, data=None):

        poissonSolver = obj.get_active_text()
        multigrid_grid = get_obj("fsp_grid_poissonSolver_multigrid")
        pointJacobi_grid = get_obj("fsp_grid_poissonSolver_pointJacobi")

        if multigrid_grid is not None and pointJacobi_grid is not None:

            if poissonSolver == "Multigrid":
                multigrid_grid.show()
                pointJacobi_grid.hide()

            elif poissonSolver == "Point Jacobi":
                multigrid_grid.hide()
                pointJacobi_grid.show()

    @staticmethod
    def on_entry_out_cflLimit(obj, event):
        data = obj.get_text()
        cfl = -1
        try:
            cfl = float(data)
        except ValueError:
            pass

        if cfl <= 0 or cfl >= 4:
            obj.set_text("")

    @staticmethod
    def on_entry_out_viscousLimit(obj, event):
        data = obj.get_text()
        viscousLimit = -1
        try:
            viscousLimit = float(data)
        except ValueError:
            pass

        if viscousLimit <= 0 or viscousLimit >= 1:
            obj.set_text("")

    @staticmethod
    def on_entry_out_iterationsLimit(obj, event):
        data = obj.get_text()
        iterationsLimit = -1
        try:
            iterationsLimit = int(data)
        except ValueError:
            pass

        if iterationsLimit < 1:
            obj.set_text("")

    @staticmethod
    def on_entry_out_hybridLimit(obj, event):
        data = obj.get_text()
        hybridLimit = -1
        try:
            hybridLimit = float(data)
        except ValueError:
            pass

        if hybridLimit < 0 or hybridLimit > 1:
            obj.set_text("")


class SolutionOutputHandler:

    def on_toggle_switch_restartSolution(self, obj, data):
        state = not obj.get_active()
        folderChooser = get_obj("so_folder_restartDirectory")
        if folderChooser is not None:
            folderChooser.set_sensitive(state)

    @staticmethod
    def on_toggle_cb_writeRestartFiles(obj, data=None):
        velocity_cb = get_obj("so_cb_outputTXT_u")
        pressure_cb = get_obj("so_cb_outputTXT_p")
        temperature_cb = get_obj("so_cb_outputTXT_t")
        state = obj.get_active()
        file_format_selected = get_obj('so_cb_rawData_binary').get_active() or \
                               get_obj('so_cb_rawData_parallelBinary').get_active() or \
                               get_obj('so_cb_rawData_ASCII').get_active() or \
                               get_obj('so_cb_rawData_parallelASCII').get_active()

        if state:
            velocity_cb.set_active(state)
            pressure_cb.set_active(state)
            temperature_cb.set_active(state)

        if file_format_selected:
            velocity_cb.set_sensitive(not state)
            pressure_cb.set_sensitive(not state)

    @staticmethod
    def on_button_clicked_rawData_selectAll(obj, data=None):
        grid = get_obj("so_grid_rawData_outputVariables")
        if grid is not None:
            gList = grid.get_children()
            for var in gList:
                var.set_active(True)

    @staticmethod
    def on_button_clicked_rawData_deselectAll(obj, data=None):
        grid = get_obj("so_grid_rawData_outputVariables")
        if grid is not None:
            gList = grid.get_children()
            for var in gList:
                if var.get_sensitive():
                    var.set_active(False)

    @staticmethod
    def on_button_clicked_formatOutput_selectAll(obj, data=None):
        grid = get_obj("so_grid_formatOutput_outputVariables")
        if grid is not None:
            gList = grid.get_children()
            for var in gList:
                var.set_active(True)

    @staticmethod
    def on_button_clicked_formatOutput_deselectAll(obj, data=None):
        grid = get_obj("so_grid_formatOutput_outputVariables")
        if grid is not None:
            gList = grid.get_children()
            for var in gList:
                if var.get_sensitive():
                    var.set_active(False)

    @staticmethod
    def on_toggle_cb_formatOutput_timeStep(obj, data=None):
        state = obj.get_active()
        # get_obj("so_text_formatOutput_timeStep_start").set_sensitive(state)
        get_obj("so_text_formatOutput_timeStep_interval1").set_sensitive(state)
        # get_obj("so_text_formatOutput_timeStep_end").set_sensitive(state)
        get_obj("so_cb_formatOutput_timeStep_switch").set_sensitive(state)
        if not state:
            get_obj("so_cb_formatOutput_timeStep_switch").set_active(state)

    @staticmethod
    def on_toggle_cb_formatOutput_timeStep_switch(obj, data=None):
        state = obj.get_active()
        get_obj("so_text_formatOutput_timeStep_switch").set_sensitive(state)
        get_obj("so_text_formatOutput_timeStep_interval2").set_sensitive(state)

    @staticmethod
    def on_toggle_cb_formatOutput_physicalTime(obj, data=None):
        state = obj.get_active()
        # get_obj("so_text_formatOutput_physicalTime_start").set_sensitive(state)
        get_obj("so_text_formatOutput_physicalTime_interval1").set_sensitive(state)
        # get_obj("so_text_formatOutput_physicalTime_end").set_sensitive(state)
        get_obj("so_cb_formatOutput_physicalTime_switch").set_sensitive(state)
        if not state:
            get_obj("so_cb_formatOutput_physicalTime_switch").set_active(state)

    @staticmethod
    def on_toggle_cb_formatOutput_physicalTime_switch(obj, data=None):
        state = obj.get_active()
        get_obj("so_text_formatOutput_physicalTime_switch").set_sensitive(state)
        get_obj("so_text_formatOutput_physicalTime_interval2").set_sensitive(state)

    @staticmethod
    def on_toggle_cb_rawData_timeStep(obj, data=None):
        state = obj.get_active()
        # get_obj("so_text_rawData_timeStep_start").set_sensitive(state)
        get_obj("so_text_rawData_timeStep_interval1").set_sensitive(state)
        # get_obj("so_text_rawData_timeStep_end").set_sensitive(state)
        get_obj("so_cb_rawData_timeStep_switch").set_sensitive(state)
        if not state:
            get_obj("so_cb_rawData_timeStep_switch").set_active(state)

    @staticmethod
    def on_toggle_cb_rawData_timeStep_switch(obj, data=None):
        state = obj.get_active()
        get_obj("so_text_rawData_timeStep_switch").set_sensitive(state)
        get_obj("so_text_rawData_timeStep_interval2").set_sensitive(state)

    @staticmethod
    def on_toggle_cb_rawData_physicalTime(obj, data=None):
        state = obj.get_active()
        # get_obj("so_text_rawData_physicalTime_start").set_sensitive(state)
        get_obj("so_text_rawData_physicalTime_interval1").set_sensitive(state)
        # get_obj("so_text_rawData_physicalTime_end").set_sensitive(state)
        get_obj("so_cb_rawData_physicalTime_switch").set_sensitive(state)
        if not state:
            get_obj("so_cb_rawData_physicalTime_switch").set_active(state)

    @staticmethod
    def on_toggle_cb_rawData_physicalTime_switch(obj, data=None):
        state = obj.get_active()
        get_obj("so_text_rawData_physicalTime_switch").set_sensitive(state)
        get_obj("so_text_rawData_physicalTime_interval2").set_sensitive(state)

    @staticmethod
    def on_toggle_cb_rawData_binary(obj, data=None):
        state = obj.get_active()
        if state and obj.get_label() == "Binary":
            get_obj("so_cb_rawData_parallelBinary").set_active(False)
        elif state and obj.get_label() == "Parallel Binary":
            get_obj("so_cb_rawData_binary").set_active(False)

    @staticmethod
    def on_toggle_cb_rawData_ascii(obj, data=None):
        state = obj.get_active()
        if state and obj.get_label() == "ASCII Dump":
            get_obj("so_cb_rawData_parallelASCII").set_active(False)
        elif state and obj.get_label() == "Parallel ASCII Dump":
            get_obj("so_cb_rawData_ASCII").set_active(False)

    @staticmethod
    def on_out_text_rawData_timeStep_switch(obj, data=None):
        try:
            switchAt = int(obj.get_text())
            #start = int(get_obj("so_text_rawData_timeStep_start").get_text())
            #if switchAt <= start:
            #    obj.set_text("")
        except:
            obj.set_text("")

    @staticmethod
    def on_out_text_rawData_timeStep_end(obj, data=None):
        try:
            end = int(obj.get_text())
            start = int(get_obj("so_text_rawData_timeStep_start").get_text())
            if end <= start:
                obj.set_text("")
        except:
            obj.set_text("")

    @staticmethod
    def on_out_text_rawData_physicalTime_switch(obj, data=None):
        try:
            switchAt = int(obj.get_text())
            #start = int(get_obj("so_text_rawData_timeStep_start").get_text())
            #if switchAt <= start:
            #   obj.set_text("")
        except:
            obj.set_text("")

    @staticmethod
    def on_out_text_rawData_physicalTime_end(obj, data=None):
        try:
            end = int(obj.get_text())
            start = int(get_obj("so_text_rawData_physicalTime_start").get_text())
            if end <= start:
                obj.set_text("")
        except:
            obj.set_text("")

    @staticmethod
    def on_out_text_formatOutput_timeStep_switch(obj, data=None):
        try:
            switchAt = int(obj.get_text())
            #start = int(get_obj("so_text_formatOutput_timeStep_start").get_text())
            #if switchAt <= start:
             #   obj.set_text("")
        except:
            obj.set_text("")

    @staticmethod
    def on_out_text_formatOutput_timeStep_end(obj, data=None):
        try:
            end = int(obj.get_text())
            start = int(get_obj("so_text_formatOutput_timeStep_start").get_text())
            if end <= start:
                obj.set_text("")
        except:
            obj.set_text("")

    @staticmethod
    def on_out_text_formatOutput_physicalTime_switch(obj, data=None):
        try:
            switchAt = int(obj.get_text())
            #start = int(get_obj("so_text_formatOutput_timeStep_start").get_text())
            #if switchAt <= start:
            #    obj.set_text("")
        except:
            obj.set_text("")

    @staticmethod
    def on_out_text_formatOutput_physicalTime_end(obj, data=None):
        try:
            end = int(obj.get_text())
            start = int(get_obj("so_text_formatOutput_physicalTime_start").get_text())
            if end <= start:
                obj.set_text("")
        except:
            obj.set_text("")

    @staticmethod
    def on_toggle_rb_screenDump_outputPlane(obj, data=None):
        state = obj.get_active()
        label = obj.get_label()
        if label == "X-Y":
            get_obj("so_text_screenDump_xy").set_sensitive(state)
        elif label == "X-Z":
            get_obj("so_text_screenDump_xz").set_sensitive(state)
        elif label == "Y-Z":
            get_obj("so_text_screenDump_yz").set_sensitive(state)

    @staticmethod
    def on_out_text_screenDump_outputPlane(obj, data=None):
        try:
            val = float(obj.get_text())
            if val < 0 or val > 1:
                obj.set_text("")
        except:
            obj.set_text("")

    @staticmethod
    def on_toggle_cb_mean_velocity(obj, data=None):
        state = obj.get_active()
        get_obj('so_file_meanprofile_u').set_sensitive(state)
        get_obj('so_file_meanprofile_v').set_sensitive(state)
        get_obj('so_file_meanprofile_w').set_sensitive(state)

    @staticmethod
    def on_toggle_cb_mean_temperature(obj, data=None):
        state = obj.get_active()
        get_obj('so_file_meanprofile_phi').set_sensitive(state)

    @staticmethod
    def on_toggle_cb_mean_pressure(obj, data=None):
        state = obj.get_active()
        get_obj('so_file_meanprofile_p').set_sensitive(state)

    @staticmethod
    def toggle_outputVariables_state(id, state):
        grid = get_obj(id)
        if grid is not None:
            gList = grid.get_children()
            for var in gList:
                var.set_sensitive(state)

    @staticmethod
    def on_toggle_cb_formatOutput(obj, data=None):
        timeStep = get_obj('so_cb_formatOutput_timeStep')
        phyTime = get_obj('so_cb_formatOutput_physicalTime')
        if get_obj('so_cb_fileFormat_hdf5').get_active() or get_obj('so_cb_fileFormat_vti').get_active():
            timeStep.set_sensitive(True)
            phyTime.set_sensitive(True)
            SolutionOutputHandler.toggle_outputVariables_state('so_grid_formatOutput_outputVariables', True)
        else:
            timeStep.set_active(False)
            timeStep.set_sensitive(False)
            timeStep.emit('toggled')

            phyTime.set_active(False)
            phyTime.set_sensitive(False)
            phyTime.emit('toggled')

            SolutionOutputHandler.toggle_outputVariables_state('so_grid_formatOutput_outputVariables', False)

    @staticmethod
    def on_toggle_cb_rawData(obj, data=None):
        timeStep = get_obj('so_cb_rawData_timeStep')
        phyTime = get_obj('so_cb_rawData_physicalTime')
        if get_obj('so_cb_rawData_binary').get_active() or get_obj('so_cb_rawData_parallelBinary').get_active() \
                or get_obj('so_cb_rawData_ASCII').get_active() or get_obj('so_cb_rawData_parallelASCII').get_active():
            timeStep.set_sensitive(True)
            phyTime.set_sensitive(True)
            SolutionOutputHandler.toggle_outputVariables_state('so_grid_rawData_outputVariables', True)
        else:
            timeStep.set_active(False)
            timeStep.set_sensitive(False)
            timeStep.emit('toggled')

            phyTime.set_active(False)
            phyTime.set_sensitive(False)
            phyTime.emit('toggled')

            SolutionOutputHandler.toggle_outputVariables_state('so_grid_rawData_outputVariables', False)

    @staticmethod
    def on_toggle_rb_order(obj, data=None):
        if obj.get_active() and obj.get_label() == 'Higher Order Statistics':
            get_obj('so_cb_stats_d').set_sensitive(True)
            get_obj('so_cb_stats_s').set_sensitive(True)
            get_obj('so_cb_stats_i').set_sensitive(True)
        else:
            get_obj('so_cb_stats_d').set_sensitive(False)
            get_obj('so_cb_stats_s').set_sensitive(False)
            get_obj('so_cb_stats_i').set_sensitive(False)

class Handler(WindowHandler, FileMenuHandler, MainParametersHandler, BoundaryConditionsHandler,
              InitialConditionsHandler, FlowSolverParametersHandler, SolutionOutputHandler):
    gtkBuilder = None
    gtkConfigWindow = None

    def __init__(self, builder):
        Handler.gtkBuilder = builder
        Handler.gtkConfigWindow = builder.get_object('MainWindow')

    @staticmethod
    def is_int_positive(obj, event):
        data = obj.get_text()
        if len(data) < 1:
            return
        if not isInt(data) or not isPositive(data):
            obj.set_text('')

    @staticmethod
    def is_float_positive(obj, event):
        data = obj.get_text()
        if len(data) < 1:
            return
        if not isFloat(data) or not isPositive(data):
            obj.set_text('')

    @staticmethod
    def is_int(obj, event):
        data = obj.get_text()
        if len(data) < 1:
            return
        if not isInt(data):
            obj.set_text('')

    @staticmethod
    def is_float(obj, event):
        data = obj.get_text()
        if len(data) < 1:
            return
        if not isFloat(data):
            obj.set_text('')
