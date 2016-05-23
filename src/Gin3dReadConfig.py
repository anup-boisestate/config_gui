#!/bin/python3
import sys
#from Gin3dGUI import Gin3dGUI

try:
    from gi.repository import Gtk
except:
    print('Gtk not available')
    sys.exit(1)

from FileChooserWindow import FileChooserWindow


class Gin3dReadConfig:
    def __init__(self, pWin):
        self.cnfg = {}
        self.win = pWin

    def on_filemenu_open_activate(self, obj, data=None):
        filters = {'Gin3D Config': '*.cfg'}
        filename = FileChooserWindow().file_chooser(self.win, filters)

        try:
            self.parse_file(filename)
        except FileExistsError as err:
            print(err)

        gladefile = "../GIN3D_ConfigGUI.glade"
        g3dGUI = Gin3dGUI(gladefile, self.cnfg)
        g3dGUI.show_gui()

        gtk = g3dGUI.get_gtk()
        gtk.main()

    def parse_file(self, filename):

        if filename is None:
            return

        lineNo = 0;
        self.cnfg = {}
        with open(filename, "r") as fp:
            for line in fp:
                lineNo += 1
                l = line.strip()
                if not (l.startswith('#') or '' == l):
                    split = l.split(' ', 1)
                    key = split[0]
                    val = ''
                    if len(split) == 2:
                        val = split[1]
                    self.cnfg[key] = val
                    # else:
                    #   raise FileExistsError('Line No.: %d - %s' %(lineNo, line))
        return
