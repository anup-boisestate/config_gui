#!/bin/python3
import sys

try:
    from gi.repository import Gtk
except:
    print('Gtk not available')
    sys.exit(1)


class FileChooserWindow:
    @staticmethod
    def folder_chooser():
        return

    @staticmethod
    def file_chooser(win, filters=None):
        dialog = Gtk.FileChooserDialog("Please choose a file", win, Gtk.FileChooserAction.OPEN,
                                       (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                        Gtk.STOCK_OPEN, Gtk.ResponseType.OK))

        if filters is not None:
            FileChooserWindow.add_filters(dialog, filters)

        response = dialog.run()
        filename = None
        if response == Gtk.ResponseType.OK:
            #print("Open clicked")
            filename = dialog.get_filename()
        #selif response == Gtk.ResponseType.CANCEL:
            #print("Cancel clicked")

        dialog.destroy()

        return filename

    @staticmethod
    def add_filters(dialog, filters):
        for k in filters:
            f = Gtk.FileFilter()
            f.set_name(k)
            f.add_pattern(filters[k])
            dialog.add_filter(f)



