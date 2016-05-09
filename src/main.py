
from Gin3dGUI import Gin3dGUI


# Create an instance of the Gin3dGUI object
# and run the GUI

# path to the glade file
gladefile = "../GIN3D_ConfigGUI.glade"
g3dGUI = Gin3dGUI(gladefile)
g3dGUI.show_gui()

gtk = g3dGUI.get_gtk()
gtk.main()



