'''
    spiral_overlay.py
    
    DESCRIPTION:        Overlay spirals onto a FITS image.
    
    USAGE:              python spiral_overlay.py <filename>
    
    INPUTS:             Filename.
                        FITS filename can be passed as a command line argument,
                        but the program will prompt for it if it isn't.
    
    OUTPUT:             Interactive plot in MatPlotLib window.
    
    NOTES:              Uses pylab (including matplotlib and numpy) and astropy.
                        All necessary libraries included in Ureka package.
    
                        Eqns used:
    
                            r = ae^(b*theta)
    
                            x = r*cos(theta)
                            y = r*sin(theta)
    
                            tan(phi) = b
                            phi = pitch angle
    
                        
                        It is recommended that you use cropped images with this program. Image
                        plotting slows down dramatically as image size increases.
                    
                        Rotation angle and pitch can be changed while the program is running
                        using the sliders below the plot. For strongly barred galaxies, the minimum 
                        polar angle can be set with the slider marked minimum theta. The reset button 
                        initializes the slider to default values (0.0 for rotation, 20.0 for pitch, and 
                        pi for min. theta).
                        
                        Scale can be swapped between logarithmic and linear modes with the radio
                        buttons on the left side of the plot. There are two colormaps for the logarithmic
                        scale: Seismic (diverging) and Greyscale (sequential). Defaults to linear scale,
                        which you may want to change if you're dealing with faint real images.
                        
                        Chirality (corresponding to the sign of the pitch angle: CCW = -, CW = +)
                        can be changed with the second set of radio buttons on the right side of
                        the plot. Defaults to positive chirality.
                        
                        Arm number (corresponding to mode) can be changed with the radio buttons
                        on the right side of the plot. Defaults to 2 arms.
                        
                        Default values can be changed in the block of constants below.
                        
                        The overlay can be saved as a .png (without all the sliders and buttons) by pressing
                        the save button below the sliders.
'''

#-----------------------------code begins----------------------------
#!/usr/bin/python
from pylab import *
import numpy as np
import glob
from astropy.io import fits
from matplotlib.colors import LogNorm, Normalize
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import sys

#-----------------------------Default/initial values-----------------------------
DEFAULT_ARMS = 2    #Between 1 and 6, inclusive
DEFAULT_PITCH = 20.0  #Between about 0.0 and 45.0 degrees.
DEFAULT_ROTATION = 0.0 #Between -180.0 and 180.0, inclusive
DEFAULT_THETA_MIN = np.pi #Between 0.0 and 3*pi
DEFAULT_COLORSCALE = 'n' #n is linear, yg is logarithmic grey scale, ys is logarithmic seismic scale
line_color = '' #Options: '' for default, 'red', 'black', 'white'

class Spiral():
    '''
        Spiral
        A class interface to hold parameters related to the spiral overlay.
    '''

    def __init__(self, num_arms, pitch, chirality, rotation, theta_min, color):
        self.num_arms = num_arms
        self.pitch = pitch
        self.chirality = chirality
        self.rotation = rotation
        self.theta_min = theta_min
        self.color = color
        return
    #end definition

    def setNumArms(self, num_arms):
        self.num_arms = num_arms
        return
    #end definition

    def setPitch(self, pitch):
        self.pitch = pitch
        return
    #end definition
    
    def setChirality(self, chirality):
        self.chirality = chirality
        return
    #end definition

    def setRotation(self, rotation):
        self.rotation = rotation
        return
    #end definition

    def setTheta(self, theta_min):
        self.theta_min = theta_min
        return
    #end definition

    def setColor(self, color):
        self.color = color
        return
    #end definition

#end class definition

class Image():
    '''
        Image
        A class interface to hold data from and related to the galaxy image being 
        overlaid.
    '''
    def __init__(self, name):
        self.name = name
        return
    #end definition

    def setData(self, data):
        self.data = data
        return
    #end definition

    def setColorscale(self, colorscale):
        self.colorscale = colorscale
        return
    #end definition

#end class definition

#-----------------------------Command line input-----------------------------
if len(sys.argv) != 1:
    if sys.argv[1][-4:] == ".fit" or sys.argv[1][-5:] == ".fits":
        gal_image = Image(sys.argv[1])
    elif sys.argv[1] in ["-h","-H","--help"]:
        print(__doc__)
        sys.exit()
    else:
        print("Could not parse command line argument")
        print(__doc__)
        sys.exit()
    #end if

else:
    print "Available images in this directory: "
    #If you're using python 3.5+ you can easily make this search subdirectories recursively
    for filename in (glob.glob("*.fit") + glob.glob("*.fits")):
        print filename
    gal_image = Image(str(raw_input("\nEnter filename for galaxy image.\n> ")))
#end if

try:
    #pull the pixel values and save them as a uint8 numpy array
    gal_image.setData(fits.getdata(gal_image.name))
except IOError:
    print "\nIOError: File not found. Check cwd or filename and try again"
    sys.exit()

#-----------------------------Setup for plot-----------------------------
fig, ax = plt.subplots()
plt.subplots_adjust(left= 0.15, bottom= 0.20) #set placement of plot
plt.tick_params(axis= 'both', which= 'both', bottom= 'off', top= 'off', left= 'off', right= 'off', labelbottom= 'off', labelleft= 'off')
plt.axis([0, int(len(gal_image.data)), 0, int(len(gal_image.data))]) #initial axis bounds set by this array
fig.text(0.52, 0.92, "Image: "+gal_image.name, fontsize= '14', ha= 'center')

#-----------------------------Image handling-----------------------------
gal_image.setColorscale(DEFAULT_COLORSCALE)

if gal_image.colorscale == 'n':
    colorscale_flag = 0
elif gal_image.colorscale == 'yg':
    colorscale_flag = 1
elif gal_image.colorscale == 'ys':
    colorscale_flag = 2
#end if

#ImagePlot
#   A function to handle image renormalization, color mapping, and display
def ImagePlot(image):
    if str(image.colorscale)=='n':
        remap = Normalize()
        remap.autoscale(image.data)
        ax.imshow(image.data, cmap='gray', norm=remap, origin='lower')
    elif str(image.colorscale) == 'yg':
        remap = LogNorm()
        remap.autoscale(image.data)
        ax.imshow(image.data, cmap='gray', norm=remap, origin='lower')
    elif str(image.colorscale) == 'ys':
        remap = LogNorm()
        remap.autoscale(image.data)
        ax.imshow(image.data, cmap='seismic', norm=remap, origin='lower')
    #end if
#end definition

#-----------------------------Spiral arms-----------------------------

# SpiralArms
#   A function to plot the spiral on top of the image.
#   Takes an object of the Spiral class and an object of the Image class (see above)
def SpiralArms(spiral, image):
    
    #Define the maximum angle that the arms will wind through
    theta_max = 3*np.pi
    
    #Convert the pitch angle to radians
    pitch_rad = spiral.pitch*np.pi/180
    b = np.tan(abs(pitch_rad))
    
    #Find the diameter and radius of the galaxy image
    gal_pixel_diameter = float(len(image.data))
    gal_pixel_radius = gal_pixel_diameter/2
    
    #Find the scale factor a by assuming a max radius slightly larger than the actual image radius
    a = (1.15*gal_pixel_radius)/np.exp(b*theta_max)
    
    #Convert the rotation angle to radians
    rotation_angle_rad = float(spiral.rotation)*np.pi/180
    
    #Create an array for theta values. For small pitch angles (<15 degrees), it may be a good
    #   idea to increase the radians to more than 6.
    theta = np.arange(spiral.theta_min,theta_max,0.05)
    
    #Initialize numpy arrays for the coordinates
    all_x = np.array([])
    all_y = np.array([])
    
    # Plot each spiral arm.
    for i in range(0,spiral.num_arms):
        #Add the rotation angle to the phase
        phase=i*2*np.pi/spiral.num_arms + rotation_angle_rad
        #Create arrays (implicit numpy arrays) for the x- & y-values using the theta array.
        if spiral.chirality == 'CCW':
            x = (-a*np.cos(theta+phase)*np.exp(b*theta) + gal_pixel_radius)
        elif spiral.chirality == 'CW':
            x = (a*np.cos(theta+phase)*np.exp(b*theta) + gal_pixel_radius)
        #end if
        
        y = (a*np.sin(theta+phase)*np.exp(b*theta) + gal_pixel_radius)
        all_x = np.append(all_x, x)
        all_y = np.append(all_y, y)
    #end loop

    return all_x, all_y

#end definition


#The flags here control which button starts out pressed
if DEFAULT_PITCH >= 0:
    chirality = 'CCW'
    chir_flag = 0
    pitch = DEFAULT_PITCH
elif DEFAULT_PITCH < 0:
    chirality = 'CW'    #If the pitch is negative, flip the chirality and drop the sign
    chir_flag = 1
    pitch = np.abs(DEFAULT_PITCH)
#end if

#Instantiate a Spiral object
overlay = Spiral(DEFAULT_ARMS, pitch, chirality, DEFAULT_ROTATION, DEFAULT_THETA_MIN, line_color)

#Plot the image and overlay, saving the lines object for future modification
ImagePlot(gal_image)
x,y = SpiralArms(overlay, gal_image)
overlay.lines, = ax.plot(x, y, '.', ms = 5.0)

if overlay.color == 'red':
    plt.setp(overlay.lines, color = 'red')
elif overlay.color == 'black':
    plt.setp(overlay.lines, color = 'black')
elif overlay.color == 'white':
    plt.setp(overlay.lines, color = 'white')
else:
    plt.setp(overlay.lines, color = 'blue')
#end if

#-----------------------------Basic GUI elements-----------------------------
#Define a consistent UI color
axcolor = 'white'

#Define a new "axis" for the slider
rot_ax = plt.axes([0.20, 0.16, 0.65, 0.03], axisbg=axcolor)

#Slider ranges from -180..180, initialized at the user's value (default 0.0).
rot_slider = Slider(rot_ax, 'Rotation angle (deg)', -180.0, 180.0, valinit=overlay.rotation)
rot_slider.label.set_size(10)
rot_slider.label.set_horizontalalignment('right')

pitch_ax = plt.axes([0.20, 0.12, 0.65,0.03], axisbg=axcolor)
pitch_slider = Slider(pitch_ax, 'Pitch angle (deg)', 0.001, 45.00, valinit=overlay.pitch)
pitch_slider.label.set_size(10)
pitch_slider.label.set_horizontalalignment('right')


theta_ax = plt.axes([0.20, 0.08, 0.65,0.03], axisbg=axcolor)
theta_slider = Slider(theta_ax, 'Minimum Theta (rad)', 0.0, 3*np.pi, valinit=overlay.theta_min)
theta_slider.label.set_size(10)
theta_slider.label.set_horizontalalignment('right')


# update
#   Define a set of actions to take when the slider is modified
#   Note that since the image and all the arms have to be redrawn, this function is fairly
#   slow. As such, it's best to click on the slider rather than click and drag.
#   Changing the rotation can also make the image jump, but the built-in cursor can fix that.
def Update(val):
    
    overlay.setRotation(rot_slider.val)
    overlay.setPitch(pitch_slider.val)
    overlay.setTheta(theta_slider.val)
    
    x,y = SpiralArms(overlay, gal_image)
    overlay.lines.set_xdata(x)
    overlay.lines.set_ydata(y)
    
    #Force the buffer to flip
    fig.canvas.draw_idle()
#end definition

rot_slider.on_changed(Update)
pitch_slider.on_changed(Update)
theta_slider.on_changed(Update)

#Define a new "axis" for the reset button
reset_ax = plt.axes([0.47, 0.03, 0.15, 0.04])
reset_button = Button(reset_ax, 'Reset sliders', color=axcolor, hovercolor='lightsteelblue')

# Reset
#   Define actions to take when  reset button is pressed
def Reset(event):
    #Reset sliders to initial values
    rot_slider.reset()
    pitch_slider.reset()
    theta_slider.reset()
#end definition

reset_button.on_clicked(Reset)

#Define a new axis for the save button
save_ax = plt.axes([0.80, 0.03, 0.1, 0.04])
save_button = Button(save_ax, 'Save', color=axcolor, hovercolor='lightsteelblue')

# Save
#   Define actions to take when save button is pressed
def Save(event):
    #Find the extent of the subplot that shows the image
    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(gal_image.name+'_'+str(overlay.pitch)+'degrees'+'.png', bbox_inches=extent)
#end definition

save_button.on_clicked(Save)

#Initialize a set of radio buttons for colorscale
log_ax = plt.axes([0.025, 0.60, 0.17, 0.15], axisbg=axcolor)
log_radio = RadioButtons(log_ax, ('linear','log grey', 'log seismic'), active=colorscale_flag)
fig.text(0.10,0.76,"Colorscale", fontsize='12',ha='center')

# ScaleChange
#   Allow the user to swap between linear/log scales on the fly
def ScaleChange(label):
    
    if label == 'linear':
        gal_image.setColorscale('n')
    elif label == 'log grey':
        gal_image.setColorscale('yg')
    elif label == 'log seismic':
        gal_image.setColorscale('ys')
    #end if
    ax.cla()

    ImagePlot(gal_image)
    x,y = SpiralArms(overlay, gal_image)
    overlay.lines, = ax.plot(x, y, '.', ms = 5.0)
    if overlay.color == 'red':
        plt.setp(overlay.lines, color = 'red')
    elif overlay.color == 'black':
        plt.setp(overlay.lines, color = 'black')
    elif overlay.color == 'white':
        plt.setp(overlay.lines, color = 'white')
    else:
        plt.setp(overlay.lines, color = 'blue')
    #end if
    fig.canvas.draw_idle()
#end definition

log_radio.on_clicked(ScaleChange)

#Initialize a set of radio buttons for chirality
chir_ax = plt.axes([0.8155, 0.30, 0.15, 0.15], axisbg=axcolor)
chir_radio = RadioButtons(chir_ax, ('CCW (-)','CW (+)'), active=chir_flag)
fig.text(0.8855,0.46,"Chirality", fontsize='12',ha='center')

# ChirChange
#   Allow the user to change the overlay's chirality on the fly
def ChirChange(label):
    
    if label == 'CCW (-)':
        overlay.setChirality('CCW')
    elif label == 'CW (+)':
        overlay.setChirality('CW')
    #end if
    
    x,y = SpiralArms(overlay, gal_image)
    overlay.lines.set_xdata(x)
    overlay.lines.set_ydata(y)

    fig.canvas.draw_idle()
#end definition

chir_radio.on_clicked(ChirChange)

#Initialize a set of radio buttons for the number of arms
arm_ax = plt.axes([0.8255, 0.50, 0.12, 0.25], axisbg=axcolor)
arm_radio = RadioButtons(arm_ax, ('1','2','3','4','5','6'), active=(DEFAULT_ARMS-1))
fig.text(0.8855,0.76,"Number of Arms", fontsize='12',ha='center')

# ArmChange
#   Allow the user to change the number of overlaid arms on the fly
def ArmChange(label):
    
    overlay.setNumArms(int(label))
    
    x,y = SpiralArms(overlay, gal_image)
    overlay.lines.set_xdata(x)
    overlay.lines.set_ydata(y)
    
    fig.canvas.draw_idle()
#end definition

arm_radio.on_clicked(ArmChange)

#Initialize a set of radio buttons for the line color
color_ax = plt.axes([0.032, 0.30, 0.12, 0.25], axisbg=axcolor)
color_radio = RadioButtons(color_ax, ('default','red','black','white'), active=0)
fig.text(0.10,0.56,"Line Color", fontsize='12',ha='center')

# ColorChange
#   Allow the user to change the color of overlaid arms on the fly
def ColorChange(label):

    overlay.setColor(label)
    
    if overlay.color == 'red':
        plt.setp(overlay.lines, color = 'red')
    elif overlay.color  == 'black':
        plt.setp(overlay.lines, color = 'black')
    elif overlay.color  == 'white':
        plt.setp(overlay.lines, color = 'white')
    else:
        plt.setp(overlay.lines, color = 'blue')
    #end if
    
    fig.canvas.draw_idle()
#end definition

color_radio.on_clicked(ColorChange)

plt.show()

#-----------------------------end of file-----------------------------
