'''
    --------------------------------
    2DFFTUtils
    --------------------------------
    
    A module containing helpful utilities for use with the 2DFFT method of
    spiral galaxy logarithmic pitch angle measurement.
    
    Dependencies:
        -NumPy: http://www.numpy.org/
        -IRAF: http://iraf.noao.edu/
        -PyRaf: http://www.stsci.edu/institute/software_hardware/pyraf
        -AstroPy: http://www.astropy.org/
        -2DFFT: https://mseigar.wixsite.com/mysite/2dfft
    These URLs are accurate as of April 2017.
    
    All of the above are inluded bundled with Python (and other useful tools like ds9) in the Space
    Telescope Science Institute's UREKA package: http://ssb.stsci.edu/ureka/
    However, UREKA has recently been deprecated, so users looking for continued support will need to 
    install the relevant libraries a la carte.
    
    These functions were written with batch automation in mind, and  for the most part they operate on 
    lists of files. However, they can easily be used on a single file by wrapping the filename in brackets 
    as in: utilFunction(["filename.fits"]).
    
    This module was developed on Mac OS 10.11 using Python 2.7.
    Individual utilities have been tested only on OS 10 and Ubuntu in Python 2.7.
    
'''


import numpy as np
import glob
#from pylab import *
#import matplotlib.pyplot as plt
import re
from pyraf import iraf
import sys
from os.path import splitext
from astropy.io import fits

__all__ = ["getIter","scripterList","getCenter","getSimRadius","getRadius","crop","autoCrop","convertText","avgPitch","getStableRegions"]

#---------------------------------Convenience Functions---------------------------------

def getIter(regex):
    ''' 
        getIter: Uses glob and regular expressions in concert to do more rigorous pattern
                 matching than glob alone to get specified files.
        Inputs:  a string <regex> defining a regular expression. See https://docs.python.org/2/library/re.html for
                 information on regular expressions in Python.
        Outputs: an iterable containing the filenames matching the pattern specified in <regex>.
    '''
    out = []
    #all_fits = glob.glob("*.fit") + glob.glob("*.fits")
    all_fits = glob.glob("*")
    for name in all_fits:
        match = re.search(regex, name)
        if match:
            out.append(name)
        #end if
    #end loop
    return out
#end definition

def scripterList(textfiles,keywords,radii,outfile='scripter-list.txt'):
    '''
        scripterList:   Creates a CSV file to feed to the 2DFFT scripter. See the readme for the 2DFFT package.
        Inputs:         <textfiles> a list of the filenames of the text versions of the images to be processed
                        by 2DFFT.
                        <keywords>, a list which specifies the names for the 2DFFT output. These are the prefixes
                        for the _m files.
                        <radii>, an array formatted like the output for getRadius or getSimRadius specifying the 
                        outer radii for the objects to be processed.
                        <outfile>, the filename for the output, can optionally be specified. The default is 
                        "scripter-list.txt".
                        All of these inputs must be the same length and must relate to the same objects, in the same
                        order, or the output file will not work as expected.
        Outputs:        a text file with the name specified by <outfile>, formatted like
                                                          [blank line]
                                            image_textfile_1,keyword_1,outer_radius_1
                                            image_textfile_2,keyword_2,outer_radius_2
                                            image_textfile_3,keyword_3,outer_radius_3
                                                              ...
                                            last_image_textfile,last_keyword,last_outer_radius
                                                          [blank line]
        Notes:          The outer radius in the output is 1 pixel less than the one specified in <radii>
                        Based on code by J.E. Berlanga Medina.
    '''
    if len(textfiles)!=len(keywords):
        print("Input lists must be the same length.")
        sys.exit(1)
    #end if
    if len(textfiles)!=len(radii):
        print("Input lists must be the same length.")
        sys.exit(1)
    #end if
    printfile = open(outfile,'w')

    printfile.write('\n')
    for i in range(0,len(textfiles)):
        printfile.write(textfiles[i])
        printfile.write(',')
        printfile.write(keywords[i])
        printfile.write(',')
        outer_radius = str(int(radii[i,1])-1)
        printfile.write(outer_radius)
        printfile.write('\n')
    #end loop
    printfile.write('\n')
#end definition

#-----------------------------------Image Preparation-----------------------------------

def getCenter(imgs):
    ''' 
        getCenter: Uses the Iraf imcntr task to find the centers of a list of objects
        Inputs:    a list of image filenames <imgs> (see getIter()) to find the centers of.
        Outputs:   a numpy array with shape (len(imgs), 3), of the form
                                   ...         ...         ...
                                <img_name>  <center_x>  <center_y>
                                   ...         ...         ...
                   where all three fields are strings.
    '''
    out = np.full((len(imgs),3),"-1",dtype=object)
    shape_re = re.compile(r"\[(?P<rows>\d+),(?P<cols>\d+)\]")
    center_re = re.compile(r"x:\s+(?P<center_x>\d+\.\d*)\s+y:\s+(?P<center_y>\d+\.\d*)")
    for i in range(0,len(imgs)):
        header_string = iraf.imheader(imgs[i],Stdout=1)[0]
        shape = shape_re.search(header_string)
        rows,cols = -1,-1
        if shape:
            rows = int(shape.group("rows"))
            cols = int(shape.group("cols"))
        else:
            sys.exit(1)
        #end if
        center_string = iraf.imcntr(imgs[i],rows/2,cols/2,cboxsize=31, Stdout=1)[0]
        center = center_re.search(center_string)
        center_x,center_y = "-1","-1"
        if center:
            center_x = str(int(float(center.group("center_x"))))
            center_y = str(int(float(center.group("center_y"))))
        else:
            sys.exit(1)
        #end if
        out[i,0] = imgs[i]
        out[i,1] = center_x
        out[i,2] = center_y
    #end loop
    return out
#end definition

def getSimRadius(imgs, centers):
    '''
       getSimRadius: Calculates an outer radius for n-body simulation galaxies (i.e. fairly simple galaxies on a
                     background of zeroes).
       Inputs:       <imgs>, a list of images to retrive a radius for, and <centers> an array of object centers
                     formatted like the output from getCenter (help(getCenter)). These two imputs must be in the same order and
                     of the same length.
       Outputs:      a two-column array of the form
                                       ...                 ...
                                    <img_name>        <outer_radius>
                                       ...                 ...
                    where both fields are strings.
        Notes:      This method scales poorly: as O(mn) where m and n are the image width and height. For most simulated
                    images it remains faster than choosing an outer radius by hand, but use caution and consider using
                    using getRadius for large images.
                    Simulated images created by some software may have scales embedded at the edges of the image in the data
                    layer of the FITS file. These will need to be removed in order to properly calculate a radius. This is easy to
                    automate with AstroPy.
    '''
    out = np.full((len(imgs),2),"-1",dtype=object)
    for k in range(len(imgs)):
        if imgs[k] == centers[k,0]:
            out[k,0] = imgs[k]
            center_x = int(centers[k,1])
            center_y = int(centers[k,2])
            img_data = fits.getdata(imgs[k], ignore_missing_end = True) #Grab the pixel values from the image
            rows,cols = img_data.shape
            radius = 1
            for i in range(rows):
                for j in range(cols):
                    if img_data[i][j] != 0: #If the pixel at (j,i) isn't black ...
                        distance = int(np.ceil(np.sqrt((j-center_x)**2+(i-center_y)**2))) #find the distance from the center to the pixel at (j,i) ...
                        if distance > radius: #and if it's greater than the current record ...
                            radius = distance #set the distance as the new record.
                        #end if
                    #end if
                #end loop
            #end loop
            out[k,1] = str(radius)
        #end if
    #end loop
    return out
#end definition

def getRadius(imgs, centers, a=0.03, f='fro'):
    '''
        getRadius: A tunable method for calculating the outer radius of galaxies. This function applies a specified function to incrementally
                   smaller square portions of the image centered on the object under study and compares the previous and current values. When 
                   the percent difference between these two values exceeds a specified value, iteration stops and the function returns half 
                   the side length of the current portion as the outer radius of the object.
        Inputs:    <imgs>, a list of image filenames to find a radius for and <centers> an array of object centers formatted like
                   the output from getCenter (help (getCenter)). These two imputs must be in the same order and of the same length.
                   A parameter <a>, required to be a float between 0 and 1, can also be specified to determine when iteration stops.
                   The default is a = 0.03.
                   A parameter <f> specifying a matrix norm: https://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.norm.html
                   The default is the Frobenius norm, which is related most closely to the total brightness.
        Outputs:   a two-column array of the form
                                       ...                 ...
                                    <img_name>        <outer_radius>
                                       ...                 ...
                where both fields are strings.
    '''
    out = np.full((len(imgs),2),"-1",dtype=object)
    for j in range(len(imgs)):
        if imgs[j] == centers[j,0]:
            out[j,0] = imgs[j]
            center_x = int(centers[j,1])
            center_y = int(centers[j,2])
            img_data = fits.getdata(imgs[j], ignore_missing_end = True)
            rows,cols = img_data.shape
            done = False
            i = 5
            last = np.linalg.norm(img_data, ord=f)
            start = np.minimum(np.minimum(center_x,center_y),np.minimum(cols-center_x,rows-center_y))
            radius = start
            while(not done and i < start):
                
                # Using the frobenius norm avoids overflow and produces a human-readable value associated with the image
                val = np.linalg.norm(img_data[center_y-radius:center_y+radius:,center_x-radius:center_x+radius:], ord=f)
                #print("norm at r=" + str(cols/2-i) + " : " + str(val) + '\n')
                
                if (np.abs(val-last) >= int(a*last)):
                    done = True
                #end if
                
                radius = radius - i
                i += 1
                last = val
            #end loop
            out[j,1] = str(radius)
        #end if
    return out
#end definition
                

def crop(img, center_x, center_y, radius, append="_crop"):
    '''
        crop:   Crops a single FITS image using the IRAF imcopy task.
        Inputs: <img>, a filename specifying the image to be cropped, <center_x> and <center_y>, the 
                coordinates of the object center, and <radius>, the radius to crop to. Optionally, 
                a tag <append> to be appended to the filename can be specified. The default appendant is "_crop".
        Output: cropped image, written to disk
    '''
    base = splitext(img)[0]
    extension = splitext(img)[1]

    upper_x = center_x + radius
    lower_x = center_x - radius
    upper_y = center_y + radius
    lower_y = center_y - radius
    try:
        # These indices look backwards, but that's how the IRAF convention works.
        iraf.imcopy(img +"["+str(lower_x)+":"+str(upper_x)+","+str(lower_y)+":"+str(upper_y)+"]",base + append + extension,Stdout="/dev/null")
    except iraf.IrafError as e:
        print "IRAF Error: " + e
    #end try
#end definition

def autoCrop(imgs, centers, radii, append="_crop"):
    ''' 
        autoCrop: Crops a list of FITS images using the IRAF imcopy task.
        Inputs:   a list of FITS images <imgs>, an array <centers> specifying the object centers in the format
                  output by getCenter (help(getCenter), and an array <radii> specifying the radii to crop to. WARNING:in
                  order to get the expected results, <imgs>, <centers>, and <radii> must be THE SAME LENGTH and IN THE
                  SAME ORDER. Optionally, a tag <append> to be appended to the filename can be specified. The default
                  appendant is "_crop".
        Outputs:  cropped images, written to disk.
    '''
    for i in range(len(imgs)):
        crop(imgs[i], int(centers[i,1]), int(centers[i,2]), int(radii[i,1]), append)
    #end loop
#end definition

def convertText(imgs):
    ''' 
        convertText: Converts a FITS image to a text file (with .txt extension) using the IRAF
                     wtextimage task.
        Inputs:      a list of strings <imgs> containging the filenames of FITS format images.
        Outputs:     a file (written to disk) with the same basename as the input and the .txt 
                     extension, containing the data from the original image.
    '''
    for filename in imgs:
        base = splitext(filename)[0]
        try:
            iraf.wtextimage(filename,base+".txt",header='no',pixels='yes')
        except iraf.IrafError as e:
            print("IRAF Error: " + e)
    #end try
#end definition

#----------------------------------2DFFT Analysis----------------------------------

def avgPitch(snapshot_file,start_pixel,end_pixel):
    '''
        avgPitch: Gives average pitch angle, standard deviation error, and 2DFFT error for
                  a single mode file. Uses a range of pixels as input.
        Inputs:   <snapshot_file>, a string specifying a 2DFFT output file, i.e. one of the
                  <filename>_mX, for X=0..6, where X is the mode.
                  <start_pixel> and <end_pixel> are integers specifying the range to average 
                  the pitch angle over.
        Outputs:  Prints the snapshot name, pixel range, average pitch, standard deviation,
                  2DFFT error, and total error to the terminal. All measurements are in degrees.
        Notes:    Written by J.E. Berlanga Medina and lightly edited by E. Monson to add extra
                  error calculation and for consistency with the rest of 2DFFTUtils.
                  This is mainly a helper function for slopeChange, but it can be used on its own.
    '''
    
    # Get list of all the pitch angles in the mode file.
    allpitch = np.loadtxt(snapshot_file, usecols=(4,), unpack=True)
    
    # Sum up all the pitch angles in the list from start_pixel to end pixel.
    # **Note: The way indexing works, we must add 1 to the second index to
    # 	include end_pixel in the calculation.
    pitch_sum = sum(allpitch[start_pixel:end_pixel+1])
    
    #print("The sum of pitch angles in the above pixel range is: "+str(pitch_sum))
    pitch_avg = pitch_sum/float(end_pixel - start_pixel + 1)
    print("\n")
    print("Snapshot data file: "+snapshot_file)
    # Stable region range given in counting numbers that start at 1, not 0.
    print("Radial range in pixels: "+str(start_pixel+1)+"-"+str(end_pixel+1))
    print("Average pitch: "+str(pitch_avg))
    
    # Compute standard deviation; ddof means that standard deviation is computed
    #	as sqrt(sum[(pitch-mean)^2]/(N - ddof)), where ddof gives degreees of
    # 	freedom instead of sqrt(sum[(pitch-mean)^2]/N).
    pitch_std_dev = np.std(allpitch[start_pixel:end_pixel+1], dtype=np.float64, ddof=1)
    print("Standard deviation: "+str(pitch_std_dev))
    
    m = snapshot_file[-1:]
    x = np.absolute(pitch_avg)
    
    # 2DFFT error is calculated using the cubic polynomials presented in Figure 7
    #   of Davis et al. 2012.
    if m == '1':
        fft_error = -0.00004*x**3 + 0.0058*x**2 - 0.0137*x + 0.0234
    elif m == '2':
        fft_error = -0.00002*x**3 + 0.0029*x**2 - 0.0084*x + 0.0222
    elif m == '3':
        fft_error = -0.00001*x**3 + 0.0020*x**2 - 0.0064*x + 0.0214
    elif m == '4':
        fft_error = -0.00001*x**3 + 0.0015*x**2 - 0.0054*x + 0.0207
    elif m == '5':
        fft_error = -0.000009*x**3 + 0.0012*x**2 - 0.0046*x + 0.0200
    elif m == '6':
        fft_error = -0.000007*x**3 + 0.0010*x**2 - 0.0041*x + 0.0191
    elif m == '0':
        fft_error = None
    #end if

    print("2DFFT error: " + str(fft_error))

    # We calculate the total error as the sum in quadrature of the standard
    #   deviation and 2DFFT error.
    print("Total error: ") + str(np.sqrt(pitch_std_dev**2 + fft_error**2))
    print("\n")

#end definition

def getStableRegions(file_list):
    '''
        getStableRegions: uses a method developed by J.E. Berlanga Medina to detect and print
                          information on ranges in 2DFFT output where the pitch angle of the
                          object is relatively constant. Then, calls avgPitch to print out the
                          average pitch angle for the region and associated errors calculations.
        Inputs:           <file_list> a list of 2DFFT output files to examine. These are the _mX files,
                          where X=0..6 is the mode.
        Ouputs:           Prints to the terminal the:
                          -name of the input file
                          -outer radius of the object
                          -upper radial limit (0.95*radius)
                          -lower radial limit (0.20*radius)
                          -minimum stable region length (0.05*radius)
                          followed by, for each stable region found, the output from avgPitch.
        Notes:            Written by J.E. Berlanga Medina and lightly edited by E. Monson for consistency
                          with the rest of 2DFFTUtils.
    '''
    for input_file in file_list:
    
        # Load the pitch angle values into the array y_pitch.
        y_pitch = np.loadtxt(input_file, usecols=(4,), unpack=True)
        
        # Make the upper limit of the pixel range used roughly 95% of the galaxy's radius.
        # Make the lower limit of the pixel range used roughly 20% of the galaxy's radius.
        # Make the minimum length of the stable region roughly 5% of the galaxy's radius.
        upper_lim_pixel = int(np.ceil(0.95*len(y_pitch)))
        lower_lim_pixel = int(np.floor(0.20*len(y_pitch)))
        min_length = int(np.ceil(0.05*len(y_pitch)))
        
        print("------" + input_file + "------")
        print("Radius of the galaxy: "+str(len(y_pitch)))
        print("Upper radial limit used: "+str(upper_lim_pixel))
        print("Lower radial limit used: "+str(lower_lim_pixel))
        print("Minimum stable region length: "+str(min_length))
        print("\n")
        
        # Set the list that holds all the slopes to empty.
        all_slopes_list = []
        # Calculate slope for each line.
        # 	The slope of line 0 connects point 0 & point 1, etc.
        for i in range(lower_lim_pixel,upper_lim_pixel-1):
            all_slopes_list.append((y_pitch[i+1] - y_pitch[i]))
        #end loop
        
        
        # Import the list of slopes into a 2D array which has the index of the slopes
        # 	in the left column (pixel), and the slopes in the right column.
        all_slopes = np.array([[i+lower_lim_pixel,all_slopes_list[i]] for i in range(0,len(all_slopes_list))])
        #np.savetxt("all_slopes--"+input_file+".txt", all_slopes, fmt='%.3f', delimiter='	')

        # Now, we need to find out which slopes are greater than a threshhold value,
        # 	and cut those out.  Here, the threshhold value is 2.00.
        to_cut_slope_list = []
        
        for i in range(0,len(all_slopes)):
            if np.absolute(all_slopes[i][1])>float(2):
                to_cut_slope_list.append(i)
            #end if
        #end loop

        # Save all the rows we still want (those that contain the stable regions)
        # 	to a new array by deleting the rows with slopes that are too large.
        stable_regions=np.delete(all_slopes,to_cut_slope_list,0)
        #np.savetxt("stable_regions--"+input_file+".txt", stable_regions, fmt='%.3f', delimiter='	')
        # Take the absolute values of the slopes.
        stable_absolute=np.absolute(stable_regions)
        #np.savetxt("stable_absolute--"+input_file+".txt", stable_absolute, fmt='%.3f', delimiter='	')
        
        # Now, look at the array containing the absolute values of the slopes wherein
        # 	lies the stable region.
        # If two consecutive indices (pixels) in the array are not bordering pairs,
        # 	then we know that there is a break here & we must extract one or more
        # 	stable subrange of indices (pixels).
        
        # Make a list of pixels that enclose breaks.
        breaks_list = []
        
        for i in range(0,len(stable_absolute)-1):
            index_diff = stable_absolute[i+1][0]-stable_absolute[i][0]
            if index_diff!=float(1):
                breaks_list.append(stable_absolute[i][0])
                breaks_list.append(stable_absolute[i+1][0])
            #end if
        #end loop

        #np.savetxt("breaks_list--"+input_file+".txt", breaks_list, fmt='%.3f', delimiter='	')

        # Find the stable subregions.  If they are longer than some threshold value
        # 	of pixels (this can be a set number, or a percentage of the galaxy's radius
        # 	in pixels), then we call avgPitch for this subregion.

        number_stable_regions=0
        
        for j in range(0,len(breaks_list)):
            #print("j = "+str(j))
            if j==0:
                #print(str(j)+" should equal exactly 0.")
                min_radius = stable_regions[0][0]
                #print("min_radius = "+str(min_radius))
                max_radius = breaks_list[j]
                #print("max_radius = "+str(max_radius))
            #end if
            if j==int(len(breaks_list)-1):
                #print(str(j)+" should equal exactly 1 less than the length of the breaks_list, or: "+str(len(breaks_list) - 1))
                min_radius = breaks_list[j]
                max_radius = stable_regions[len(stable_regions)-1][0]
                #print("min_radius = "+str(min_radius)+" and max_radius = "+str(max_radius))
            #end if
            if j>0 and j%2==0 and j!=len(breaks_list)-1:
                #print(str(j)+" should be an even number, and not equal to len(breaks_list)-1, or: "+str(len(breaks_list)-1))
                min_radius = breaks_list[j-1]
                max_radius = breaks_list[j]
            #end if
            pixel_difference = max_radius - min_radius
            if pixel_difference>min_length:
                number_stable_regions = number_stable_regions + 1
                print("Stable region #"+str(number_stable_regions))
                avgPitch(input_file, min_radius, max_radius)
                max_radius = 0
                min_radius = 0
            #end if
        #end loop

        if number_stable_regions==0:
            print("No stable regions longer than 5 percent of the outer radius were found.")
            print("\n")
            #end if
    #end loop
#end definition

#----------------------------------end of file----------------------------------
