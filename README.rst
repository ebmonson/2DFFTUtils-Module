*****************
2DFFTUtils Module
*****************

The 2DFFTUtils Python module (``_2DFFTUtils.py``) implements most of the 2DFFT utilities in one place, to make
preparing images for 2DFFT and dealing with 2DFFT data interactively or in scripts even easier.

The overlay application is implemented in a separate script file, ``spiral_overlay.py``.

Scripts for plotting pitch angle as a function of radius from 2DFFT output are not included in ``_2DFFTUtils.py``. Scripts by Jazmin Berlanga Medina are available in the ``plots/`` section of ``2dfft_utils``.

Scripts for deprojecting galaxy images or converting simulation data from postscript format to FITS are also not included in this module.

More detailed documentation for the functions included in the 2DFFTUtils module is available using Python's builtin ``help`` command.

Included Functions
################

Image Preparation:
  - ``getCenter()``: find the centers of the objects in a list of galaxy images.
  - ``getSimRadius()``: given an array of centers, find an exact outer radius for a list of relatively simple objects on a background of zeroes (e.g., n-body simulations).
  - ``getRadius()``: given an array of centers, find an approximate outer radius for a list of objects on noisier backgrounds (e.g. observational data) using a method based on matrix norms.
  - ``crop()``: given an outer radius and a center, crop a galaxy image into a square centered on the galaxy.
  - ``autoCrop()``: given an array of outer radii and an array of centers, crop a list of galaxy images into squares centered on the objects using the ``crop()`` function.
  - ``convertText()``: given a list of images, convert them into text files (older versions of 2DFFT require input files to be in text form).
2DFFT Analysis:
  - ``averagePitch()``: given a 2DFFT output file and a range of pixels, calculate the average pitch angle over that range and associated errors.
  - ``getStableRegions()``: using a method developed by Jazmin Berlanga Medina, find regions in 2DFFT output where the change in pitch is relatively small. Then, call ``averagePitch()`` on the regions and present them to the user.
Convenience Functions:
  - ``getIter()``: given a regular expression pattern, return a list of files matching that pattern. Useful for grouping together simulations, observation runs, or 2DFFT output files.
  - ``scripterList()``: given a list of text files, object names, and outer radii, generate a CSV file used to generate the 2DFFT script.
