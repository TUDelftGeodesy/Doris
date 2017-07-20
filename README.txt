===========
Doris v5 Beta
===========

The new Doris version, Doris 5, is developed to process a stack of Sentinel-1 images, as well as all the already familiar
functionality of Doris 4.

This is a beta version. Therefore, you still may experience some problems. Please report them to us. But even better,
try to fix them! We are very happy to discuss with you how you can contribute to this project!

The new Doris version consists of 2 parts:
-       The doris_core directory cantaining the Doris core code, which is similar to the original Doris code and is
        written in C. This code is mainly used to create individual interferograms based on different steps.
-       The doris_stack directory containing scripts written in Python. These scripts automate the processing of a
        single master stack of Sentinel-1 images. The scripts manage the processing of the bursts of a sentinel 1 image,
        contain algorithms specific to processing sentinel 1 images and support parallelisation of the processing of the
        bursts. The functionality of these scripts can be further extended to support more sensors and modes.

Note that the python code is developed in python 2.7, so be sure you are not using python 3.

In addition, you will find a stack preparation script, to automatically download the burst you need for your Area of
Interest which you defined by a shape file, automatically download the SRTM DEM associated with this area, and setup
your processing structure.


Installation
===========

See the INSTALL file in the install directory. This file descripes the installation of the C libraries, python libraries
and some utility software.


Creating Sentinel-1 datastacks
=============================


Create a folder structure
-----------------------------

After installing the software you can create your first doris datastack. To do so you have to prepare the following:
- Create folders to download radar data and orbit files. In a further stage these files can be downloaded automatically,
    but it is also possible to do it yourself manually.
- Create a folder where you can store intermediate DEM results. Data will be downloaded automatically, so you only have
    to create the folder itself. Note that these automatic downloads are based on SRTM data and are therefore limited to
    60 degrees south and north of the equator.
- Create a .shp file with your area of interest. You can use different software packages, but ArcGIS and QGIS (free) are
    the most convenient for this purpose. Or you could download from one of the websites that offer free shapefiles for
    administrative boundaries (for example: http://www.diva-gis.org/Data)
- Finally, create the folder where you want to process your datastack. Be aware that to process your data you will need
    at least 100 GB of free space on your disc.


Register for Sentinel and SRTM downloads
----------------------------------------

Additionally, you will need an account for downloading Sentinel-1 and SRTM data. You can use the following links to
create an account. (How to link them to your code is described in the INSTALL file)
- To register for Sentinel-1 data download use: https://scihub.copernicus.eu/dhus/#/self-registration
- To register for SRTM download use: https://urs.earthdata.nasa.gov/users/new/


Run the stack preparation script
----------------------------------------

Move to the prepare_stack directory:
cd prepare_stack
Run the python script:
python prepare_datastack_main.py

This code will ask you to define the different folders you created before. The script will ask you whether you want
to run your code in parallel. Generally, this is recommended as it speeds up your processing speed. Note that either the
number of cores and your RAM can be limiting (one process will use about 4GB of RAM). Because it is not possible to mix
different orbits in one datastack it will also ask you which orbit you want to use and whether it is ascending or
descending. Please check this beforehand on the ESA website (https://scihub.copernicus.eu)
Finally, the code will ask you the start date, end date and master date:
- start date    > What is the first image (in time) you want to process?
- end date      > What is the last image (in time) you want to process? (Tip: This date can be in the far future if you
                    just want to download all images till now)
- master data   > This image will be used as the master of your stack. Other images will be resampled
                    to the geometry of this master image.
After finishing this script, the new datastack is automatically created together with a DEM of the area. This can take
a while in case the download speeds are low or your area is large.


Editing the data stack settings (generally not needed)
----------------------------------------------------

You can enter the folder to find, the newly created DEM, your .shp file, configuration files (inputfiles) and the stack.
Further there is the doris_input.xml file where all configuration settings for your datastack are stored.
This file is created in the folder where you will process your datastack. So, if you want to change this configuration 
afterwards, you can make some changes there.


Processing
=========================================

In the main folder of your datastack you will have three bash files:
create_dem.sh           > To create a DEM for your area. This is already done if you used the automatic DEM generation
download_sentinel.sh    > This will run the a download of sentinel images for the specified track over your area of
                            interest. Only dates between your start and end date are considered. This script will also
                            download the needed precise or restituted orbit files.
You can call this scripts using bash <script_name>

After downloading your DEM, radar data and orbit files you can start your processing by the following command:
bash doris_stack.sh

or, if your server uses qsub (for parallel processing)

qsub doris_stack.sh

If you want to extend your datastack later on, you can run the scripts again for the same datastack. It will check which
files are new and only process them. This software is therefore perfectly fit for continues monitoring.
Be sure that you do not change your master image in between, as this will break your code.



Enjoy,

TUDELFT RADAR GROUP 2017
doris_users@tudelft.nl