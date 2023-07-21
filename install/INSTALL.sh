cd ../doris_core
sudo apt-get install tcsh
./configure
make
make install

# --------------------------------------------
# - COMPILATION OF THE DORIS CORE ------------
# --------------------------------------------
cd ../doris_core
# 8. Read the README file
./configure             #(creates "Makefile")  # requires tcsh shell to run, to install type "" at shell prompt sudo apt-get install
                            #tcsh on Ubuntu platform.
                            #( +answer the questions about libraries, etc.)
make                  #  (compiles the software)
make install          #  (installs doris and bin scripts)


#--------------------------------------------
# - COMPILATION OF DORIS UTILITIES -----------
# --------------------------------------------
cd ../sartools
make
# 14. Review/edit the Makefile if this does not work
#     (for example if you do not want to use GNU gcc/g++ as compiler)
make install           # (installs in /usr/local/bin unless you edit the Makefile)


cd ../envisat_tools       # on 64-bit system requires libc-dev-i386 library ex: "sudo apt-get install libc-dev-i386"
make
# 18. Review/edit the Makefile if this does not work
#     (for example if you do not want to use gcc as compiler)
make install

# --------------------------------------------
# - INSTALLATION OF Doris - python part
# --------------------------------------------

# To use the Doris Python scripts you will need to install the following Python packages:
# -       numpy, scipy (for calculations)
# -       matplotlib (visualization)
# -       requests (data download)
# -       gdal, gdal-dev, shapely, fiona, pyproj, fastkml, osr (GIS)


# After installation of the python packages you have to make a link to the doris_core, cpxfiddle and snaphu executables.
# You can do so by executing the install script, but you will have to create user accounts first to be able to download
# Sentinel-1 and SRTM DEMs.

# 25. create an account for the downloading of SRTM DEMs at
# https://urs.earthdata.nasa.gov/users/new
# 26. create an account for downloading Sentinel data at
# https://urs.earthdata.nasa.gov/users/new/

# move to the install directory
cd ../install
python init_cfg.py
# and fill in the different paths and user accounts
