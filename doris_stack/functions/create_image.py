# This function creates images based on complex or real input data.
# Input parameters are:
# - Input matrix
# - datatype
# - Image scaling (are pixels for example 1x2 km or 3x1 km?) scaling is azimuth / range.
# - Use of logscaling?
# - Plot amplitude / phase / both (not relevant for real values...)

# If you want to save your data as a geotiff or netcdf file use the read_write_data.py script. This will enable
# visualization in for example QGIS or ArcGIS
# If you want to do multilooking first apply the multilook.py script.

def create_image():
