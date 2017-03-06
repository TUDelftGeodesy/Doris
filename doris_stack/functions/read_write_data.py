import numpy as np
import os


def freadbk(path_file, line_start=1, pixel_start=1, nofLines=None, nofPixels=None, dt='float32', lines=0, pixels=0):
    # First use memmap to get a memory map of the full file, than extract the requested part.

    if dt == 'cpxint16':
        dtype = np.dtype([('re', np.int16), ('im', np.int16)])
        file_dat = np.memmap(path_file, dtype=dtype, mode='r', shape=(lines, pixels)).view(np.int16).astype(np.float32).view(np.complex64)
        data = file_dat[line_start - 1:line_start + nofLines - 1, pixel_start - 1:pixel_start + nofPixels - 1].astype(
            'complex64', subok=False)
    elif dt == 'cpxshort':
        dtype = np.dtype([('re', np.float16), ('im', np.float16)])
        file_dat = np.memmap(path_file, dtype=dtype, mode='r', shape=(lines, pixels)).view(np.float16).astype(np.float32).view(np.complex64)
        data = file_dat[line_start - 1:line_start + nofLines - 1, pixel_start - 1:pixel_start + nofPixels - 1].astype(
            'complex64', subok=False)
    else:
        dt = np.dtype(dt)
        file_dat = np.memmap(path_file, dtype=dt, mode='r', shape=(lines, pixels))
        data = file_dat[line_start - 1:line_start + nofLines - 1, pixel_start - 1:pixel_start + nofPixels - 1].astype(
            dt, subok=False)

    return data, file_dat


def fwritebk(path_file, data, dt):
    # First define dtype and write to file using memmap.

    if dt == 'cpxint16':
        dtype = np.dtype([('re', np.int16), ('im', np.int16)])
        data = np.memmap(path_file, dtype=dtype, mode='w', shape=data.shape)
        data[:, :] = data.view(np.float32).astype(np.int16).view(dtype)
    elif dt == 'cpxshort':
        dtype = np.dtype([('re', np.float16), ('im', np.float16)])
        data = np.memmap(path_file, dtype=dtype, mode='w', shape=data.shape)
        data[:, :] = data.view(np.float32).astype(np.float16).view(dtype)
    else:
        data = np.memmap(path_file, dtype=dt, mode='w', shape=data.shape)
        data[:, :] = data

    return data, file_dat

def read_tiff(path_file, line_start=1, pixel_start=1, nofLines=None, nofPixels=None, dt='float32'):
    print('under construction')


def write_tiff(path_file, line_start=1, pixel_start=1, nofLines=None, nofPixels=None, dt='float32'):
    print('under construction')


def read_nc(path_file, line_start=1, pixel_start=1, nofLines=None, nofPixels=None, dt='float32'):
    print('under construction')


def write_nc(path_file, line_start=1, pixel_start=1, nofLines=None, nofPixels=None, dt='float32'):
    print('under construction')



def python_gdal_convert():
    print('under construction')