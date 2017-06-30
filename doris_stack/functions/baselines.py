import os
import numpy as np
import warnings
from shutil import copyfile
from doris.doris_stack.main_code.resdata import ResData
import datetime
import subprocess


def baselines(dir_in,inputfile,start_date='2014-01-01',end_date='2018-01-01',doris=''):
    # This function calculates the baselines and plots a baseline plot.

    # Define doris path
    if not doris:
        doris = doris_path

    if not os.path.exists(dir_in):
        warnings.warn('The input directory does not exist!')
        return

    os.chdir(dir_in)
    process_folder = os.path.join(dir_in, 'baseline_process')
    if not os.path.exists(process_folder):
        os.mkdir(process_folder)
    os.chdir(process_folder)

    try:
        first = np.datetime64(start_date)
        last = np.datetime64(end_date)
    except:
        warnings.warn('Input dates could not be converted, use "yyyy-mm-dd"')
        return

    # Search for folders and take only the first burst.
    folders = next(os.walk(dir_in))[1]
    folders = sorted(folders)

    # Initialize... (Search for folders / resfiles / dates)
    n = 0
    res = []; date = []
    for fold in folders:
        # Select only the folders which has a name like yyyymmdd and are within
        if len(fold) == 8:
            # define date of folder
            date_prod = np.datetime64((fold[:4] + '-' + fold[4:6] + '-' + fold[6:]))

            if date_prod >= first and date_prod <= last:
                # Select the first swath
                date_fold = os.path.join(dir_in,fold)
                swath_fold = os.path.join(date_fold,next(os.walk(date_fold))[1][0])
                # Select the first burst
                prod_files = next(os.walk(swath_fold))[2]
                for file in prod_files:
                    if file.endswith('1.res'):
                        res.extend([os.path.join(swath_fold,file)])
                        date.extend([date_prod])
                        n = n + 1
                        break

    # Now create a set of baselines

    baselines = np.zeros([len(res),len(res)])
    resfiles = dict()

    # First create the ifgs.res files and store the data in a res data class.
    master = res[0]
    copyfile(master,os.path.join(process_folder,'master.res'))

    for resultfile, dat in zip(res, date):
        copyfile(resultfile,os.path.join(process_folder,'slave.res'))
        subprocess.call([doris + ' ' + inputfile], shell=True)

        dat = dat.astype(datetime.datetime).strftime('%Y-%m-%d')
        resfiles[dat] = ResData(type='interferogram',filename='ifgs.res')
        resfiles[dat].read()
        os.remove(os.path.join(process_folder,'ifgs.res'))

    # Then gather the baselines
    for dat, j in zip(date, range(len(date))):
        dat = dat.astype(datetime.datetime).strftime('%Y-%m-%d')
        baselines[j,0] = resfiles[dat].processes['coarse_orbits']['Bperp'][1]



    # Create figure of baselines.
    days = (date[0] - date).astype(float)
    plt.figure(111)
    plt.plot(baselines[:,0], days, marker='o')

    # Annotate
    for dat, x, y in zip(date, baselines[:,0], days):
        dat = dat.astype(datetime.datetime).strftime('%Y-%m-%d')
        plt.annotate(
            dat,
            xy = (x, y), xytext = (0, 0),
            textcoords = 'offset points', size = 8)

    plt.savefig('baseline_plot.pdf')

