#!/bin/python3

import mpl2nc
from pathlib import Path
import numpy as np

import datetime
import argparse

from functools import reduce
from netCDF4 import Dataset

from mpl2nc import NC_TYPE, FILL_VALUE, __version__


def join_ds(first, second):


    times_first = first['time'].shape[0]
    times_second = second['time'].shape[0]

    second['file_index_start'] = times_first

    print(times_first, times_second)


    out = {}
    for k in first.keys():
        #print(k)
        if isinstance(first[k], np.ndarray) and isinstance(second[k], np.ndarray):
            #print('an array', first[k].shape, second[k].shape)
            if len(second[k].shape) == 0:
                arr = np.hstack((first[k], second[k]))
            elif second[k].shape[0] == times_second:
                #print('an array with time')
                arr = np.concatenate((first[k], second[k]), axis=0)
                #print('arr.shape', arr.shape)
            elif len(first[k].shape) == 1 and len(second[k].shape) == 1:
                #print(first[k], second[k])
                arr = np.stack((first[k], second[k]), axis=1).T
                # averaged profiles per file
                #print('arr.shape', arr.shape)
            elif len(second[k].shape) == 1:
                arr = np.concatenate((first[k], second[k][None,:]))
            else:
                print('still something missed')
        elif isinstance(first[k], np.ndarray):
            arr = np.append(first[k], second[k])
        else:
            arr = np.hstack((first[k], second[k]))
        #print(arr.shape)
        out[k] = arr
    return out



def write(d, filename, nc_header):
    f = Dataset(filename, 'w')
    f.createDimension('profile', None)
    f.createDimension('range', None)
    f.createDimension('ap_range', None)
    f.createDimension('ol_range', None)
    f.createDimension('no_files', None)
    if 'dt_coeff' in d:
        f.createDimension('dt_coeff_degree', None)
    if 'dt_count' in d:
        f.createDimension('dt_count', None)
    for k, v in d.items():
        #print(k)
        h = nc_header[k]
        #print(h)
        var = f.createVariable(k, NC_TYPE[h['dtype']], h['dims'],
            fill_value=FILL_VALUE[h['dtype']])
        #print(var.shape)
        var[::] = v
        if h['units'] is not None: var.units = h['units']
        if h['long_name'] is not None: var.long_name = h['long_name']
        if h['comment'] is not None: var.comment = h['comment']
    f.created = datetime.datetime.utcnow().strftime('%Y-%m-%dT:%H:%M:%SZ')
    f.software = 'mpl2nc (https://github.com/peterkuma/mpl2nc)'
    f.version = __version__
    f.close()
            

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--date', help='date in the format YYYYMMDD')
    args = parser.parse_args()
    date = datetime.datetime.strptime(args.date, '%Y%m%d')

    #mpl2nc -a MMPL5081_Afterpulse_202203100900.bin -d deadtime_coeffs.csv -o MMPL5081_Overlap_202203111400.bin mpl_neumayer3/2023/09/11 converted/2023/

    calibration_coeff = 1
    calibration_coeff = 1.04e-5
    #outfile = Path(f'converted_daily/{date:%Y%m%d}_mpl.nc')
    outfile = Path(f'mpl_converted_daily/neumayer/{date:%Y}/{date:%Y%m%d}_mpl.nc')
    outfile.parent.mkdir(parents=True, exist_ok=True)
    print('outfile', outfile)
    
    file_afterpulse = 'MMPL5081_Afterpulse_202203100900.bin'
    file_deadtime = 'deadtime_coeffs.csv'
    file_overlap = 'MMPL5081_Overlap_202203111400.bin'
    #file_overlap = 'MMPL5081_Overlap_202312302200.bin'
    file_overlap = 'MMPL5081_Overlap_202312302200_pollyNR_correct_first_bin.bin'
    
    p = sorted(Path(f"mpl_neumayer3/{date:%Y}/").glob(f"nm{date:%y%m%d}*"))

    #p = sorted(Path(f"mpl_neumayer3/{date:%Y}/{date:%m}/{date:%d}/").glob(f"{date:%Y%m%d}*.mpl"))
    
    
    afterpulse, overlap, deadtime = {}, {}, {}
    print(p)
    if file_afterpulse is not None:
        with open(file_afterpulse, 'rb') as f:
            afterpulse = mpl2nc.read_afterpulse(f)
    
    if file_overlap is not None:
        with open(file_overlap, 'rb') as f:
            overlap = mpl2nc.read_overlap(f)
    
    if file_deadtime is not None:
        if file_deadtime.endswith('.csv'):
            dead_time = mpl2nc.read_dt_csv(file_deadtime)
        else:
            with open(file_deadtime, 'rb') as f:
                dead_time = mpl2nc.read_dt(f)
    
    d = {}
    d.update(afterpulse)
    d.update(overlap)
    d.update(dead_time)
    
    
    mpl_all = []
    for filename in p[:]:
    
        mpl = mpl2nc.read_mpl(filename)
        mpl.update(d)
        mpl2nc.process_nrb(mpl)
        mpl_all.append(mpl)

    no_files = len(mpl_all)
    # make the joints traceable
    mpl_all[0]['file_index_start'] = 0
    
    mpl = reduce(join_ds, mpl_all)

    rgs = np.full_like(mpl['channel_1'], np.nan, np.float64)
    for i in range(rgs.shape[0]):
        rgs[i,:] = 0.5*mpl['bin_time'][i]*mpl2nc.C*(np.arange(rgs.shape[1]) + 0.5)
    assert (rgs == rgs[0]).all()
    mpl['range'] = rgs[0,:]
    mpl['ap_range'] = mpl['ap_range'][0,:]
    if 'ol_range' in mpl:
        mpl['ol_range'] = mpl['ol_range'][0,:]
    mpl['dt_count'] = mpl['dt_count'][0,:]

    mpl['backscatter'] = (mpl['nrb_copol'][:] + 2.*mpl['nrb_crosspol'][:])*calibration_coeff
    mpl['voldepolratio'] = mpl['nrb_crosspol'][:] / (mpl['nrb_copol'][:] + mpl['nrb_crosspol'][:])

    nc_header = mpl2nc.NC_HEADER
    for k in nc_header.keys():
        #print(k, nc_header[k]['dims'])
        if nc_header[k]['dims'] == []:
            nc_header[k]['dims'] = ['no_files']
        if nc_header[k]['dims'] == ['ap_range']:
            nc_header[k]['dims'] = ['no_files', 'ap_range']
        if nc_header[k]['dims'] == ['ol_range']:
            nc_header[k]['dims'] = ['no_files', 'ol_range']
        if nc_header[k]['dims'] == ['dt_count']:
            nc_header[k]['dims'] = ['no_files', 'dt_count']
        if nc_header[k]['dims'] == ['dt_coeff_degree']:
            nc_header[k]['dims'] = ['no_files', 'dt_coeff_degree']
        #print(k, nc_header[k]['dims'])

    nc_header['ap_range']['dims'] = ['ap_range']
    nc_header['ol_range']['dims'] = ['ol_range']
    nc_header['dt_count']['dims'] = ['dt_count']
    nc_header['file_index_start'] = {
            'dtype': 'uint32', 'long_name': 'profile_index_file_start',
            'comment': 'profile index where the original .mpl file started',
            'units': 'index', 'dims': ['no_files']}
    nc_header['range'] = {
            'dtype': 'float32', 'long_name': 'range',
            'comment': 'range claculated from bin_time',
            'units': 'm', 'dims': ['range']}
    nc_header['backscatter'] = {
            'dtype': 'float32', 'long_name': 'backscatter_coefficient',
            'comment': 'calculated from NRB with (copol + 2*crosspol)*calibration_coeff',
            'units': '', 'dims': ['profile', 'range']}
    nc_header['voldepolratio'] = {
            'dtype': 'float32', 'long_name': 'volume_depolarization_ratio',
            'comment': 'calculated from NRB with crosspol/(copol+crosspol)',
            'units': '', 'dims': ['profile', 'range']}

    write(mpl, outfile, nc_header)
