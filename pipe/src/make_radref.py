#! /usr/bin/env python3

import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from datetime import datetime
import os
import yaml
import sys
# import tracemalloc

def print_error_message(message):
    import time
    print('{}: ERROR: {}'.format(time.asctime(),message))
    sys.exit(1)

def print_message(message):
    import time
    print('{}: {}'.format(time.asctime(),message))

# Read yaml control file (provided as a command line argument)
# If not present print usage and stop
if len(sys.argv) > 1:
    try:
        control = yaml.load(open(sys.argv[1]),Loader=yaml.BaseLoader)
    except Exception as e:
        print_error_message('loading {}'.format(sys.argv[1]))
        print_error_message(e)
else:
    print_message('Usage: >make_radref.py <make_radref.yml>')
    sys.exit(1)

# Check how is the memory usage
# tracemalloc.start()

# Get list of L1 RAD files to be used in calculation
L1RAD_files = control['L1RAD Input']
# Get list of L2 CLOUD files
L2CLD_files = control['L2CLD Input']
# Get minimum and maximum cloud fraction to be used
mincfr = float(control['mincfr'])
maxcfr = float(control['maxcfr'])
xtstep = int(control['xtstep'])
nsigma = float(control['nsigma'])
print_message('mincfr: {}; maxcfr: {}; nsigma: {}'.format(mincfr,maxcfr,nsigma))

# Attributes that will go into radiance reference file
time_reference = []
scan_num = []
granule_num = []
processing_version = int(control['processing_version'])
scan_type = []
time_coverage_start = []
time_coverage_end = []
time_coverage_start_since_epoch = []
time_coverage_end_since_epoch = []
production_date_time = datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
radiance_granules = []
cloud_granules = []
begin_date = []
begin_time = []
end_date = []
end_time = []
access_description = control['access_description']
title = control['title']

# Loop over files and extract xtrack and spectral_channel dimensions
# If xtrack and spectral_channel is not equal in all scans we need
# to think a more elaborate way of calculating the radiance reference
#  nxs = []; nws = []
print_message('Collecting input granules metadata information')
for fprad,fpcld in zip(L1RAD_files,L2CLD_files):
    try:
        with Dataset(fprad,'r') as src:
            # gf = src['granule_flag'][:]
            # If granule flag is telemetry only (4) skip this file
            # Other flag meanings: 0 (nominal), 1(first scan), 2 (last scan)
            # if (gf in [4]):
            #     print_message('Warning!!! skipping file {0}'.format(fp))
            #     print_message('           Granule flag {0}'.format(gf))
            #     continue
            # Read granule dimensions
            # nxs.append(src['band_290_490_nm'].dimensions['xtrack'].size)
            # nws.append(src['band_290_490_nm'].dimensions['spectral_channel'].size)
            # Read granule attributes
            time_reference = src.time_reference
            scan_num.append(src.scan_num)
            granule_num.append(src.granule_num)
            scan_type.append(src.scan_type)
            time_coverage_start.append(src.time_coverage_start)
            time_coverage_end.append(src.time_coverage_end)
            time_coverage_start_since_epoch.append(src.time_coverage_start_since_epoch)
            time_coverage_end_since_epoch.append(src.time_coverage_end_since_epoch)
            radiance_granules.append(src.local_granule_id)
            begin_date.append(src.begin_date)
            begin_time.append(src.begin_time)
            end_date.append(src.end_date)
            end_time.append(src.end_time)
            # Read granule xtrack and spectral channel dimensions
            # They should never change
            nx = src['band_290_490_nm'].dimensions['xtrack'].size
            nw = src['band_290_490_nm'].dimensions['spectral_channel'].size
            # Read spectrum and wavelength units
            spectrum_units = src['band_290_490_nm']['radiance'].units
            wavelength_units = src['band_290_490_nm']['nominal_wavelength'].units
        with Dataset(fpcld,'r') as src:
            cloud_granules.append(src.local_granule_id)
    except Exception as e:
        print_message(e)
        print_error_message('getting dimensions and metadata from L1 RAD file {}'.format(fprad))
# if (np.size(np.unique(np.array(nxs))) != 1):
#     print_error_message('Radiance reference calculation can not continue since xtrack dimension changes with file')
# else:
#     nx = nxs[0]
#     del nxs
# if (np.size(np.unique(np.array(nws))) != 1):
#     print_error_message('Radiance reference calculation can not continue since spectral_channel dimension changes with file')
# else:
#     nw = nws[0]
#     del nws

# python netcdf output concatenates string attributes without delimiters,
# so we explicitly concatenate here to add delimiters.
time_coverage_start = ','.join(time_coverage_start)
time_coverage_end = ','.join(time_coverage_end)
begin_date = ','.join(begin_date)
begin_time = ','.join(begin_time)
end_date = ','.join(end_date)
end_time = ','.join(end_time)
radiance_granules = ','.join(radiance_granules)
cloud_granules = ','.join(cloud_granules)

local_granule_id = control['output_file']

# Create arrays to hold radiance reference, wavelength and associated flags
radref = ma.zeros([nx,nw],dtype=np.float32)
wvlref = ma.zeros([nx,nw],dtype=np.float32)
mdqref = ma.zeros([nx,nw],dtype=np.short)
numref = ma.zeros([nx],dtype=np.int16)

# To keep memory usage under control loop over xtracks
for ix in range(0,nx,xtstep):
    # size, peak = tracemalloc.get_traced_memory()
    # print_message('size: {0}; peak: {1}'.format(size/1024.0/1024.0/1024.0,peak/1024.0/1024.0/1024.0))
    # Create temporal radiance array in which to concatenate rads
    # considering the latest iteration
    fx = ix+xtstep
    if fx > nx:
        fx = nx
        radtmp = ma.array(np.zeros([1,nx-ix,nw],dtype=np.float32),mask=True)
    else:
        radtmp = ma.array(np.zeros([1,xtstep,nw],dtype=np.float32),mask=True)
    print_message('Processing xtracks {} to {}'.format(ix,fx))
    for fprad, fpcld in zip(L1RAD_files,L2CLD_files):
        # Read L1 radiance files
        try:
            with Dataset(fprad,'r') as radsrc:
                gf = radsrc['granule_flag'][:]
                # If granule flag is telemetry only (4) skip this file
                # Other flag meanings: 0 (nominal), 1(first scan), 2 (last scan)
                if (gf in [4]):
                    print_message('WARNING: skipping file {}'.format(fp))
                    print_message('          Granule flag {}'.format(gf))
                    continue
                # Read granule dimensions
                nt = radsrc.dimensions['mirror_step']
                # Read variables of interest
                nominal_wvl = radsrc['band_290_490_nm']['nominal_wavelength'][:,:]
                rad = radsrc['band_290_490_nm']['radiance'][:,ix:fx,:]
                pqf = radsrc['band_290_490_nm']['pixel_quality_flag'][:,ix:fx,:]
                # Read L2 cloud files
                try:
                    with Dataset(fpcld,'r') as cldsrc:
                        # Read cloud fraction
                        cfr = cldsrc['product']['cloud_fraction'][:,ix:fx]
                except Exception as e:
                    print_message(e)
                    print_error_message('reading L2 cloud file {}'.format(fpcld))
                # For each xtrack position keep only pixels that fill the cloud filter
                cld_mask = (cfr < mincfr) & (cfr > maxcfr)
                # Only use spectral pixels with pqf equal to 0
                pqf_mask = (pqf != 0)
                # Apply this mask to the radiances
                rad[cld_mask,:] = ma.masked
                rad[pqf_mask] = ma.masked
                radtmp = ma.concatenate([rad,radtmp],axis=0)
                # size, peak = tracemalloc.get_traced_memory()
                # print_message('size: {0}; peak: {1}'.format(size/1024.0/1024.0/1024.0,peak/1024.0/1024.0/1024.0))
                del rad, pqf, cld_mask, pqf_mask
        except Exception as e:
            print_message(e)
            print_error_message('reading L1 RAD file {}'.format(fprad))

    # Calculate the median of radiance means (for the whole band) for each xtrack position
    radmeans = ma.mean(radtmp,axis=2)
    radmed   = ma.median(radmeans,axis=0)
    mad      = ma.median( ma.abs(radmeans - radmed), axis = 0)
    # Calculate the range of standard deviations of radiances
    k = 1.4826
    radsig = k * ma.median( ma.abs( radmeans - radmed), axis=0) * nsigma
    # Now loop over cross track positions and calculate the final radiance
    # Find out for each cross track position keep (radmeans > radmed - radsig) & (radmeans < radmed + radsig)
    for i,j in enumerate(range(ix,fx)):
        rad_mask = (radmeans[:,i] > radmed[i] - radsig[i]) & (radmeans[:,i] < radmed[i] + radsig[i])
        radref[j,:] = ma.mean(radtmp[:,i,:][rad_mask],axis=0)
        wvlref[j,:] = nominal_wvl[j,:]
        mdqref[j,:][radref[j,:].mask] = 1
        numref[j]   = ma.sum(rad_mask)

    del radmeans, radmed, radsig, mad, radtmp

# Save calculated radiance reference and associated variables to radref file.
# First form a unique file name using date,scan,start granule and end granule.
with Dataset(local_granule_id,'w',clobber=True) as dst:
    # First create global attributes
    dst.time_reference = time_reference
    dst.scan_num = scan_num
    dst.granule_num = granule_num
    dst.processing_version = processing_version
    dst.scan_type = scan_type
    dst.time_coverage_start = time_coverage_start
    dst.time_coverage_end = time_coverage_end
    dst.time_coverage_start_since_epoch = time_coverage_start_since_epoch
    dst.time_coverage_end_since_epoch = time_coverage_end_since_epoch
    dst.production_date_time = production_date_time
    dst.local_granule_id = local_granule_id
    dst.radiance_granules = radiance_granules
    dst.cloud_granules = cloud_granules
    #dst.shortname = shortname
    dst.begin_date = begin_date
    dst.begin_time = begin_time
    dst.end_date = end_date
    dst.end_time = end_time
    dst.access_description = access_description
    dst.title = title
    dst.number_standard_deviations = nsigma
    dst.minimum_cloud_fraction = mincfr
    dst.maximum_cloud_fraction = maxcfr
    # Create dimensions
    dst_nw = dst.createDimension('spectral_channel',nw)
    dst_nx = dst.createDimension('xtrack',nx)
    # Create variables
    dst_spc = dst.createVariable('spectrum',np.float64,('xtrack','spectral_channel'),fill_value=-1.0e+30,zlib=True,complevel=4)
    dst_spc.title = 'radiance spectrum'
    dst_spc.units = spectrum_units
    dst_spc[:] = radref

    dst_wav = dst.createVariable('wavelength',np.float64,('xtrack','spectral_channel'),fill_value=-1.0e+30,zlib=True,complevel=4)
    dst_wav.title = 'wavelengths'
    dst_wav.units = wavelength_units
    dst_wav[:] = wvlref

    dst_qfg = dst.createVariable('quality_flag',np.short,('xtrack','spectral_channel'),fill_value=-999,zlib=True,complevel=4)
    dst_qfg.title = 'spectrum quality flag'
    dst_qfg.flag_meanings='0: good, non-zero: bad'
    dst_qfg[:] = mdqref

    dst_num = dst.createVariable('number_coadd',np.int16,('xtrack'),fill_value=-999,zlib=True,complevel=4)
    dst_num.title='number of co-added spectra'
    dst_num.description='number of spectra considered in the mean to derive the radiance reference'
    dst_num[:] = numref
