#! /usr/bin/env python3
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from datetime import datetime
import os
import yaml
import sys

Deflate_Level=1

def print_message(message, error=False):
    import time
    if error:
        print('{}: ERROR {}'.format(time.asctime(),message))
        sys.exit(1)
    else:
        print('{}: {}'.format(time.asctime(),message))

def str_to_bool(s):
    ''' Convert string (s) to boolean
        ARGS:
         s: string
        RETURNS:
         boolean
    '''
    if s in ('True','true','T','t'):
        return True
    else:
        return False

# Read yaml control file (provided as a command line argument)
# If not present print usage and stop
if len(sys.argv) > 1:
    try:
        control = yaml.load(open(sys.argv[1]),Loader=yaml.BaseLoader)
    except Exception as e:
        print_message(e)
        print_message('loading {0}'.format(sys.argv[1]), error=True)
else:
    print_message('Usage: make_destripe.py <make_destripe.yml>')
    sys.exit(1)

# Get list of l2 files
input_files = control['Input files']

# Get number of standard deviation considered for RMS filter
nsigma = float(control['nsigma'])
# Get cloud fraction limits
mincfr = float(control['mincfr'])
maxcfr = float(control['maxcfr'])
# Get good main data quality flag values
mqfval = [int(val) for val in control['mqfval']]

# Logical to save diagnostic fields
yn_diag = str_to_bool(control['yn_diagnostic'])

# Attributes that will go into correction file
time_reference = []
scan_num = []
granule_num = []
processing_version = int(control['processing_version'])
time_coverage_start = []
time_coverage_end = []
time_coverage_start_since_epoch = []
time_coverage_end_since_epoch = []
begin_date = []
begin_time = []
end_date = []
end_time = []
production_date_time = datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
input_granules = []
access_description = control['access_description']
title = control['title']
cirrected_product = []

# Get number of xtrack positions to initialize masked arrays
# Get product name
try:
    print_message('getting number of xtrack positions from {}'.format(input_files[0]))
    with Dataset(input_files[0],'r') as src:
        nx = src.dimensions['xtrack'].size
        corrected_product = src.product_type
        time_reference = src.time_reference
except Exception as e:
    print_message(e)
    print_message('getting number of xtrack positions from {}'.format(input_files[0]), error=True)

lon = ma.array(np.zeros([1,nx],dtype=np.float32),mask=True)
lat = ma.array(np.zeros([1,nx],dtype=np.float32),mask=True)
mqf = ma.array(np.zeros([1,nx],dtype=np.int16),mask=True)
rms = ma.array(np.zeros([1,nx],dtype=np.float32),mask=True)
cfr = ma.array(np.zeros([1,nx],dtype=np.float32),mask=True)
scd = ma.array(np.zeros([1,nx],dtype=np.float32),mask=True)
vcd = ma.array(np.zeros([1,nx],dtype=np.float32),mask=True)
amf = ma.array(np.zeros([1,nx],dtype=np.float32),mask=True)

# Loop over L2 files
for fp in input_files:
    print_message('reading from {}'.format(fp))
    try:
        with Dataset(fp,'r') as src:
            lo = src['geolocation']['longitude'][:]
            la = src['geolocation']['latitude'][:]
            m = src['product']['main_data_quality_flag'][:]
            r = src['qa_statistics']['fit_rms_residual'][:]
            a = src['support_data']['amf'][:]
            c = src['support_data']['amf_cloud_fraction'][:]
            s = src['support_data']['fitted_slant_column'][:]
            p = src['support_data']['gas_profile'][:]
            units = src['support_data']['fitted_slant_column'].units
            v = ma.sum(p,axis=2); del p

            # Read granule attributes
            scan_num.append(src.scan_num)
            granule_num.append(src.granule_num)
            time_coverage_start.append(src.time_coverage_start)
            time_coverage_end.append(src.time_coverage_end)
            time_coverage_start_since_epoch.append(src.time_coverage_start_since_epoch)
            time_coverage_end_since_epoch.append(src.time_coverage_end_since_epoch)
            beg = src.time_coverage_start.split('T')
            end = src.time_coverage_end.split('T')
            begin_date.append(beg[0])
            begin_time.append(beg[1].strip('Z'))
            end_date.append(end[0])
            end_time.append(end[1].strip('Z'))
            input_granules.append(src.local_granule_id)
    except Exception as e:
        print_message(e)
        print_message('reading from {}'.format(fp),error=True)

    # Append to empty arrays
    lon = ma.concatenate([lon,lo],axis=0)
    lat = ma.concatenate([lat,la],axis=0)
    mqf = ma.concatenate([mqf,m],axis=0)
    rms = ma.concatenate([rms,r],axis=0)
    amf = ma.concatenate([amf,a],axis=0)
    cfr = ma.concatenate([cfr,c],axis=0)
    scd = ma.concatenate([scd,s],axis=0)
    vcd = ma.concatenate([vcd,v],axis=0)

    del m,r,c,s,v

# calculate limit for outliers of fitting rms
# Here, 3*k*MAD is considered, but it could change.
k = 1.4826
msig = k*ma.median(ma.abs(rms - ma.median(rms)))
rmslim = ma.median(rms) + nsigma*msig

# Get RMS mask
rms_mask = (rms > rmslim)
# Get cloud mask
cfr_mask = (cfr < mincfr) & (cfr > maxcfr)
# Get main data quality flag mask
mqf_mask = mqf.mask
for val in mqfval:
    tmp_mask = mqf != val
    mqf_mask = tmp_mask | mqf_mask

# Mask pixels not to be used
scd[cfr_mask | mqf_mask | rms_mask] = ma.masked
vcd[cfr_mask | mqf_mask | rms_mask] = ma.masked
amf[cfr_mask | mqf_mask | rms_mask] = ma.masked
cfr[cfr_mask | mqf_mask | rms_mask] = ma.masked
rms[cfr_mask | mqf_mask | rms_mask] = ma.masked

# median value of scds and model scds for de-striping and background corrections
# It's not stable, we need to update it.
# Calculate medians for each xtrack position
medval = ma.median(scd,axis=0)
bgrcor = ma.median(vcd.astype(np.float)*amf.astype(np.float),axis=0)
bgrcor = bgrcor.astype(np.float32)
vcd_med = ma.median(vcd,axis=0)
amf_med = ma.median(amf,axis=0)
cfr_med = ma.median(cfr,axis=0)
rms_med = ma.median(rms,axis=0)

# python netcdf output concatenates string attributes without delimiters,
# so we explicitly concatenate here to add delimiters.
time_coverage_start = ','.join(time_coverage_start)
time_coverage_end = ','.join(time_coverage_end)
begin_date = ','.join(begin_date)
begin_time = ','.join(begin_time)
end_date = ','.join(end_date)
end_time = ','.join(end_time)
input_granules = ','.join(input_granules)

local_granule_id = control['output_filename']

# Now Create a file to save the results
print_message('writing correction results to {}'.format(local_granule_id))
try:
    with Dataset(local_granule_id,'w',clobber=True) as dst:
        dst.time_reference = time_reference
        dst.scan_num = scan_num
        dst.granule_num = granule_num
        dst.processing_version = processing_version
        dst.time_coverage_start = time_coverage_start
        dst.time_coverage_end = time_coverage_end
        dst.time_coverage_start_since_epoch = time_coverage_start_since_epoch
        dst.time_coverage_end_since_epoch = time_coverage_end_since_epoch
        dst.begin_date = begin_date
        dst.begin_time = begin_time
        dst.end_date = end_date
        dst.end_time = end_time
        dst.production_date_time = production_date_time
        dst.local_granule_id = local_granule_id
        dst.input_granules = input_granules
        #dst.shortname = shortname
        dst.access_description = access_description
        dst.title = title
        dst.corrected_product = corrected_product
        dst.number_standard_deviations = nsigma
        dst.minimum_cloud_fraction = mincfr
        dst.maximum_cloud_fraction = maxcfr
        dst.used_quality_flags = mqfval
        # Create dimensions
        dst_nx = dst.createDimension('xtrack',nx)
        dst_nm = dst.createDimension('mirror_step',amf.shape[0])
        # Create variables
        if yn_diag:
            dst_lon = dst.createVariable('longitude',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_lat = dst.createVariable('latitude',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_vcd = dst.createVariable('vcd',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_scd = dst.createVariable('scd',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_amf = dst.createVariable('amf',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_cfr = dst.createVariable('cfr',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_rms = dst.createVariable('rms',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_msk = dst.createVariable('mask',np.float32,('mirror_step','xtrack'),fill_value=False,zlib=True,complevel=Deflate_Level)

            dst_vcd_med = dst.createVariable('median_vcd',np.float32,('xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_amf_med = dst.createVariable('median_amf',np.float32,('xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_cfr_med = dst.createVariable('median_cfr',np.float32,('xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_rms_med = dst.createVariable('median_rms',np.float32,('xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)

            dst_lon[:] = lon
            dst_lat[:] = lat
            dst_vcd[:] = vcd
            dst_scd[:] = scd
            dst_amf[:] = amf
            dst_cfr[:] = cfr
            dst_rms[:] = rms
            msk = amf * 0.0
            msk[rms.mask] = 1.0
            dst_msk[:] = msk

            dst_vcd_med[:] = vcd_med
            dst_amf_med[:] = amf_med
            dst_cfr_med[:] = cfr_med
            dst_rms_med[:] = rms_med

        dst_bgr = dst.createVariable('background_correction',np.float32,('xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
        dst_bgr.title = '{} background correction'.format(corrected_product)
        dst_bgr.units = units
        dst_bgr[:] = bgrcor

        dst_des = dst.createVariable('destriping_correction',np.float32,('xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
        dst_des.title = '{} destriping correction'.format(corrected_product)
        dst_des.units = units
        dst_des[:] = medval

except Exception as e:
    print_message(e)
    print_message('writing correction results to {}'.format(local_granule_id),error=True)
