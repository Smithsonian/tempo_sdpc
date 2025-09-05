#! /usr/bin/env python3
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from datetime import datetime, timezone
from scipy.ndimage import generic_filter
from scipy.signal import savgol_filter
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
# Get angle limits
maxsza = float(control['maxsza'])
maxvza = float(control['maxvza'])
# Limits to filter for polluted pixels in a priori from model
# Valid values for type are 'none', 'troposphere' or 'total'
vcd_apriori_type = control['vcd_apriori_type']
maxvcd_apriori = float(control['maxvcd_apriori'])
if vcd_apriori_type == 'troposphere':
    filter_vcd_trop_ap = True
else:
    filter_vcd_trop_ap = False
# Limits to filter for VCD value from observation
# Valid values for type are 'none' or 'total'
vcd_obs_type = control['vcd_obs_type']
minvcd_obs = float(control['minvcd_obs'])
maxvcd_obs = float(control['maxvcd_obs'])
# Get good main data quality flag values
mqfval = [int(val) for val in control['mqfval']]
# Minimum number of pixels in xtrack required to compute correction
minpixels = int(control['minpixels'])
# Number of xtrack pixels and polynomial order to use for smoothing filter
smooth_nxtrack = int(control['smooth_nxtrack'])
smooth_polyorder = int(control['smooth_polyorder'])
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
production_date_time = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
input_granules = []
title = control['title']
corrected_product = []

max_mirror_step = 0
min_mirror_step = 0

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

lon = ma.array(np.zeros([1,nx],dtype=np.float64),mask=True)
lat = ma.array(np.zeros([1,nx],dtype=np.float64),mask=True)
mqf = ma.array(np.zeros([1,nx],dtype=np.int16),mask=True)
rms = ma.array(np.zeros([1,nx],dtype=np.float64),mask=True)
cfr = ma.array(np.zeros([1,nx],dtype=np.float64),mask=True)
sza = ma.array(np.zeros([1,nx],dtype=np.float64),mask=True)
vza = ma.array(np.zeros([1,nx],dtype=np.float64),mask=True)
scd = ma.array(np.zeros([1,nx],dtype=np.float64),mask=True)
amf = ma.array(np.zeros([1,nx],dtype=np.float64),mask=True)
vcd = ma.array(np.zeros([1,nx],dtype=np.float64),mask=True)
vcdtot_ap = ma.array(np.zeros([1,nx],dtype=np.float64),mask=True)
if filter_vcd_trop_ap:
    vcdtrp_ap = ma.array(np.zeros([1,nx],dtype=np.float64),mask=True)

# Loop over L2 files
for fp in input_files:
    print_message('reading from {}'.format(fp))
    try:
        with Dataset(fp,'r') as src:
            lo = src['geolocation']['longitude'][:]
            la = src['geolocation']['latitude'][:]
            m = src['product']['main_data_quality_flag'][:]
            r = src['qa_statistics']['fit_rms_residual'][:]
            grp_prod = src['product']
            grp_supp = src['support_data']
            if 'amf_total' in grp_supp.variables:
                a = src['support_data']['amf_total'][:]
            elif 'amf' in grp_supp.variables:
                a = src['support_data']['amf'][:]
            c = src['support_data']['amf_cloud_fraction'][:]
            sz = src['geolocation']['solar_zenith_angle'][:]
            vz = src['geolocation']['viewing_zenith_angle'][:]
            # Check to see if fitted_slant_column_uncorrected exists (this
            # means that the file has already been destriped). If so, use this
            # for the SCD.
            grp = src['support_data']
            if 'fitted_slant_column_uncorrected' not in grp.variables:
                s = src['support_data']['fitted_slant_column'][:]
            else:
                s = src['support_data']['fitted_slant_column_uncorrected'][:]
            p = src['support_data']['gas_profile'][:]
            sp = src['support_data']['surface_pressure'][:]
            if 'vertical_column_total' in grp_supp.variables:
                vt = src['support_data']['vertical_column_total'][:]
            elif 'vertical_column' in grp_prod.variables:
                vt = src['product']['vertical_column'][:]
            units = src['support_data']['fitted_slant_column'].units
            vtot = ma.sum(p,axis=2)  # vcd used here is a priori vcd

            if filter_vcd_trop_ap:
            # Calculate a priori VCD trop which is needed to filter polluted pixels
                tp = src['support_data']['tropopause_pressure'][:]
                eta_a = src['support_data']['surface_pressure'].Eta_A
                eta_b = src['support_data']['surface_pressure'].Eta_B
                pz = eta_a + eta_b * sp[:,:,np.newaxis]
                pz_diff = pz - tp[:,:,np.newaxis]
                tp_idx = np.where(pz_diff > 0, pz_diff, np.inf).argmin(axis=2) - 1
                vtrp = ma.array(np.zeros(tp_idx.shape))
                for iy, ix in np.ndindex(tp_idx.shape):
                    vtrp[iy,ix] = ma.sum(p[iy,ix,:tp_idx[iy,ix]])
                del pz, tp_idx, pz_diff

            del p

            # find the range of mirror step indices
            mirror_step = src.variables['mirror_step'][:]
            max_mirror_step = max(max_mirror_step, mirror_step.max())
            min_mirror_step = min(min_mirror_step, mirror_step.min())

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
            input_granules.append(os.path.basename(fp))
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
    sza = ma.concatenate([sza,sz],axis=0)
    vza = ma.concatenate([vza,vz],axis=0)
    scd = ma.concatenate([scd,s],axis=0)
    vcd = ma.concatenate([vcd,vt],axis=0)
    vcdtot_ap = ma.concatenate([vcdtot_ap,vtot],axis=0)
    if filter_vcd_trop_ap:
        vcdtrp_ap = ma.concatenate([vcdtrp_ap,vtrp],axis=0)
        del vtrp

    del m,r,c,s,vt,vtot

##############################
# Calculate destriping factors
##############################
# calculate limit for outliers of fitting rms
# Here, 3*k*MAD is considered, but it could change.

k = 1.4826
msig = k*ma.median(ma.abs(rms - ma.median(rms)))
rmslim = ma.median(rms) + nsigma*msig

# Get RMS mask
rms_mask = (rms > rmslim)
# Get cloud mask
cfr_mask = (cfr < mincfr) | (cfr > maxcfr)
# Geometry mask
ang_mask = (sza > maxsza) | (vza > maxvza)
# Filter out polluted pixels from model a priori
if filter_vcd_trop_ap:
    vcdap_mask = (vcdtrp_ap > maxvcd_apriori)
elif vcd_apriori_type == 'total':
    vcdap_mask = (vcdtot_ap > maxvcd_apriori)
elif vcd_apriori_type == 'none':
    vcdap_mask = False
else:
    print_message('Invalid vcd_apriori_type threshold type: do not filter.')
    vcdap_mask = False
# Filter out pixels that have VCD outside defined range
if vcd_obs_type == 'total':
    vcd_mask = (vcd > maxvcd_obs) | (vcd < minvcd_obs)
elif vcd_obs_type == 'none':
    vcd_mask = False
else:
    print_message('Invalid option for vcd_obs_type outlier type: do not filter.')
    vcdap_mask = False
# Get main data quality flag mask
use_idx = False
for val in mqfval:
    use_idx = (mqf == val) |  use_idx
mqf_mask = ~use_idx | mqf.mask

# Mask pixels not to be used
dstr_mask = cfr_mask | ang_mask | vcdap_mask | vcd_mask | mqf_mask | rms_mask
scd[dstr_mask] = ma.masked
vcdtot_ap[dstr_mask] = ma.masked
vcd[dstr_mask] = ma.masked
amf[dstr_mask] = ma.masked
cfr[dstr_mask] = ma.masked
sza[dstr_mask] = ma.masked
vza[dstr_mask] = ma.masked
rms[dstr_mask] = ma.masked
mqf[dstr_mask] = ma.masked
if filter_vcd_trop_ap:
    vcdtrp_ap[dstr_mask] = ma.masked

# Find number of valid samples in each cross track row
number_samples = scd.count(axis=0)

# Calculate a priori SCD from model
scd_ap = vcdtot_ap * amf

# Difference between fitted and a priori SCD
scd_obs_ap_diff = scd - scd_ap
scd_obs_ap_diff_xtrack = np.nanmean(scd_obs_ap_diff,axis=0)
np.ma.masked_invalid(scd_obs_ap_diff_xtrack)
# # Smooth difference to extract high frequency stripe component
# # This works but is not good at tracking some anomalies due to fires, so
# # implementing a Savitzky-Golay filter which can handle these better.
# scd_diff_smooth = (generic_filter(scd_obs_ap_diff_xtrack,np.nanmean, \
#                                   size=smooth_nxtrack).astype(np.float32))

# The Savitzky-Golay filter cannot handle NaN or masked values, so
# interpolate any missing values in the middle, and do not use the
# most extreme North/South pixels where not enough valid data (5 positions
# in the North and 7 positions in the South are always bad, but also check
# to see if there are others bordering these as these can mess up filtering).
first_north_index = np.argmax(number_samples >= minpixels)
try:
    last_south_index = np.where(number_samples >= minpixels)[0][-1]

    scd_obs_ap_diff_xtrack_tmp = scd_obs_ap_diff_xtrack[first_north_index:last_south_index+1]
    nans = np.isnan(scd_obs_ap_diff_xtrack_tmp.data)
    x = np.arange(len(scd_obs_ap_diff_xtrack_tmp))
    if nans.any():
        scd_obs_ap_diff_xtrack_tmp[nans] = np.interp(x[nans], x[~nans], scd_obs_ap_diff_xtrack_tmp[~nans])

    scd_diff_smooth_tmp = savgol_filter(scd_obs_ap_diff_xtrack_tmp, smooth_nxtrack,
                                        smooth_polyorder, mode='mirror')
    scd_diff_smooth = np.zeros(nx)
    scd_diff_smooth[first_north_index:last_south_index+1] = scd_diff_smooth_tmp
    # Subtract the smoothed difference from the original to get the final values
    # of the high frequency stripes
    stripe_val_tmp = scd_obs_ap_diff_xtrack - scd_diff_smooth
    # Set values where there not enough samples to zero so these positions are not destriped
    stripe_val = np.where(number_samples < minpixels, 0, stripe_val_tmp)
except Exception as e:
    print_message(e)
    print_message('Not enough valid input values: Set destriping correction to zero.')
    scd_diff_smooth = np.zeros(nx)
    stripe_val = np.zeros(nx)

# Calculate medians for each xtrack position
vcdtot_ap_med = ma.median(vcdtot_ap,axis=0)
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

# Now create a file to save the results
print_message('writing destriping correction to {}'.format(local_granule_id))
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
        dst.title = title
        dst.corrected_product = corrected_product
        dst.number_standard_deviations = nsigma
        dst.minimum_cloud_fraction = mincfr
        dst.maximum_cloud_fraction = maxcfr
        dst.maximum_solar_zenith_angle = maxsza
        dst.maximum_viewing_zenith_angle = maxvza
        dst.used_quality_flags = mqfval
        dst.num_mirror_pos = max_mirror_step - min_mirror_step + 1
        # Create dimensions
        dst_nx = dst.createDimension('xtrack',nx)
        dst_nm = dst.createDimension('mirror_step',amf.shape[0])
        # Create variables
        if yn_diag:
            dst_lon = dst.createVariable('longitude',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_lat = dst.createVariable('latitude',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_vcdtot_ap = dst.createVariable('vcd_apriori',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_scd_ap = dst.createVariable('scd_apriori',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_scd = dst.createVariable('scd',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_amf = dst.createVariable('amf',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_vcd = dst.createVariable('vcd',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_cfr = dst.createVariable('cfr',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_sza = dst.createVariable('sza',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_vza = dst.createVariable('vza',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_rms = dst.createVariable('rms',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_mqf = dst.createVariable('mdqf',np.int16,('mirror_step','xtrack'),fill_value=-1,zlib=True,complevel=Deflate_Level)
            dst_msk = dst.createVariable('mask',np.float32,('mirror_step','xtrack'),fill_value=False,zlib=True,complevel=Deflate_Level)
            if filter_vcd_trop_ap:
                dst_vcdtrp_ap = dst.createVariable('vcd_trop_apriori',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)

            dst_vcdtot_ap_med = dst.createVariable('median_vcd_apriori',np.float32,('xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_amf_med = dst.createVariable('median_amf',np.float32,('xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_cfr_med = dst.createVariable('median_cfr',np.float32,('xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_rms_med = dst.createVariable('median_rms',np.float32,('xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)

            dst_lon[:] = lon
            dst_lat[:] = lat
            dst_vcdtot_ap[:] = vcdtot_ap
            dst_scd_ap[:] = scd_ap
            dst_scd[:] = scd
            dst_amf[:] = amf
            dst_vcd[:] = vcd
            dst_cfr[:] = cfr
            dst_sza[:] = sza
            dst_vza[:] = vza
            dst_rms[:] = rms
            dst_mqf[:] = mqf
            msk = amf * 0.0
            msk.mask = False
            msk[dstr_mask] = 1.0
            dst_msk[:] = msk

            if filter_vcd_trop_ap:
                dst_vcdtrp_ap[:] = vcdtrp_ap

            dst_vcdtot_ap_med[:] = vcdtot_ap_med
            dst_amf_med[:] = amf_med
            dst_cfr_med[:] = cfr_med
            dst_rms_med[:] = rms_med

            # Save difference between observed and a priori SCD (all pixels)
            dst_scddiff = dst.createVariable('scd_difference',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_scddiff.title = '{} SCD difference (observed - a priori)'.format(corrected_product)
            dst_scddiff.units = units
            dst_scddiff[:] = scd_obs_ap_diff

        # Save difference between observed and a priori SCD (averaged over xtrack)
        dst_scddiff_xtrack = dst.createVariable('scd_difference_xtrack',np.float32,('xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
        dst_scddiff_xtrack.title = '{} SCD difference by xtrack (observed - a priori)'.format(corrected_product)
        dst_scddiff_xtrack.units = units
        dst_scddiff_xtrack[:] = scd_obs_ap_diff_xtrack

        # Save smoothed background of SCD difference
        dst_smooth = dst.createVariable('scd_difference_smoothed',np.float32,('xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
        dst_smooth.title = '{} SCD difference smoothed (observed - a priori)'.format(corrected_product)
        dst_smooth.units = units
        dst_smooth[:] = scd_diff_smooth

        # Save final destriping correction
        dst_des = dst.createVariable('destriping_correction',np.float32,('xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
        dst_des.title = '{} destriping correction'.format(corrected_product)
        dst_des.units = units
        dst_des[:] = stripe_val

        # Save number of samples used in each cross track stripe calculation
        dst_num_samples = dst.createVariable('number_samples',np.int16,('xtrack'),fill_value=-1,zlib=True,complevel=Deflate_Level)
        dst_num_samples.title = 'number of ground pixels included in each across track stripe calculation'
        dst_num_samples[:] = number_samples

except Exception as e:
    print_message(e)
    print_message('writing destriping correction to {}'.format(local_granule_id),error=True)
