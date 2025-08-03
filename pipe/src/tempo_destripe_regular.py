#! /usr/bin/env python3
######################################
# purpose: destriping TEMPO L2 CLDO4
#    run between O2O2 fitting and CLDO4 retrieval
#
# a simplified version of tempo_cldo4/tempo_destripe_cldo4.py
#
# history:
#    hwang Oct 2023 initial core code 1D descor
#    hwang Mar 2025 2D desarr based on scd
#
######################################
import argparse
import sys
import os
import random
import numpy as np
import numpy.ma as ma
from pathlib import Path
from netCDF4 import Dataset
#from matplotlib import pyplot as plt
#import pdb

######################
# global parameters for O2-O2
######################
nXtrack = 2048 # number of cross-track pixels
max_scd = 1.0e+44 # max o2o2 scd to include for descor
min_scd = 5.e+41 # min o2o2 scd to include for descor
tiny_scd = 1.e+42 # scd<tinyscd only retained when rms<rms4tiny
big_scd = 4.e+43 #  scd>bigscd only retained when rms<rms4tiny
rms4tiny = 0.01 # rms threshold for excluding tinyscd & bigscd during application

max_rms = 0.0025 #0.005 # maxrms for sza=40 to include when deriving descor
                # will change with granule median sza
# poly3 coeffs for 75th rms percentile derived from 20230819
# used to adjust maxrms when derive descor
c0 = 0.001892
c1 = -3.34829e-06
c2 = -4.87580e-07
c3 = 7.66683e-09
r40 = c0 + c1*40. + c2*1600.+c3*64000. # rms at sza=40 from coeffs

min_sza = 0. # minimum solar_zenith_angle to include for descor
max_sza = 89. # maximum solar zenith angle to include for descor

normcol = 1.e43 # normalization factor for o2o2
missval = -1.0e+30 # missing values
width_savgol = 81 # wider stripes need larger width_savgol, 35 works for most
order_savgol = 3 # polynomial order for savgol
# of descor pixels to set to 1 at end1 & end2
# these end points are not destriped
n_end1 = 10  # south end
n_end2 = n_end1 # north end

# pixels to mirror extend for savgol fitting
n_mirror = width_savgol*2

# limit for descor amplitude
# larger descor_nsigma allow more pixels to be corrected
# huge descor_nsigma sends control to deslimit_low & deslimit_high
# smaller descor_nsigma provides tighter control on correction limit
descor_nsigma = 10.0 # huge
deslimit_low = 0.90 # 0.50
deslimit_high = 1.10 # 2.0

# descor outside [deslimit_low,deslimit_high] will be set
# to (oob_void=True) 1 or (oob_void=False) deslimit
# True:no correction, False: correct up to deslimit
oob_void = False #True

# correction method: divide//subtract//smaller//larger//average
# smaller & larger refers to the amount of correction
# average use the mean of divide & subtract
c_method = 'average'

# whether remove residual outlier from descor
rm_outlier = 'y'

#######################
# derive descor vector (core program)
#######################
def derive_descor(input_files,nx=nXtrack,fout=None,
                 maxrms=max_rms,minscd=min_scd,maxscd=max_scd,
                 maxsza=max_sza,minsza=min_sza,nmask_pct=0.95,
                 rm_suspicious_column=True,rms_nsigma = 3.,
                 savgol_width=width_savgol,savgol_order=order_savgol,
                 nend1=n_end1,nend2=n_end2,rm_outlier=rm_outlier,
                 missing=missval,descor_nsigma=descor_nsigma,
                 deslimit_low=deslimit_low,deslimit_high=deslimit_high,
                 nmirror=n_mirror,oob_void=oob_void):

    from scipy.signal import savgol_filter
    from scipy.interpolate import interp1d
    from tempo_cldo4 import tempo_cldo4_util as tcu

    maxscdnorm = maxscd/normcol
    minscdnorm = minscd/normcol

    #initialize masked arrays
    rms = ma.array(np.zeros([1,nx],dtype=float),mask=True)
    scd = ma.array(np.zeros([1,nx],dtype=float),mask=True)

    nfiles = 0 # count number of files in input_files

    for fp in input_files:
        #print(fp)
        try:
            # read data from fp
            with Dataset(fp,'r') as src:
                # note scd is normalized here
                s00 = src['support_data']['fitted_slant_column'][:]/normcol
                q = src['qa_statistics']['fit_convergence_flag'][:] #1=good
                r = src['qa_statistics']['fit_rms_residual'][:]
                theta = src['geolocation']['solar_zenith_angle'][:]

            # remove suspicious scd column regardless of rm_suspicious_column
            # to ensure they do not mess up statistics
            s = np.copy(s00)
            if rm_suspicious_column:
               s = remove_suspicious_scdcolumn(s00,missing=missing,percent=0.4,
                   max_value=maxscdnorm,silent=True)

            # remove out of range values
            s[s00<minscdnorm] = missing
            s[s00>maxscdnorm] = missing
            # remove sza outside the requested range
            s[(theta > maxsza)|(theta < minsza)] = missing

            # remove non-convergence
            s[q != 1] = missing

            # rms depends on sza
            # find max rms limit for granule median sza
            # using fitted poly3 from tempo_exam_rmssza for 20230819
            thmed = ma.median(theta)
            rmed = c0 + c1 * thmed + c2 * thmed * thmed + \
                    c3 * thmed * thmed * thmed
            rmax = rmed + (maxrms - r40) # maxrms is for sza=40.
            rmedian = ma.median(r)
            rsigma = ma.median(ma.abs(rms - rmedian))
            rlim = rmedian + rsigma*rms_nsigma
            if rlim > rmax:
               rlim = rmax
            # remove large rms
            s[(r > rlim) | (r <= 0.)] = missing

            nfiles = nfiles + 1
        except Exception as e:
            print('reading from {0}'.format(fp),error=True)
            continue # to next file

        # Append arrays
        rms = ma.concatenate([rms,r],axis=0)
        scd = ma.concatenate([scd,s],axis=0)

        del q,r,s,s00,theta,thmed,rmed,rmax,rmedian,rsigma,rlim

    print('scd shape:',scd.shape)

    # Get RMS mask
    # calculate limit for outliers of fitting rms
    # larger of rms_nsigma * msig & maxrms is used to set rmslim
    # though rms_mask may improve results for low & moderate sza
    # large sza will be preferentially filtered out below
    # note rms filtering was applied during data accumulation
    # thus, rms_mask serves as an overall filtering when needed
    msig = ma.median(ma.abs(rms - ma.median(rms)))
    rmslim = ma.median(rms) + rms_nsigma*msig
    rms_mask = (rms > rmslim) | (rms <= 0.) # anomalous rms

    # Get scd mask
    scd_mask = (scd < minscdnorm) | (scd > maxscdnorm) # bad scd

    # Mask out bad overall pixels, uncomment if needed
    # scd[rms_mask] = ma.masked # bad overall
    # as rlim is considered before, we only apply scd_mask
    scd[scd_mask] = ma.masked  # mask out scd

    # number of masked points which are not used for ma.median
    count_masked = ma.count_masked(scd, axis=0)
    count_allpix = scd.shape[0]

    # maximum number of masked out pixels for each column
    # to ensure that at least 5% or 100 are left for statistics
    # note 1-0.95=5% left for calc works for 1d descor
    # for 2d desarr, 5% may cause many xtrack to filter out
    # in this case, increase nmask_pct, e.g., from 0.95 to 0.99
    max_number_masked0 = count_allpix * nmask_pct #percent masked out
    # ensure at least 100 data points for subsequent calculation
    max_number_masked1 = count_allpix - 100 # 100 points left
    # max_number_masked will be used for filtering
    max_number_masked = max_number_masked1 # max_number_masked0

    print('max_number_masked =',max_number_masked)
    print('count_masked:',count_masked)

    # median value of scds for de-striping corrections.
    # Calculate medians for each xtrack position
    scd_med  = ma.median(scd,axis=0)
    print('scd_med shape:',scd_med.shape)

    # mask pixel with too many filtered values
    medval = scd_med
    medval = np.ma.masked_where(count_masked>max_number_masked,scd_med)

    # number of masked medval
    count_medmasked = ma.count_masked(medval)
    if (count_medmasked > 512): # too many masked medval
       # cannot derive useful descor, assign all to 1, return
       descorfinal = np.ones(nx)
       return descorfinal

    #save a copy for use later
    medval00 = medval.copy()

    #mirror extend nmirror pixels on both ends
    #to prevent large anomaly at edge
    #1st & last few TEMPO pixels are always bad
    medval_flip = np.flip(medval00)
    medval_ext = ma.append(ma.append(medval_flip[0:nx-1],medval00),medval_flip[1:nx])

    # replace masked value in medval_ext with interpolated value
    # otherwise, masked value can mess up savgol_filter result
    filled_array = np.copy(ma.getdata(medval_ext))
    valid_indices = np.where(~medval_ext.mask)[0]
    valid_data = medval_ext[valid_indices]
    vfunct = interp1d(valid_indices,valid_data,kind='linear',
                      fill_value='extrapolate')
    masked_indices = np.where(medval_ext.mask)[0]
    filled_array[masked_indices] = vfunct(masked_indices)

    nmirror1 = nmirror + 1
    filled_use = filled_array[nx-nmirror1:nx-nmirror1+nx+2*nmirror]
    nx_use = len(filled_use)
    x_use = np.arange(nx_use)

    # savgol filter
    fillsavgol = savgol_filter(filled_use,savgol_width,savgol_order)
    #pdb.set_trace()

    #select the non-extension part of medvalfit
    fill_savgol = fillsavgol[nmirror:nx+nmirror]

    # initial descor through ratio
    descor11 = np.divide(medval00, fill_savgol)

    # set filled_value to 1.0
    descor11 = descor11.filled(fill_value=1.0)

    #further adjust correction limits
    descor_mean = np.mean(descor11)
    descor_stdev = np.std(descor11)
    # anything < des1 will be set to -1 below
    des1 = descor_mean - descor_stdev*descor_nsigma
    if des1 < deslimit_low:
       des1 = deslimit_low
    # anything > des2 will be set to -1 below
    des2 = descor_mean + descor_stdev*descor_nsigma
    if des2 > deslimit_high:
        des2 = deslimit_high
    #print('des1,des2=',des1,des2)

    descorsavg = descor11
    if (oob_void): # set extremes to -1
       descorsavg[descor11<des1] = -1
       descorsavg[descor11>des2] = -1
    else: # set extreme to des1 or des2
       descorsavg[descor11<des1] = des1
       descorsavg[descor11>des2] = des2

    # negative values set to 1.0 (no correction)
    descorsavg[descorsavg<0.] = 1.0

    # first & last nend points set to 1.0 (no correction)
    # they may be affected by edges
    descorsavg[0:nend1] = 1.0
    descorsavg[-nend2:] = 1.0

    # remove missed outliers
    # descorsavg is typically at 5% level here
    # outliers will leave sigular line in result
    if (rm_outlier == 'y'):
       descorfinal = remove_descor_outlier(descorsavg,
           specialvalue=1.0,nsigma=2.,nptmax=2,nptmin=2)
    else:
       descorfinal = np.copy(descorsavg)

    #print(descorfinal)
    #pdb.set_trace()

    # write descor to file according to fout extension
    if (fout is not None):
        write_descor1d_result(fout,descorfinal,medval00,fill_savgol)

    return descorfinal

######
# derive descor from fnms in list
def derive_descor_from_list(listnm,minscd=min_scd,maxscd=max_scd,
    minsza=min_sza,maxsza=max_sza,fout='descor.nc',
    missing=missval,nmirror=n_mirror,nx=nXtrack,maxrms=max_rms,
    descor_nsigma=descor_nsigma,
    deslimit_low=deslimit_low,deslimit_high=deslimit_high,
    nend1=n_end1,nend2=n_end2,maxnfile=500):

    input_files = get_fnms_from_list(listnm)
    n1 = len(input_files)
    if n1 ==  0:
        # all ones, no destriping due to empty listnm
        descor = np.ones(nXtrack)
    else:
        if n1 > maxnfile: # for memory and speed
           input_fnms = random.sample(input_files,maxnfile)
        else:
           input_fnms = input_files

        # derive overall descor using all input_files in listnm
        descor = derive_descor(input_fnms,maxrms=maxrms,
           minscd=minscd,maxscd=maxscd,nend1=nend1,nend2=nend2,
           maxsza=maxsza,minsza=minsza,descor_nsigma=descor_nsigma,
           deslimit_low=deslimit_low,deslimit_high=deslimit_high,
           nmirror=nmirror,nx=nx,missing=missing,fout=fout)

    return descor

#####################################################################
# application
#####################################################################
#
# apply single destriping factor in descorfnm to a L2 file
def apply_single_descor(l2fnm,descorfnm,missing=missval,maxscd=max_scd,
        tinyscd=tiny_scd,bigscd=big_scd,rms4tiny=rms4tiny,
        minsza=min_sza,maxsza=max_sza,method=c_method,
        write_new='n',dirout='./',write_update='n'):

    # method = divide//subtract//smaller//larger
    # only minsza<sza<maxsza will be corrected
    # for tiny & big scd, only rms<rms4tiny will be retained

    # read l2fnm
    with Dataset(l2fnm,'r') as src:
       scd = src['support_data']['fitted_slant_column'][:]
       rms = src['qa_statistics']['fit_rms_residual'][:]
       sza = src['geolocation']['solar_zenith_angle'][:]

    (nt, nx) = scd.shape

    # read descor file
    descor, medval00, medval_fit = read_descor_result(descorfnm,
                                   nx=nx)
    # medval_fit has been normalized, scd is not normalized

    # init scddes,scddes1,scddes2
    # multiplicative
    scddes1 = np.zeros(scd.shape,dtype='float64')
    scddes1[:,:] = missing
    # additive
    scddes2 = np.zeros(scd.shape,dtype='float64')
    scddes2[:,:] = missing
    # result
    scddes = np.copy(scd)

    # if max(descor) == 0, nothing to be done
    if np.max(descor) == 0:
        print('invalid descor')
        return scddes

    #pdb.set_trace()
    # apply descor
    for ix in np.arange(nx):
        thisfact = descor[ix]
        thisadd = medval_fit[ix]*(1.-thisfact)
        if (thisfact > 0.):
            scddes1[:,ix] = scd[:,ix]/thisfact
            scddes2[:,ix] = scd[:,ix]-scd[:,ix]*(thisfact-1.)
            # assign scddes to the one with smaller change
            for it in np.arange(nt):
                thisscd = scd[it,ix]
                if (~ma.is_masked(thisscd)):
                   change1 = np.abs(scddes1[it,ix]-thisscd)
                   change2 = np.abs(scddes2[it,ix]-thisscd)
                   if (method == 'divide'):
                      scddes[it,ix] = scddes1[it,ix]
                   elif (method == 'subtract'):
                      scddes[it,ix] = scddes2[it,ix]
                   elif (method == 'smaller'):
                      if (change1 > change2):
                         scddes[it,ix] = scddes2[it,ix]
                      else:
                         scddes[it,ix] = scddes1[it,ix]
                   elif (method == 'larger'):
                      if (change1 > change2):
                         scddes[it,ix] = scddes1[it,ix]
                      else:
                         scddes[it,ix] = scddes2[it,ix]
                   else: # average
                      scddes[it,ix] = 0.5*(scddes1[it,ix]+scddes2[it,ix])

    # set out-of-range szas to scd
    ind1 = sza < minsza
    ind2 = sza > maxsza
    scddes[ind1] = scd[ind1]
    scddes[ind2] = scd[ind2]

    #set non-physical values to missing
    indbad=(scd<0.)|(scd>maxscd)
    scddes[indbad]=missing
    scd[indbad] = missing

    #restrict tinyscd & bigscd only to acceptable rms
    indbad2=(scd<tinyscd) & (rms>rms4tiny)
    scddes[indbad2] = missing
    scd[indbad2] = missing
    indbad2=(scd>bigscd) & (rms>rms4tiny)
    scddes[indbad2] = missing
    scd[indbad2] = missing

    # void pixels with invalid sza
    indbad3 = (sza < 0.) | (sza > 90.)
    scddes[indbad3] = missing
    scd[indbad3] = missing

    # add or update l2fnm (for sdpc)
    if (write_update == 'y'):
       tempo_addorupdate_2dfloat(l2fnm,'support_data','scddes',scddes,
           valid_min=0.,valid_max=1.e45,long_name='destriped column amount',
           units='molec^2 cm^-5',comment='reglar destriping')
       #tempo_addorupdate_1dfloat(l2fnm,'support_data','descor',descor,
       #    valid_min=-2.,valid_max=2.)

    # write to new nc file (for testing)
    if (write_new == 'y'):
        newfnm = dirout+'/'+Path(l2fnm).stem+'_new.nc'
        create_scddes_file(newfnm,scd,scddes,descor)

    return scddes, scd

######
# apply 1d descor in descorfnm to l2 files in l2list
def apply_descor1d_list(l2list,descorfnm,missing=missval,maxscd=max_scd,
        tinyscd=tiny_scd,bigscd=big_scd,rms4tiny=rms4tiny,
        minsza=min_sza,maxsza=max_sza,method=c_method,
        dirout='./',write_new='y',write_update='n'):
    # write_update='y' will update original l2 file
    # write_new ='y' will write new file
    with open(l2list,'r') as src:
        for line in src:
            fnm = line.strip()
            print(fnm)
            apply_single_descor(fnm,descorfnm,missing=missing,
                maxscd=maxscd,tinyscd=tinyscd,dirout=dirout,method=method,
                minsza=minsza,maxsza=maxsza,bigscd=bigscd,rms4tiny=rms4tiny,
                write_new=write_new,write_update=write_update)

    return

########################################################################
# utility
########################################################################
def remove_suspicious_scdcolumn(scd,percent=0.4,min_value=0.,max_value=max_scd,
                  missing=missval,silent=False):
# columns with count(scd<min_value | scd>max_value) > percent will be set to missing
# columns with median beyond normal median threshold will also be set to missing
# the columns are likely affected by bad spectral pixels at critical wavelength
# which results in abnormal fit that is not caught by fitting rms
# because the fitting performs outlier rejection, e.g. 20050710
    ss = scd.shape
    nx = ss[1]
    ny = ss[0]
    nthresh = int(ny * percent) # threshold for number of invalid points

    scdmed = np.zeros((nx))
    carr = np.zeros((nx))

    for ix in np.arange(nx):
        acol = scd[:,ix].flatten()
        # replace masked data with fill_value
        bcol = acol.filled(fill_value=missing)
        # number of masked invalid points for each column
        carr[ix] = ((bcol < min_value)|(bcol > max_value)).sum()

        # median scd of each column
        # using only unmasked points
        am = ma.median(acol)
        if ma.is_masked(am):
           am = missing
        # if all masked, am=masked will be converted to nan
        scdmed[ix] = am

    # in case of nan, replace with missing
    trash = scdmed
    trash[np.isnan(scdmed)]=missing
    scdmed = trash

    # anything < scdthresh1 or > scdthresh2 are unlikely
    scdnormal = np.median(scdmed)
    scdthresh1 = scdnormal* 0.1
    scdthresh2 = scdnormal* 10.0
    if scdthresh2 > max_value:
        scdthresh2 = max_value

    # set invalid column to missing
    brrout = scd
    missing_idx = []
    for ix in np.arange(nx):
        cc = carr[ix]
        yy = scdmed[ix]
        if ((yy < scdthresh1) | (yy>scdthresh2) | (cc >nthresh)):
            brrout[:,ix] = missing
            missing_idx.append(ix) # append invalid column index

    if silent != True:
       print('column idx  set to missing:',missing_idx)
       print('scdthresh:',scdthresh1,scdthresh2)
       print('scdmed:',scdmed[missing_idx])
       print('nthresh=',nthresh)
       print('carr:',carr[missing_idx])

    return brrout

######
# remove outlier from a descor vector
# if abs(max-5thmax) > nsigma or abs(5thmin - min) > nsigma
# the max or min position will be set to special value
def remove_descor_outlier(arrin,nsigma=2,specialvalue=1.,
    nptmax=3,nptmin=3):

    # make a copy of arrin
    arrout = np.copy(arrin)
    nx = len(arrin)

    sigma = np.std(arrin)
    dthresh = nsigma * sigma

    # index for sorted arrin increasing order
    sindex = arrin.argsort()

    # steep drop between within 5 pixels of sorted indicate outlier
    # test nptmax points
    for i in np.arange(nptmax):
       id1 = sindex[nx-i-1]
       id5 = sindex[nx-i-5]
       dmax = np.abs(arrin[id1] - arrin[id5])
       if (dmax > dthresh):
          arrout[id1] = specialvalue

    # test nptmin points:
    for i in np.arange(nptmin):
       id1 = sindex[i]
       id5 = sindex[i+5]
       dmin = np.abs(arrin[id5] - arrin[id1])
       if (dmin > dthresh):
          arrout[id1] = specialvalue

    return arrout
#
########################################################################
# reading & writing
########################################################################
#
######
# get filenames from list
def get_fnms_from_list(listnm):
    filenames = []
    try:
       with open(listnm) as f:
           fnms = f.readlines()
       for fnm in fnms:
           filenames.append(fnm.strip())
    except:
       print('cannot find '+listnm)

    return filenames

######
# create ncdf variable
def create_variable(ncid,name,kind,dimensions,attributes,fill_value=None):
    ''' Create variable in ncid netCDF file
        ARGS:
            ncid (netCDF4 object)
            name (string): variable name
            kind: variable type
            dimensions (list): variable dimensions
            attributes (dictionary): attributes of the variable
    '''

    try:
        print('adding variable {0}'.format(name))
        if fill_value:
            vid = ncid.createVariable(name,kind,dimensions=dimensions,\
                  fill_value=fill_value)
            vid.setncatts(attributes)
        else:
            vid = ncid.createVariable(name,kind,dimensions=dimensions)
            vid.setncatts(attributes)
    except Exception as e:
        print('error creating {0} variable; exception {1}'.format(name,e))

    return

######
# write 1D variable to l2 ncdf file
def tempo_addorupdate_1dfloat(filename,groupname,varname,varval,
        comment='None',long_name='None',units='None',missing=missval,
        valid_min=-1.e45,valid_max=1.e45):

    with Dataset(filename,'r+') as dst:
        al = dst.dimensions['mirror_step'].size
        ax = dst.dimensions['xtrack'].size
        if varname not in dst[groupname].variables:
           print('add '+varname+' to '+groupname)
           attributes={'long_name':long_name,'comment':comment,\
                'valid_min':valid_min,'units':units,\
                'valid_max':valid_max,\
                }
           create_variable(dst.groups[groupname],varname,np.float64,\
                ('xtrack'),attributes,
                fill_value=missing)
        else:
           print('update '+varname+' in '+groupname)

        # write data
        dst[groupname][varname][:] = varval

        # write attributes
        dst[groupname][varname].long_name = long_name
        dst[groupname][varname].comment = comment
        dst[groupname][varname].units = units
        dst[groupname][varname].valid_min = valid_min
        dst[groupname][varname].valid_max = valid_max
    return

######
# write 2D variable to level 2 ncdf file
def tempo_addorupdate_2dfloat(filename,groupname,varname,varval,
        missing=missval,comment='None',long_name='None',units='None',
        valid_min=float(-1.e45),valid_max=float(1.e45)):

    with Dataset(filename,'r+') as dst:
        al = dst.dimensions['mirror_step'].size
        ax = dst.dimensions['xtrack'].size
        if varname not in dst[groupname].variables:
           print('add '+varname+' to '+groupname)
           attributes={'long_name':long_name,'comment':comment,\
                'valid_min':valid_min,'units':units,\
                'valid_max':valid_max,\
                }
           create_variable(dst.groups[groupname],varname,np.float64,\
                ('mirror_step','xtrack'),attributes,
                fill_value=missing)
        else:
           print('update '+varname+' in '+groupname)

        # write attributes
        dst[groupname][varname].long_name = long_name
        dst[groupname][varname].comment = comment
        dst[groupname][varname].units = units
        dst[groupname][varname].valid_min = valid_min
        dst[groupname][varname].valid_max = valid_max

        # write data
        dst[groupname][varname][:,:] = varval

    return

######
# create diagnostic descor file
def create_descor_diagfile(diagfnm,medval00,medval_fit,descor,nx=nXtrack):

    #Create a file to save the results of the calculation
    print('writting correction results to {0}'.format(diagfnm))
    try:
        with Dataset(diagfnm,'w',clobber=True) as dst:

            # Create dimensions
            dst_nx = dst.createDimension('cross_track',nx)

            # Create variables
            print('descor')
            dst_des = dst.createVariable('descor',np.float64,('cross_track'),
                      fill_value=-1.0e30,zlib=True,complevel=4)
            dst_des.title = 'final destriping vector (divide)'
            dst_des.units = '1'
            dst_des[:] = descor

            print('medval00')
            dst_medval = dst.createVariable('medval00',np.float64,('cross_track'),
                      fill_value=-1.0e30,zlib=True,complevel=4)
            dst_medval.title = 'median scd used for descor00'
            dst_medval.units = 'molec cm^-2'
            dst_medval[:] = medval00

            print('medval_fit')
            dst_medfit = dst.createVariable('medval_fit',np.float64,('cross_track'),
                      fill_value=-1.e30,zlib=True,complevel=4)
            dst_medfit.title = 'fit to medval00'
            dst_medfit.units = 'molec cm^-2'
            dst_medfit[:] = medval_fit

    except Exception as e:
        print('writting correction results to {0}'.format(diagfnm),error=True)

    return

######
# write txt or nc descor file
def write_descor1d_result(fout,descor,medval00,medval_fit):
    print('writing '+fout)

    # convert masked array to numpy array
    descor2 = ma.getdata(descor)
    medval002 = ma.getdata(medval00)
    medval_fit2 = ma.getdata(medval_fit)

    # find filename extension
    extension = os.path.splitext(fout)[1][1:]

    if (extension == 'txt'): # txt
       with open(fout,'w') as f:
          for value in descor2:
             f.write(str(value)+'\n')
          for a in medval002:
             f.write(str(a)+'\n')
          for b in medval_fit2:
             f.write(str(b)+'\n')
    else: # nc
       create_descor_diagfile(fout,medval002,medval_fit2,descor2)

    return

######
#
def read_descor_result(fin,nx=nXtrack):
    extension = os.path.splitext(fin)[1][1:]

    descor = np.ones(nx)
    medval00 = np.ones(nx)
    medval_fit = np.ones(nx)

    if (extension == 'txt'): #txt
       try:
           with open(fin,'r') as f:
              for i in np.arange(nx):
                 thisline = f.readline()
                 descor[i] = float(thisline.strip())
              for i in np.arange(nx):
                 thisline = f.readline()
                 medval00[i] = float(thisline.strip())
              for i in np.arange(nx):
                 thisline = f.readline()
                 medval_fit[i] = float(thisline.strip())
       except:
           print('error reading '+fin)
    else : #nc
       try:
           with Dataset(fin,'r') as f:
              descor = f['descor'][:]
              medval00 = f['medval00'][:]
              medval_fit = f['medval_fit'][:]
       except:
           print('error reading '+fin)

    return descor,medval00,medval_fit

######
#
def create_scddes_file(fnm,scd,scddes,descor):

    #Create a file to save scd & scddes
    print('create {0}'.format(fnm))
    (nt, nx) = scd.shape
    ndimdes = len(descor.shape)
    if ndimdes == 1:
        d1 = descor.shape
        ninterval = 1
    else:
        (d1,ninterval) = descor.shape
    if (d1 != nx):
        print('dimension difference between scd and descor')

    try:
        with Dataset(fnm,'w',clobber=True) as dst:

            # Create dimensions
            dst_nx = dst.createDimension('cross_track',nx)
            dst_nt = dst.createDimension('mirror_step',nt)
            if (ndimdes == 2):
                dst_nint = dst.createDimension('interval',ninterval)

            # Create variables
            dst_des = dst.createVariable('fitted_slant_column',np.float64,
                ('mirror_step','cross_track'),fill_value=-1.0e30,
                zlib=True,complevel=4)
            dst_des.title = 'fitted_slant_column'
            dst_des.units = 'molec^2 cm^-5'
            dst_des[:] = scd

            dst_des2 = dst.createVariable('scddes',np.float64,
                ('mirror_step','cross_track'),fill_value=-1.0e30,
                zlib=True,complevel=4)
            dst_des2.title = 'scddes'
            dst_des2.units = 'molec^2 cm^-5'
            dst_des2[:] = scddes

            if (ndimdes == 1):
               dst_des3 = dst.createVariable('descor',np.float64,
                   ('cross_track'),fill_value=-1.0e30)
               dst_des3.title = 'descor'
               dst_des3.units = '1'
               dst_des3[:] = descor
            else:
               dst_des3 = dst.createVariable('descor',np.float64,
                   ('cross_track','interval'),fill_value=-1.0e30,
                   zlib=True,complevel=4)
               dst_des3.title = 'descor'
               dst_des3.units = '1'
               dst_des3[:] = descor

    except Exception as e:
        print('writting {0}'.format(fnm),error=True)

    return

########################################################################
# MAIN progrm
########################################################################
# 1d destriping a list of L2 files
def main1d_list():
    # parse arguments
    parser = argparse.ArgumentParser(description='destripe o2o2 scd')
    parser.add_argument('--list_to_destripe',default="unknown",
                        help="list of fnms to destripe")
    parser.add_argument('--list4descor',default="unknown",
                        help="list of files for deriving descor")
    parser.add_argument('--write_new',default="y",
                        help="whether to write new output file y/n")
    parser.add_argument('--outdir',default='./',
                        help="output dir for destriped files")
    parser.add_argument('--desfnm',default='descor.txt',
                        help="descor txt filename")

    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    list_to_destripe = args.list_to_destripe
    list4descor = args.list4descor
    write_new = args.write_new
    dirout = args.outdir # only used when write_new == 'y'
    desfnm = args.desfnm
    if write_new == 'y':
        write_update = 'n'
    else:
        write_update = 'y'

    if os.path.exists(desfnm):
        # will use desfnm if it already exists
        print('will use descor in '+desfnm)
    else:
        # will derive descor
        print('will derive descor, save to '+desfnm)
        descor = derive_descor_from_list(list4descor,fout=desfnm)

    # apply descor to files in list_to_destripe
    print('apply descor to each file in '+list_to_destripe)
    fnms = get_fnms_from_list(list_to_destripe)
    n = len(fnms)
    for fnm in fnms:
        print(fnm)
        apply_single_descor(fnm,desfnm,dirout=dirout,
            write_new=write_new,write_update=write_update)

    return

########
# 1d destriping a single L2 file
# mainly for use with sdpc
def main1d_l2file():
    # parse arguments
    parser = argparse.ArgumentParser(description='destripe o2o2 scd for a L2 file')
    parser.add_argument('--l2fnm',default="none",
                        help="L2 filename to destripe")
    parser.add_argument('--list4descor',default="unknown",
                        help="list of files for deriving descor, used only when descor txt needs to be created")
    parser.add_argument('--desfnm',default='descor.nc',
                        help="descor filename to create from list4descor or use")
    parser.add_argument('--mode',default='both',
                        help="derive/apply/both")

    write_new = 'n' # hardcoded to not write new file
    dirout = './' # used only when write_new == 'y', thus, not used
    write_update = 'y' # hardcoded to update l2fnm
    fillvalue = -1.e30
    min_nfiles = 30 # minimum number of files in list4descor for it to be used'

    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    l2fnm = args.l2fnm
    list4descor = args.list4descor # used only when desfnm does not exist
    desfnm = args.desfnm
    mode = args.mode
    # in derive or both mode, list4descor will be used to derive descor
    # result will be saved to desfnm
    # in derive mode, delete existing desfnm first
    if (mode == 'derive'):
        if os.path.exists(desfnm):
            os.remove(desfnm)
            print('existing '+desfnm+' is deleted.')

    # apply mode will not run the block below
    if (mode == 'derive') | (mode == 'both'):
        print('deriving descor from '+list4descor)
        fffs = get_fnms_from_list(list4descor)
        mmm = len(fffs)
        print(f"    n_files in list4descor = {mmm}")
        if (mmm < min_nfiles): # not enough filenames in list4descor
            print(list4descor+' is too short to derive descor')
            print('   I will not create new '+desfnm)
        else: # enough filenames in list4descor
            descor = derive_descor_from_list(list4descor,fout=desfnm)
            print('   save descor to '+desfnm)

    #---
    if (mode == 'derive'):
        print('return from derive mode.')
        return

    # derive mode stops here
    # apply or both mode will continue below

    #---
    if os.path.exists(desfnm):
        # will use desfnm if it exists
        print('apply descor to '+l2fnm)
        print('    will use descor in '+desfnm)
        apply_single_descor(l2fnm,desfnm,dirout=dirout,missing=fillvalue,
            write_new=write_new,write_update=write_update)
    else: # desfnm does not exist
        print('   cannot find '+desfnm)
        # will add invalid support_data/scddes with fillvalue
        print('   add support_data/scddes with fillvalues in '+l2fnm)
        with Dataset(l2fnm,'r') as src:
            tmparr=src['support_data']['fitted_slant_column'][:]
        tmparr[:,:] = fillvalue
        tempo_addorupdate_2dfloat(l2fnm,'support_data','scddes',tmparr,
            valid_min=0.,valid_max=1.e45,units='molec^2 cm^-5',
            long_name='destriped column amount',
            comment = 'no destriping')

    return

#=====================================================================
if __name__ == '__main__':
    main1d_l2file()
