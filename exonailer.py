# -*- coding: utf-8 -*-
import sys
import os
import numpy as np

sys.path.append('utilities')
import data_utils as dutils
import general_utils as gutils

################# OPTIONS ######################

# Define the target name, detrending method and parameters of it:
target = '9792'
phot_noise_model = 'white'
phot_detrend = None
window = 41

# Define if you want to perform automatic outlier removal (sigma-clipping):
phot_get_outliers = None

# Define which transits you want to ommit (counted from first transit):
n_ommit = []

# Define if you want to perform the resampling technique and in 
# which phase range you want to perform such resampling. Additionally, 
# define how many samples you want to resample:
resampling = False
phase_max = 0.01
N_resampling = 10

# Limb-darkening law to be used:
ld_law = 'quadratic'

# Define the mode to be used:
# mode = 'transit' 
mode = 'full'

# Define noise properties:
rv_jitter = True

# Define emcee parameters:
nwalkers = 50
njumps = 100
nburnin = 100

# Define time conversions:
transit_time_def = 'utc->utc'
rv_time_def = 'utc->utc'

################################################

# ---------- DATA PRE-PROCESSING ------------- #

print('Target is: %s' % target)
# First, get the transit and RV data:
output = gutils.read_data(target, mode, transit_time_def, rv_time_def)
t_tr, f, f_err, transit_instruments, t_rv, rv, rv_err, rv_instruments = output
print('Loaded data')

# Initialize the parameters:
parameters = gutils.read_priors(target, transit_instruments, rv_instruments, mode)

# Pre-process the transit data if available:
if mode != 'rvs':
    print('Pre-processing transit data...')
    t_tr, phases, f, f_err = dutils.pre_process(t_tr, f, f_err, phot_detrend,
                                                phot_get_outliers, n_ommit,
                                                window, parameters, ld_law, mode)
    if resampling:
        # Define indexes between which data will be resampled:
        idx_resampling = np.where((phases>-phase_max) & (phases<phase_max))[0]
    else:
        idx_resampling = []


# Create results folder if not already created:
if not os.path.exists('results'):
  print ('Creating "results" folder')
  os.mkdir('results')


# If chains not ran, run the MCMC and save results:
if os.path.exists('results/'+target+'_'+mode+'_'+phot_noise_model+'_'+ld_law) \
  and not '-f' in sys.argv:
    parameters = gutils.read_results(target, mode, phot_noise_model, ld_law,
                                     transit_instruments, rv_instruments)
else:
    print('Sampling...')
    dutils.exonailer_mcmc_fit(t_tr, f, f_err, transit_instruments, t_rv, rv, rv_err, rv_instruments,
                              parameters, ld_law, mode, rv_jitter=rv_jitter,
                              njumps=njumps, nburnin=nburnin,
                              nwalkers=nwalkers, noise_model=phot_noise_model,
                              resampling=resampling, idx_resampling=idx_resampling,
                              N_resampling=N_resampling)

    gutils.save_results(target, mode, phot_noise_model, ld_law, parameters)

# Get plot of the transit-fit:
if mode == 'transit':
    dutils.plot_transit(t_tr, f, parameters, ld_law, transit_instruments)
elif mode == 'full':
    dutils.plot_transit_and_rv(t_tr,f,t_rv,rv,rv_err,parameters,ld_law,rv_jitter,
                               transit_instruments, rv_instruments,
                               resampling=resampling, phase_max=phase_max,
                               N_resampling=N_resampling)

print 
posterior_file = os.path.join('results', 
                              target+'_'+mode+'_'+phot_noise_model+'_'+ld_law,
                              'posterior_parameters.dat')
with open(posterior_file) as f:
  print f.read()