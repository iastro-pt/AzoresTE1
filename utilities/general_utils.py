# -*- coding: utf-8 -*-
import numpy as np

def read_priors(target,all_transit_instruments,all_rv_instruments,mode,filename = None):
    def generate_parameter(values):
        out_dict = {}
        out_dict['type'] = values[1]
        if values[1] == 'Normal':
           out_dict['object'] = normal_parameter(np.array(values[2].split(',')).astype('float64'))
        elif values[1] == 'Uniform':
           out_dict['object'] = uniform_parameter(np.array(values[2].split(',')).astype('float64'))
        elif values[1] == 'Jeffreys':
           out_dict['object'] = jeffreys_parameter(np.array(values[2].split(',')).astype('float64'))
        elif values[1] == 'FIXED':
           out_dict['object'] = constant_parameter(np.array(values[2].split(',')).astype('float64')[0])
        if len(values)>=4:
           out_dict['object'].set_value(np.float(values[3]))
        return out_dict
    # Check instrument names of transit and rv measurements:
    if mode != 'rvs':
        transit_instruments = get_instruments(all_transit_instruments)
    if mode != 'transit':
        rv_instruments = get_instruments(all_rv_instruments)

    # Open the file containing the priors:
    if filename is None:
        f = open('priors_data/'+target+'_priors.dat','r')
    else:
        f = open(filename)
    # Generate dictionary that will save the data on the priors:
    priors = {}
    while True:
        line = f.readline()
        if line == '':
            break
        elif line[0] != '#':
            # Extract values from text file: [0]: parameter name,
            #                                [1]: prior type,
            #                                [2]: hyperparameters,
            #                                [3]: starting value (optional)
            values = line.split()
            # If more than one instrument on the RVs, and if not defined by the user, 
            # create a different parameter for both the jitter and the center-of-mass velocity:
            if mode != 'transit':
                if (values[0] == 'mu' or values[0] == 'sigma_w_rv') and len(rv_instruments)>1:
                    for instrument in rv_instruments:
                        priors[values[0]+'_'+instrument] = generate_parameter(values)
                else:
                    priors[values[0]] = generate_parameter(values)
            else:
                priors[values[0]] = generate_parameter(values)
    f.close()
    return priors

def get_instruments(instrument_list):
    all_instruments = []
    for instrument in instrument_list:
        if instrument not in all_instruments:
            all_instruments.append(instrument)
    return all_instruments

from astropy.time import Time as APYTime
def convert_time(conv_string,t):
    input_t,output_t = conv_string.split('->')
    if input_t != output_t:
        tobj = APYTime(t, format = 'jd', scale = input_t)
        exec 'new_t = tobj.'+output_t+'.jd'
        return new_t
    else:
        return t

def read_data(target,mode,transit_time_def,rv_time_def):
    t_tr,f,f_err,transit_instruments = None,None,None,None
    t_rv,rv,rv_err,rv_instruments = None,None,None,None
    if mode != 'rvs':
        # Read in transit data:
        transit_data = np.genfromtxt('transit_data/'+target+'_lc.dat',dtype='|S15')
        # Get times and fluxes:
        t_tr = transit_data[:,0].astype('float')
        f = transit_data[:,1].astype('float')
        # If more than four columns, assume fourth is instrument name, if only three, 
        # assume only errors are passed. If only two, fill instrument name with generic name:
        if transit_data.shape[1]>=4:
            f_err = transit_data[:,2].astype('float')
            transit_instruments = transit_data[:,3]
        elif transit_data.shape[1] == 3:
            f_err = transit_data[:,2].astype('float')
            transit_instruments = np.array(len(t_tr)*['instrument'])
        else:
            transit_instruments = np.array(len(t_tr)*['instrument'])
        # Convert transit times (if input and output are the same, does nothing):
        t_tr = convert_time(transit_time_def,t_tr)
    if 'transit' not in mode:
        # Read in RV data:
        rv_data = np.genfromtxt('rv_data/'+target+'_rvs.dat',dtype='|S15')
        # Get times and RVs:
        t_rv = rv_data[:,0].astype('float')
        rv = rv_data[:,1].astype('float')
        # If more than four columns, assume fourth is instrument name, if only three,
        # assume only errors are passed. If only two, fill instrument name with generic name:
        if rv_data.shape[1]>=4:
            rv_err = rv_data[:,2].astype('float')
            rv_instruments = rv_data[:,3]
        elif rv_data.shape[1] == 3:
            rv_err = rv_data[:,2].astype('float')
            rv_instruments = np.array(len(t_rv)*['instrument'])
        else:
            rv_instruments = np.array(len(t_rv)*['instrument'])
        # Convert RV times:
        t_rv = convert_time(rv_time_def,t_rv)
    return t_tr,f,f_err,transit_instruments,t_rv,rv,rv_err,rv_instruments

import pickle,os
def save_results(target,mode,phot_noise_model,ld_law,parameters):
    out_dir = 'results/'+target+'_'+mode+'_'+phot_noise_model+'_'+ld_law+'/'
    os.mkdir(out_dir)
    # Copy used prior file to the results folder:
    os.system('cp priors_data/'+target+'_priors.dat '+out_dir+'priors.dat')
    out_posterior_file = open(out_dir+'posterior_parameters.dat','w')
    out_posterior_file.write('# This file has the final parameters obtained from the MCMC chains.\n')

    # Generate an output dictionary with the posteriors:
    out_dict = {}
    for parameter in parameters.keys():
        if parameters[parameter]['type'] != 'FIXED' and len(parameters[parameter]['object'].posterior)>0:
            # Save parameter values in posterior file:
            param = parameters[parameter]['object'].value
            up_error = parameters[parameter]['object'].value_u-param
            low_error = param-parameters[parameter]['object'].value_l
            out_dict[parameter] = parameters[parameter]['object'].posterior
        else:
            param = parameters[parameter]['object'].value
            up_error = 0
            low_error = 0

        out_posterior_file.write('{0:10}  {1:10.10f}  {2:10.10f}  {3:10.10f}\n'.format(\
                                   parameter, param, up_error, low_error))
    # Save posterior dict:
    f = open(out_dir+'posteriors.pkl','w')
    pickle.dump(out_dict,f)
    f.close()

def read_results(target,mode,phot_noise_model,ld_law,all_transit_instruments,all_rv_instruments):
    out_dir = 'results/'+target+'_'+mode+'_'+phot_noise_model+'_'+ld_law+'/'
    parameters = read_priors(target,all_transit_instruments,all_rv_instruments,mode,filename = out_dir+'priors.dat')
    thefile = open(out_dir+'posteriors.pkl','r')
    posteriors = pickle.load(thefile)
    for parameter in parameters.keys():
        if parameters[parameter]['type'] != 'FIXED':
            try:
                parameters[parameter]['object'].set_posterior(posteriors[parameter])  
            except:
                print 'No posterior for parameter '+parameter
    thefile.close()
    return parameters

def get_quantiles(dist,alpha = 0.68, method = 'median'):
    """
    get_quantiles function

    DESCRIPTION

        This function returns, in the default case, the parameter median and the error% 
        credibility around it. This assumes you give a non-ordered 
        distribution of parameters.

    OUTPUTS

        Median of the parameter,upper credibility bound, lower credibility bound

    """
    ordered_dist = dist[np.argsort(dist)]
    param = 0.0
    # Define the number of samples from posterior
    nsamples = len(dist)
    nsamples_at_each_side = int(nsamples*(alpha/2.)+1)
    if(method == 'median'):
       med_idx = 0
       if(nsamples%2 == 0.0): # Number of points is even
          med_idx_up = int(nsamples/2.)+1
          med_idx_down = med_idx_up-1
          param = (ordered_dist[med_idx_up]+ordered_dist[med_idx_down])/2.
          return param,ordered_dist[med_idx_up+nsamples_at_each_side],\
                 ordered_dist[med_idx_down-nsamples_at_each_side]
       else:
          med_idx = int(nsamples/2.)
          param = ordered_dist[med_idx]
          return param,ordered_dist[med_idx+nsamples_at_each_side],\
                 ordered_dist[med_idx-nsamples_at_each_side]

class normal_parameter:
      """
      Description
      -----------

      This class defines a parameter object which has a normal prior. It serves 
      to save both the prior and the posterior chains for an easier check of the parameter.

      """   
      def __init__(self,prior_hypp):
          self.value = prior_hypp[0]
          self_value_u = 0.0
          self_value_l = 0.0
          self.prior_hypp = prior_hypp
          self.posterior = []

      def get_ln_prior(self):
          return np.log(1./np.sqrt(2.*np.pi*(self.prior_hypp[1]**2)))-\
                 0.5*(((self.prior_hypp[0]-self.value)**2/(self.prior_hypp[1]**2)))

      def set_value(self,new_val):
          self.value = new_val

      def set_posterior(self,posterior_chain):
          self.posterior = posterior_chain
          param, param_u, param_l = get_quantiles(posterior_chain)
          self.value = param
          self.value_u = param_u
          self.value_l = param_l
class uniform_parameter:
      """
      Description
      -----------

      This class defines a parameter object which has a uniform prior. It serves 
      to save both the prior and the posterior chains for an easier check of the parameter.

      """
      def __init__(self,prior_hypp):
          self.value = (prior_hypp[0]+prior_hypp[1])/2.
          self_value_u = 0.0
          self_value_l = 0.0
          self.prior_hypp = prior_hypp
          self.posterior = []

      def get_ln_prior(self):
          return np.log(1./(self.prior_hypp[1]-self.prior_hypp[0]))

      def check_value(self,x):
          if x > self.prior_hypp[0] and  x < self.prior_hypp[1]:
              return True
          else:
              return False  
 
      def set_value(self,new_val):
          self.value = new_val

      def set_posterior(self,posterior_chain):
          self.posterior = posterior_chain
          param, param_u, param_l = get_quantiles(posterior_chain)
          self.value = param
          self.value_u = param_u
          self.value_l = param_l

log1 = np.log(1)
class jeffreys_parameter:
      """
      Description
      -----------

      This class defines a parameter object which has a Jeffreys prior. It serves 
      to save both the prior and the posterior chains for an easier check of the parameter.

      """
      def __init__(self,prior_hypp):
          self.value = np.sqrt(prior_hypp[0]*prior_hypp[1])
          self_value_u = 0.0
          self_value_l = 0.0
          self.prior_hypp = prior_hypp
          self.posterior = []

      def get_ln_prior(self):
          return log1 - np.log(self.value*np.log(self.prior_hypp[1]/self.prior_hypp[0]))

      def check_value(self,x):
          if x > self.prior_hypp[0] and  x < self.prior_hypp[1]:
              return True
          else:
              return False

      def set_value(self,new_val):
          self.value = new_val

      def set_posterior(self,posterior_chain):
          self.posterior = posterior_chain
          param, param_u, param_l = get_quantiles(posterior_chain)
          self.value = param
          self.value_u = param_u
          self.value_l = param_l

class constant_parameter:
      """
      Description
      -----------

      This class defines a parameter object which has a constant value. It serves 
      to save both the prior and the posterior chains for an easier check of the parameter.

      """
      def __init__(self,val):
          self.value = val

