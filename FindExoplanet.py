import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import emcee
import batman
from utilities.ajplanet import pl_rv_array as rv_curve

fix_ecc = True


def get_tp(P, ecc, omega, tt):
    """ 
    Return the epoch of periastron from other orbital parameters
        P: orbital period
        ecc: eccentricity
        omega: argument of periastron
        tt: transit epoch
     """
    omega = np.deg2rad(omega)
    E0 = np.arctan2(np.sqrt(1. - ecc**2)*np.cos(omega), np.sin(omega) + ecc)
    Tp = tt - P/(2.*np.pi)*(E0 - ecc * np.sin(E0))
    return Tp


def get_tt(P, ecc, omega, tp):
    """
    Return the transit epoch from other orbital parameters
        P: orbital period
        ecc: eccentricity
        omega: argument of periastron
        tp: epoch of periastron
     """
    E0 = np.arctan2(np.sqrt(1. - ecc**2)*np.cos(omega), np.sin(omega) + ecc)
    Tt = tp + P/(2.*np.pi)*(E0 - ecc*np.sin(E0))
    return Tt




class Planet(object):
    def __init__(self):
        # transit parameters
        self.rp = 0. # planet radius (in units of stellar radii)
        self.a = 0. # semi-major axis (in units of stellar radii)
        self.inc = 0. # orbital inclination (in degrees)
        self.t0 = 0. # time of inferior conjunction
        self.jitter_lc = 0.
        
        # these are fixed
        self.limb_dark = "quadratic" # limb darkening model
        self.u = [.4983, 0.2042]

        # rv parameters
        self.rvsys = 0.
        self.K = 0.
        self.T0 = 0. #function of self.t0  # time of periastron
        self.jitter_rv = 0.

        # parameters shared between transit and rv
        self.period = 0. # orbital period
        self.ecc = 0. # eccentricity
        self.w = 0. # longitude of periastron (in degrees)


        self.N_free_parameters = 11


        ## build the batman model
        self.params = batman.TransitParams()
        self.params.limb_dark = self.limb_dark        #limb darkening model
        self.params.u = self.u     #limb darkening coefficients

        self.labels = [r'$R_p / R_*$',
                       r'$a/R_*$',
                       r'$i$',
                       r'$t_0$',
                       r'$s_{\rm lc}$',
                       'vsys',
                       r'$K$',
                       r'$s_{\rm RV}$',
                       r'$P$',
                       r'$e$',
                       r'$\omega$']
        if fix_ecc:
            self.labels.pop(9)

    # def get_rv_parameters(self):
        # pass
    def set_rv_parameters(self, pars):
        """ pars should be [rvsys, K, jitter_rv] """
        assert len(pars) == 3
        self.rvsys, self.K, self.jitter_rv = pars

    # def get_transit_parameters(self):
        # pass
    def set_transit_parameters(self, pars):
        """ pars should be [rp, a, inc, t0, jitter_lc] """
        assert len(pars) == 5
        self.rp, self.a, self.inc, self.t0, self.jitter_lc = pars

    def set_shared_parameters(self, pars):
        """ pars should be [period, ecc, w] """
        assert len(pars) in (2, 3)
        if len(pars) == 2:
            self.period, self.w = pars
        else:
            self.period, self.ecc, self.w = pars

    def get_rv_curve(self, time, debug=False):
        if debug: print self.w, self.t0
        self.T0 = self.t0 #get_tp(self.period, self.ecc, self.w, self.t0)
        return rv_curve(time,
                        self.rvsys, self.K, self.w, 
                        self.ecc, self.T0, self.period)

    def get_transit_curve(self, time, debug=False):
        self.params.t0 = self.t0
        self.params.per = self.period
        self.params.rp = self.rp
        self.params.a = self.a
        self.params.inc = self.inc
        self.params.ecc = self.ecc
        self.params.w = self.w

        if debug: print self.params.__dict__
        self.batman_model = batman.TransitModel(self.params, time)
        light_curve = self.batman_model.light_curve(self.params)
        return light_curve


    def get_priors(self):
        
        self.prior_rp = stats.uniform(0.1, 0.5) # planet radius (in units of stellar radii)
        self.prior_a = stats.uniform(3, 11) # semi-major axis (in units of stellar radii)
        self.prior_inc = stats.uniform(70, 20) # orbital inclination (in degrees)
        self.prior_t0 = stats.norm(3219.0128, 0.001) # time of inferior conjunction
        self.prior_jitter_lc = stats.reciprocal(0.001, 0.1)

        # rv parameters
        self.prior_rvsys = stats.uniform(32, 2)
        self.prior_K = stats.uniform(50, 500)
        # prior_T0 = get_tp(prior_period, prior_ecc, prior_w, prior_t0)
        self.prior_jitter_rv = stats.reciprocal(1, 50)
 
        # parameters shared between transit and rv
        self.prior_period = stats.uniform(3.2, 0.1) # orbital period
        self.prior_ecc = stats.uniform(0, 0) # eccentricity
        self.prior_w = stats.uniform(0, 360) # longitude of periastron (in degrees)

        return [self.prior_rp,self.prior_a, self.prior_inc, self.prior_t0, self.prior_jitter_lc,
                self.prior_rvsys, self.prior_K, self.prior_jitter_rv,
                self.prior_period, self.prior_ecc, self.prior_w]


    def get_from_prior(self, nwalkers):

        self.get_priors()

        pars_from_prior = []
        for i in range(nwalkers):
            random_rp = self.prior_rp.rvs() # planet radius (in units of stellar radii)
            random_a = self.prior_a.rvs() # semi-major axis (in units of stellar radii)
            random_inc = self.prior_inc.rvs() # orbital inclination (in degrees)
            random_t0 = self.prior_t0.rvs() # time of inferior conjunction
            random_jitter_lc = self.prior_jitter_lc.rvs()
            
            # rv parameters
            random_rvsys = self.prior_rvsys.rvs()
            random_K = self.prior_K.rvs()
            # random_T0 = get_tp(random_period, random_ecc, random_w, random_t0)
            random_jitter_rv = self.prior_jitter_rv.rvs()

            # parameters shared between transit and rv
            random_period = self.prior_period.rvs() # orbital period
            random_ecc = self.prior_ecc.rvs() # eccentricity
            random_w = self.prior_w.rvs() # longitude of periastron (in degrees)

            pars_from_prior.append([random_rp,random_a, random_inc, random_t0, random_jitter_lc,
                                    random_rvsys, random_K, random_jitter_rv,
                                    random_period, random_ecc, random_w])

        return pars_from_prior
        # [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

class Data(object):
    def __init__(self, rv_file, lc_file, skip_rv_rows=0, skip_lc_rows=0):
        self.rv_file = rv_file
        self.lc_file = lc_file

        # read RVs
        self.RVtime, self.RV, self.RVerror = np.loadtxt(rv_file, 
                                                        unpack=True, skiprows=skip_rv_rows)

        # read light curve
        self.LCtime, self.LC, self.LCerror = np.loadtxt(lc_file,
                                                        unpack=True, skiprows=skip_lc_rows)


        # cut last part
        inds = self.LCtime < 57112
        self.LCtime = self.LCtime[inds]
        self.LC = self.LC[inds]
        self.LCerror = self.LCerror[inds]

        self.N_rvs = self.RVtime.size
        self.N_lc = self.LCtime.size



def lnlike(pars, planet, data, debug=False):
    """ pars should be
    # [rp, a, inc, t0, jitter_lc,
    #  rvsys, K, jitter_rv,
    #  period, ecc, w]
    """
    log2pi = np.log(2*np.pi)


    # set the transit params
    planet.set_transit_parameters(pars[:5])
    # set the RV params
    planet.set_rv_parameters(pars[5:8])
    # set the shared params
    planet.set_shared_parameters(pars[8:11])

    # calculate the lnlike for transit
    transit_model = planet.get_transit_curve(data.LCtime)

    sigma = data.LCerror**2 + planet.jitter_lc**2
    chi2 = np.log(sigma)/2. + (data.LC - transit_model)**2 / (2. * (sigma))
    log_like_transit = - data.N_lc/2.0 * np.log(2*np.pi) - np.sum(chi2)

    if debug: print 'log_like_transit', log_like_transit

    # calculate the lnlike for RVs
    rv_model = planet.get_rv_curve(data.RVtime)
    sigma = data.RVerror**2 + planet.jitter_rv**2
    chi2 = np.log(sigma)/2. + (data.RV - rv_model)**2 / (2. * (sigma))
    log_like_rv = - data.N_rvs/2.0 * np.log(2*np.pi) - np.sum(chi2)

    if debug: print 'log_like_rv', log_like_rv

    # sum lnlikes
    log_like = log_like_transit + log_like_rv
    if debug: print log_like

    if not np.isfinite(log_like):
        return -np.inf
    else:
        return log_like


def lnprior(pars, planet, data, debug=False, fix_ecc=fix_ecc):
    """ pars should be
    # [rp, a, inc, t0, jitter_lc,
    #  rvsys, K, jitter_rv,
    #  period, ecc, w]
    """

    # transit parameters
    prior_rp = planet.prior_rp.logpdf(pars[0]) # planet radius (in units of stellar radii)
    prior_a = planet.prior_a.logpdf(pars[1]) # semi-major axis (in units of stellar radii)
    prior_inc = planet.prior_inc.logpdf(pars[2]) # orbital inclination (in degrees)
    prior_t0 = planet.prior_t0.logpdf(pars[3]) # time of inferior conjunction
    prior_jitter_lc = planet.prior_jitter_lc.logpdf(pars[4])
    
    # rv parameters
    prior_rvsys = planet.prior_rvsys.logpdf(pars[5])
    prior_K = planet.prior_K.logpdf(pars[6])
    # prior_T0 = 0.
    prior_jitter_rv = planet.prior_jitter_rv.logpdf(pars[7])

    # parameters shared between transit and rv
    prior_period = planet.prior_period.logpdf(pars[8]) # orbital period
    if fix_ecc:
        prior_ecc = 0.
    else:
        prior_ecc = planet.prior_ecc.logpdf(pars[9]) # eccentricity
    prior_w = planet.prior_w.logpdf(pars[10]) # longitude of periastron (in degrees)

    if debug:
        print prior_rp, '\n', prior_a, '\n', prior_inc, '\n', prior_t0, '\n', prior_jitter_lc \
               , '\n', prior_rvsys, '\n', prior_K, '\n', prior_jitter_rv \
               , '\n', prior_period, '\n', prior_ecc, '\n', prior_w

    ln_prior = prior_rp + prior_a + prior_inc + prior_t0 + prior_jitter_lc \
                + prior_rvsys + prior_K + prior_jitter_rv \
                + prior_period + prior_ecc + prior_w

    if debug: print ln_prior

    return ln_prior

def lnprob(pars, planet, data, debug=True):
    lp = lnprior(pars, planet, data)
    ll = lnlike(pars, planet, data)
    return lp + ll





data = Data(rv_file='EPIC-9792_SOPHIE.rdb', skip_rv_rows=2, 
            lc_file='ktwo211089792c04_FILlpd_LC.txt')
planet = Planet()

ndim, nwalkers = planet.N_free_parameters, 30

pos = planet.get_from_prior(nwalkers)
# pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

import emcee
import corner

run_MCMC = True

if run_MCMC:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(planet, data))
    sampler.run_mcmc(pos, 100)

    samples = sampler.chain[:, 0:, :].reshape((-1, ndim))

    if fix_ecc:
        samples = np.delete(samples, 9, 1)

fig = corner.corner(samples, labels=planet.labels)
fig.savefig('samples.png')
plt.show()

print ['%7s' % s.replace('$', '').replace('_', '').replace('\\rm', '').replace('\\','') for s in planet.labels]
print ['%7.2f' % s for s in np.median(samples, axis=0)]



median_pars = np.median(samples, axis=0)

# set the transit params
planet.set_transit_parameters(median_pars[:5])
# set the RV params
planet.set_rv_parameters(median_pars[5:8])
# set the shared params
planet.set_shared_parameters(median_pars[8:11])


fig = plt.figure()
ax = fig.add_subplot(211)

time = np.linspace(data.LCtime.min(), data.LCtime.max(), 1000)
lc = planet.get_transit_curve(time)
ax.plot(time, lc)
ax.plot(data.LCtime, data.LC)


ax = fig.add_subplot(212)

time = np.linspace(data.RVtime.min(), data.RVtime.max(), 1000)
rv = planet.get_rv_curve(time)
ax.errorbar(data.RVtime, data.RV - median_pars[5], data.RVerror, fmt='o')
ax.plot(time, rv*1e-3)

plt.show()