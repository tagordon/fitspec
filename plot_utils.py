import numpy as np
import batman
import exotic_ld
from exotic_ld import StellarLimbDarkening
from scipy.ndimage import gaussian_filter1d
from astropy.stats import sigma_clip
import celerite2

from distributions import *
from transit_parameterization import *
from wl_utils import _unpack_params, get_trend_model

def get_model_from_sample(time, p, detrending_vectors=None, polyorder=1, gp=False, flux=None):

    if detrending_vectors is None:
        len_det = 0
    else:
        len_det = len(detrending_vectors.T)

    ncoeffs = 1 + polyorder + len_det

    if gp:

        err, lsigma, lw0, coeffs, f, p = _unpack_params(p, ncoeffs, gp=True)
        trend = get_trend_model(time, detrending_vectors, coeffs[1:], polyorder)
        term = celerite2.terms.SHOTerm(sigma=np.exp(lsigma), w0=np.exp(lw0), Q=1/np.sqrt(2))
        gp = celerite2.GaussianProcess(term)
        mu, jac = reparam(time, p)

        if flux is not None:
            gp.compute(time, diag=err**2)
            pred = gp.predict(flux - mu * f - trend)
            return mu * f, trend, pred
            
        if flux is None:
            return mu * f, trend
        
    else:

        err, coeffs, f, p = _unpack_params(p, ncoeffs, gp=False)
        trend = get_trend_model(time, detrending_vectors, coeffs[1:], polyorder)
        mu, jac = reparam(time, p)
        
    return mu * f, trend

def get_wl_models(wl_results, nsamples=100):

    time = wl_results['time']
    spec = wl_results['spec']
    stellar_params = wl_results['stellar_params']
    cube = wl_results['cube']
    detector = wl_results['detector']
    gp = wl_results['gp']
    wl_params = wl_results['wl_params']
    wavs = wl_results['wavs']
    mask = wl_results['mask']
    detrending_vectors = wl_results['detrending_vectors']

    flux = np.sum(spec, axis=1)
    flat_samples = np.concatenate(wl_results['chains'], axis=0)

    ret = get_models(
        time[~mask], 
        flat_samples, 
        nsamples=100, 
        detrending_vectors=detrending_vectors, 
        gp=gp, 
        flux=flux[~mask], 
        return_idx=True
    )

    idx = ret[-1]
    if gp:
        f = flat_samples[idx, 3]
    else:
        f = flat_samples[idx, 1]

    return *ret[:-1], f

def get_models(time, flat_samples, nsamples=100, detrending_vectors=None, gp=False, flux=None, return_idx=False):

    idx = np.random.randint(len(flat_samples), size=nsamples)
    if gp & (flux is not None):
        trans = np.zeros((nsamples, len(time)))
        sys = np.zeros((nsamples, len(time)))
        pred = np.zeros((nsamples, len(time)))
        for i, j in enumerate(idx):
            trans[i], sys[i], pred[i] = get_model_from_sample(
                time, 
                flat_samples[j], 
                detrending_vectors=detrending_vectors, 
                gp=gp, 
                flux=flux
            )
        if return_idx:
            return trans, sys, pred, idx
        else:
            return trans, sys, pred

    else:
        trans = np.zeros((nsamples, len(time)))
        sys = np.zeros((nsamples, len(time)))
        for i, j in enumerate(idx):
            trans[i], sys[i] = get_model_from_sample(
                time, 
                flat_samples[j], 
                detrending_vectors=detrending_vectors, 
                gp=gp, 
                flux=flux
            )
        if return_idx:
            return trans, sys, idx
        else:
            return trans, sys