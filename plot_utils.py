import numpy as np
import batman
import exotic_ld
from exotic_ld import StellarLimbDarkening
from scipy.ndimage import gaussian_filter1d
from astropy.stats import sigma_clip
import celerite2

from distributions import *
from transit import *
from wl_utils import _unpack_params, get_trend_model

def get_canonical_params(result, polyorder):

    detrending_vectors = result['detrending_vectors']
    polyorder = result['polyorder']
    gp = result['gp']

    if detrending_vectors is None:
        ncoeffs = 1 + polyorder
    else:
        ncoeffs = 1 + polyorder + len(detrending_vectors.T)

    if gp:
        p = result['chains'][:, :, 5 + ncoeffs:]
    else:
        p = result['chains'][:, :, 3 + ncoeffs:]

    return np.array(np.vectorize(canon_from_eastman)(*p.T))

def _get_model_from_sample(time, p, detrending_vectors=None, polyorder=1, gp=False, flux=None):

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
    '''
    Get posterior samples for the the transit model, 
    systematics model, and GP model (if applicable).

    Args:
        wl_results (dict): results dictionary returned from fit_wlc
        nsamples (int, default=100): number of posterior samples 
        to return. 

    Returns: 
        transit models (2D array): the transit model
        systematics models (2D array): polynomial + linear combination 
        detrending vectors and PCA vectors as specified in fit_wlc
        GP prediction (2D array): only returned if gp=True in fit_wlc
        f0 (1D array): the constant term of the polynomial in the systematics 
        model, which is useful for plotting the systematics model and GP 
        prediction.

    Note: the systematics model and GP prediction vectors will have zero flux 
    offset, so f0 should be added to these if they are to be plotted over the 
    observations. 
    '''

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
    polyorder = wl_results['polyorder']

    if wl_results['chains'] is None:
        raise AttributeError(
            '''
            Result dictionary does not contain MCMC chains. 
            Run fit_wlc with return_chains=True.
            '''
        )

    flux = np.sum(spec, axis=1)
    flat_samples = np.concatenate(wl_results['chains'], axis=0)

    ret = _get_models(
        time[~mask], 
        flat_samples, 
        nsamples=100, 
        detrending_vectors=detrending_vectors, 
        gp=gp, 
        polyorder=polyorder,
        flux=flux[~mask], 
        return_idx=True
    )

    idx = ret[-1]
    if gp:
        f = flat_samples[idx, 3]
    else:
        f = flat_samples[idx, 1]

    return *ret[:-1], f

def _get_models(time, flat_samples, polyorder=1, nsamples=100, detrending_vectors=None, gp=False, flux=None, return_idx=False):

    idx = np.random.randint(len(flat_samples), size=nsamples)
    if gp & (flux is not None):
        trans = np.zeros((nsamples, len(time)))
        sys = np.zeros((nsamples, len(time)))
        pred = np.zeros((nsamples, len(time)))
        for i, j in enumerate(idx):
            trans[i], sys[i], pred[i] = _get_model_from_sample(
                time, 
                flat_samples[j], 
                detrending_vectors=detrending_vectors, 
                gp=gp, 
                polyorder=polyorder,
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
            trans[i], sys[i] = _get_model_from_sample(
                time, 
                flat_samples[j], 
                detrending_vectors=detrending_vectors, 
                gp=gp, 
                polyorder=polyorder,
                flux=flux
            )
        if return_idx:
            return trans, sys, idx
        else:
            return trans, sys