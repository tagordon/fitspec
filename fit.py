import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d
from astropy.stats import sigma_clip
from wl_utils import run as run_wl
from spec_utils import run as run_spec
import pca

import warnings
warnings.filterwarnings("ignore")

pl_param_names = [
    'period', 'radius', 't0', 
    'semimajor_axis', 'inclination', 
    'eccentricity', 'periastron'
]
st_param_names = ['mh', 'logg', 'teff']

def validate_priors(priors_dict):

    return all(item in priors_dict.keys() for item in pl_param_names)

def validate_st_params(stellar_params_dict):
    
    return  all(item in stellar_params_dict.keys() for item in st_param_names)

def validate_options(wl_options_dir, options_dir):

    return all(wl_options_dir[key] == options_dir[key] for key in ['gp', 'detector', 'n_pca_components'])

def fit_wlc(
    time, 
    spec, 
    wavs,
    priors,
    stellar_params,
    detector,
    polyorder=1,
    cube=None, 
    out_dir='./', 
    samples=30000, 
    burnin=10000, 
    thin=10,
    nproc=1, 
    gp=False,
    out_filter_width=50,
    out_sigma=4,
    n_pca_components=0,
    detrending_vectors=None,
    save_chains=True,
    return_chains=False
):

    if not validate_priors(priors):
        
        raise ValueError(
            '''
                priors Dictionary must contain the following keys: 
                period, radius, t0, semimajor_axis, inclination, 
                eccentricity, periastron.'''
            ) 
        
    if not validate_st_params(stellar_params):
        
        raise ValueError(
            '''
                stellar_params Dictionary must contain the following keys: 
                mh, logg, teff.'''
            )

    time -= time[0]
    time = np.array(time, dtype=np.float64)

    if len(spec.shape) == 2:
        flux = np.sum(spec, axis=1)
    elif len(spec.shape) == 1:
        flux = spec
    else:
        raise ValueError('Flux or spectrum has wrong shape.')

    if len(time) != len(flux):
        raise ValueError('Time array should have same length as flux or spectrum.')
        
    out_mask = sigma_clip(flux - gaussian_filter1d(flux, out_filter_width), sigma=out_sigma).mask

    if (cube is not None) & (n_pca_components > 0):

        if cube.shape[0] != len(time):
            raise ValueError(
                'First dimension of data cube should be the same length as time array.'
            )
        
        pca_components = pca.get_pca_components(cube, n_components=n_pca_components)
        pca_components = pca.replace_outliers_all_components(pca_components, 3, radius=5)
        pca_components = pca.mask_pca(pca_components, ~out_mask)
        pca_components = gaussian_filter1d(pca_components, 5, axis=0)

        if detrending_vectors is not None: 
            detrending_vectors = pca.mask_pca(detrending_vectors, ~out_mask)
            detrending_vectors = np.hstack([detrending_vectors, pca_components])
        else:
            detrending_vectors = pca_components

    elif detrending_vectors is not None:
        detrending_vectors = pca.mask_pca(detrending_vectors, ~out_mask)
        pca_components = None
    else:
        pca_components = None

    sampler = run_wl(
        time[~out_mask],
        flux[~out_mask], 
        detector, 
        priors,
        stellar_params,
        detrending_vectors=detrending_vectors, 
        polyorder=polyorder, 
        samples=samples,
        progress=True,
        nproc=nproc,
        gp=gp
    )

    chains = sampler.get_chain()[burnin::thin, :, :]

    result_dir = {
        'pca_components': pca_components,
        'n_pca_components': n_pca_components,
        'gp': gp,
        'detector': detector,
        'chains': [chains if return_chains else None][0],
        'wl_params': np.median(chains, axis=(0, 1)),
        'mask': out_mask,
        'time': time,
        'spec': spec,
        'wavs': wavs,
        'stellar_params': stellar_params,
        'cube': cube,
        'detrending_vectors': detrending_vectors,
        'polyorder': polyorder
    }

    if save_chains:
        np.save(out_dir + 'wl_mcmc_chains', sampler.get_chain())

    return result_dir

def fit_spec(
    wl_results,
    wav_per_bin=0.02,
    polyorder=None, 
    out_dir='./', 
    samples=10000, 
    burnin=5000, 
    nproc=1, 
    out_filter_width=50,
    out_sigma=4,
    n_pca_components_spec=0,
    n_pca_components=None,
    detrending_vectors=None,
    save_chains=True,
    return_chains=False,
    progress=False,
    gp=None
):

    if polyorder is None:
        polyorder = wl_results['polyorder']

    if (not wl_results['gp']) & gp:
        raise AttributeError(
            '''
            Cannot use GP model for spectral 
            fitting if a GP was not used for 
            white light fitting.
            '''
        )
    if gp is None:
        gp = wl_results['gp']

    time = wl_results['time']
    spec = wl_results['spec']
    stellar_params = wl_results['stellar_params']
    cube = wl_results['cube']
    detector = wl_results['detector']
    wl_params = wl_results['wl_params']
    wavs = wl_results['wavs']

    if n_pca_components is not None:
        n_pca_components = len(wl_results['pca_components'].T)
    else:
        n_pca_components = 0

    time -= time[0]
    time = np.array(time, dtype=np.float64)

    if detector == 'nrs1':

        spec = spec[:, wavs > 2.87]
        wavs = wavs[wavs > 2.87]

    if detector == 'nrs2':
        
        spec = spec[:, wavs < 5.17692]
        wavs = wavs[wavs < 5.17692]

    nbands = np.int64((wavs[-1] - wavs[0]) // wav_per_bin)
    wav_bin_edges = np.linspace(wavs[0], wavs[-1], nbands + 1)
    binned_wavs = wav_bin_edges[:-1] + 0.5 * np.diff(wav_bin_edges)
    
    binned_spec = np.array([
        np.sum(
            spec[:, np.where(
                (wavs >= wav_bin_edges[i]) & 
                (wavs <= wav_bin_edges[i+1])
            )[0]],
            axis=1
        )
        for i in range(nbands)
    ])
    binned_spec = binned_spec.T

    filt = gaussian_filter1d(binned_spec, out_filter_width, axis=0)
    mask = sigma_clip(binned_spec - filt, sigma=out_sigma).mask

    if (cube is not None) & (n_pca_components_spec > 0):

        detrending_cube = np.zeros(
            (binned_spec.shape[1], binned_spec.shape[0], n_pca_components_spec)
        )
    
        for i in range(binned_spec.shape[1]):
            detrending_cube[i, :, :] = pca.get_pca_cube(
                wavs, wav_bin_edges, cube, i, n_components=n_pca_components_spec
            )

        if n_pca_components > 0:
            pca_components = pca.get_pca_components(cube, n_components=n_pca_components)
            pca_components = pca.replace_outliers_all_components(pca_components, 3, radius=5)
            pca_components = gaussian_filter1d(pca_components, 5, axis=0)

            if detrending_vectors is not None: 
                detrending_vectors = np.hstack([detrending_vectors, pca_components])
            else:
                detrending_vectors = pca_components

    else:
        detrending_cube = None

    post = run_spec(
        time,
        binned_spec, 
        wl_params,
        stellar_params,
        wav_bin_edges,
        out_mask=mask,
        detrending_cube=detrending_cube,
        detrending_vectors=detrending_vectors,
        nproc=nproc,
        samples=samples,
        polyorder=polyorder,
        progress=progress,
        gp=gp
    )

    if save_chains:
        for i in range(len(post)):
            np.save(out_dir + 'spec_mcmc_chains_{0}'.format(i), post[i].get_chain()[burnin:, :, :])
            np.save(out_dir + 'wavs_{0}'.format(i), binned_wavs)
    else:
        return post, binned_wavs

    