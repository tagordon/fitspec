from distributions import *

stellar_params_dict = {

    
    '134.01' : {
        'mh': 0.04,
        'logg': 4.04,
        'teff': 3800,
        'reference': 'https://ui.adsabs.harvard.edu/abs/2020A%26A...636A..58A/abstract'
    },

    '175.01': {
        'mh': -0.46,
        'logg': 4.86,
        'teff': 3415,
        'reference': 'https://ui.adsabs.harvard.edu/abs/2021A%26A...653A..41D/abstract'  
    },

    '260.01': {
        'mh': -0.211,
        'logg': 4.7,
        'teff': 4050,
        'reference': 'https://ui.adsabs.harvard.edu/abs/2023AJ....166...33M/abstract'  
    },

    '402.01': {
        'mh': 0.03,
        'logg': 4.37,
        'teff': 5131,
        'reference': 'https://ui.adsabs.harvard.edu/abs/2019A&A...627A..43D/abstract'
    },

    '402.02': {
        'mh': 0.03,
        'logg': 4.37,
        'teff': 5131,
        'reference': 'https://ui.adsabs.harvard.edu/abs/2019A&A...627A..43D/abstract'
    },

    '455.01': {
        
    },

    '562.01': {
        'mh': -0.12,
        'logg': 4.94,
        'teff': 3505,
        'reference': 'https://ui.adsabs.harvard.edu/abs/2023AJ....165..134O/abstract'
    },

    '776.01': {
        'mh': -0.21,
        'logg': 4.8,
        'teff': 3725,
        'reference': 'https://ui.adsabs.harvard.edu/abs/2024A&A...684A..12F/abstract'
    },

    '776.02': {
        'mh': -0.21,
        'logg': 4.8,
        'teff': 3725,
        'reference': 'https://ui.adsabs.harvard.edu/abs/2024A&A...684A..12F/abstract'
    },

    '836.01': {
        'mh': -0.284,
        'logg': 4.743,
        'teff': 4552
    },

    '836.02': {
        'mh': -0.284,
        'logg': 4.743,
        'teff': 4552
    },

}

priors_dict = {

    '134.01': {
        'period': normal_prior(1.4015, 0.00018),
        'radius': normal_prior(0.0212, 0.001),
        't0': uniform_prior(0.05, 0.125, init=0.085),
        'semimajor_axis': normal_prior(7.61, 0.31),
        'inclination': normal_prior(85.8 * np.pi / 180, 0.8 * np.pi / 180),
        'eccentricity': uniform_prior(0.0, 0.21, init=0.01),
        'periastron': uniform_prior(0, np.pi, init=np.pi/2),
        'reference': 'https://ui.adsabs.harvard.edu/abs/2020A%26A...636A..58A/abstract'
    },

    '175.01': {
        'period': normal_prior(3.6906777, 0.0000026),
        'radius': normal_prior(0.04088, 0.00068),
        't0': uniform_prior(0.07, 0.1, init=0.08),
        'semimajor_axis': normal_prior(19, 1.2),
        'inclination': normal_prior(88.11 * np.pi / 180, 0.36 * np.pi / 180),
        'eccentricity': trunc_normal_prior(0.103, 0.058, 0, 1),
        'periastron': uniform_prior(0.0, np.pi, init=np.pi/2),
        'reference': 'https://ui.adsabs.harvard.edu/abs/2021A%26A...653A..41D/abstract'
    },

    '260.01': {
        'period': normal_prior(13.475815, 0.000047),
        'radius': normal_prior(0.02469, 0.002),
        't0': uniform_prior(0.1, 0.15),
        'semimajor_axis': normal_prior(35.32, 0.9), # ExoFOP
        'inclination': normal_prior(np.pi / 2, 2.0 * np.pi / 180, init=88.7 * np.pi / 180),
        'eccentricity': uniform_prior(0, 1, init=0.01),
        'periastron': uniform_prior(0.0, np.pi, init=np.pi/2),
        'reference': 'https://ui.adsabs.harvard.edu/abs/2024ApJS..272...32P/abstract'
    },

    '402.01': {
        'period': normal_prior(4.75615, 0.00017),
        'radius': normal_prior(0.01761, 0.0006),
        't0': uniform_prior(0.15, 0.2),
        'semimajor_axis': normal_prior(13.11, 0.16), 
        'inclination': normal_prior(88.5 * np.pi / 180, 0.6 * np.pi / 180),
        'eccentricity': trunc_normal_prior(0.09, 0.05, 0, 1),
        'periastron': uniform_prior(0.0, np.pi, init=np.pi/2),
        'reference': 'https://ui.adsabs.harvard.edu/abs/2019ApJ...876L..24G/abstract'
    },

    '402.02': {
        'period': normal_prior(17.1784, 0.00016),
        #'radius': normal_prior(0.0256, 0.0011),
        'radius': normal_prior(0.0256, 0.01), # inflate uncertainty on prior 
        't0': uniform_prior(0.0, 0.15),
        #'semimajor_axis': normal_prior(31.87, 0.7), 
        'semimajor_axis': normal_prior(31.87, 5.0), # inflate uncertainty here also 
        'inclination': trunc_normal_prior(88.5 * np.pi / 180, 1.0 * np.pi / 180, 88.2 * np.pi / 180, np.pi / 2, init=88.35 * np.pi / 180),
        'eccentricity': trunc_normal_prior(0.05, 0.06, 0, 1),
        'periastron': uniform_prior(0.0, np.pi, init=np.pi/2),
        'reference': 'https://ui.adsabs.harvard.edu/abs/2019ApJ...876L..24G/abstract'
    },

    '455.01': {
        
    },

    '562.01': {
        'period': normal_prior(3.930600, 0.000031),
        'radius': normal_prior(0.0309, 0.0010),
        't0': uniform_prior(0.075, 0.125),
        'semimajor_axis': normal_prior(22.89, 1.211), 
        'inclination': normal_prior(89.228 * np.pi / 180, 0.483 * np.pi / 180),
        'eccentricity': uniform_prior(0, 1, init=0.01),
        'periastron': uniform_prior(0.0, np.pi, init=np.pi/2),
        'reference': 'https://ui.adsabs.harvard.edu/abs/2023AJ....165..134O/abstract'
    },

    '776.02': {
        #'period': normal_prior(8.246620, 0.000002),
        'period': normal_prior(8.246620, 0.5),
        #'radius': normal_prior(0.0316, 0.0008),
        'radius': normal_prior(0.03, 0.01), # trying this because other prior seemed bad. 
        #'t0': uniform_prior(0.1, 0.2, init=0.17),
        't0': uniform_prior(0.1, 0.5, init=0.17),
        #'semimajor_axis': normal_prior(27.87, 1.02), 
        'semimajor_axis': normal_prior(27.87, 10.0), # seemed too restrictive
        #'inclination': normal_prior(89.65 * np.pi / 180, 0.37 * np.pi / 180),
        'inclination': normal_prior(89.65 * np.pi / 180, 5.0 * np.pi / 180),
        #'eccentricity': trunc_normal_prior(0.052, 0.037, 0, 1),
        'eccentricity': trunc_normal_prior(0.052, 0.3, 0, 1),
        #'periastron': uniform_prior(0.0, np.pi, init=np.pi/2),
        'periastron': normal_prior(45 * np.pi / 180, 1.0 * np.pi / 180),
        'reference': [
            'https://ui.adsabs.harvard.edu/abs/2024A&A...684A..12F/abstract',
            'https://ui.adsabs.harvard.edu/abs/2023ApJS..265....4K/abstract'
        ]
    },

    '776.01': {
        'period': normal_prior(15.665323, 0.000075),
        'radius': normal_prior(0.0344, 0.0009),
        't0': uniform_prior(0.1, 0.2, init=0.125),
        'semimajor_axis': normal_prior(39.4, 1.8), 
        'inclination': normal_prior(89.51 * np.pi / 180, 0.25 * np.pi / 180),
        'eccentricity': trunc_normal_prior(0.0890, 0.0540, 0, 1),
        'periastron': uniform_prior(0.0, np.pi, init=np.pi/2),
        'reference': [
            'https://ui.adsabs.harvard.edu/abs/2024A&A...684A..12F/abstract',
            'https://ui.adsabs.harvard.edu/abs/2023ApJS..265....4K/abstract'
        ]
    },

    '836.02': {
        'period': normal_prior(3.81673, 0.00001),
        #'radius': normal_prior(0.0235, 0.0013),
        'radius': normal_prior(0.0235, 0.01), # inflate uncertainty on prior
        't0': uniform_prior(0.05, 0.15),
        'semimajor_axis': normal_prior(13.6456, 5.0), # also calculated from a (au) with guessed error 
        'inclination': normal_prior(87.57 * np.pi / 180, 0.44 * np.pi / 180),
        'eccentricity': trunc_normal_prior(0.053, 0.042, 0, 1),
        'periastron': uniform_prior(0.0, np.pi, init=np.pi/2),
        'reference': 'https://ui.adsabs.harvard.edu/abs/2023MNRAS.tmp..458H/abstract'
    },

    '836.01': {
        'period': normal_prior(8.59545, 0.00001),
        'radius': normal_prior(0.0357, 0.0018),
        't0': uniform_prior(0.15, 0.2),
        'semimajor_axis': normal_prior(24.2517, 5.0), # also calculated from a (au) with guessed error 
        'inclination': normal_prior(88.7 * np.pi / 180, 0.15 * np.pi / 180),
        'eccentricity': trunc_normal_prior(0.0780, 0.056, 0, 1),
        'periastron': uniform_prior(0.0, np.pi, init=np.pi/2),
        'reference': 'https://ui.adsabs.harvard.edu/abs/2023MNRAS.tmp..458H/abstract'
    },
}