#! /usr/bin/env python

"""
Creates a table showing the values of the various telemetry parameters at the
given MJDs. Also searches for telemetry correlations with user-defined
parameters.

The following outputs are created:
1. telemetry_table.txt
3. correlations.txt
2. correlation_matrix.png

Authors
-------
    Ben Sunnquist, 2018

Use
---
    This script can be run via the command line as such:
        
        python wfc3_telemetry.py --t <t>

    --t [Required]: The name of the table containing the MJDs of interest.

    This script can also be run within python:

        >>> import wfc3_telemetry as w
        >>> w.wfc3_telemetry(t)
    
    t [Required]: str
        The name of the table containing the MJDs of interest.

Notes
-----
    Written in Python 3.
"""

import glob
import os

import argparse
from astropy.coordinates import get_sun
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from bisect import bisect
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

# -----------------------------------------------------------------------------

def distcorr(X, Y):
    """ Compute the distance correlation function
    
    >>> a = [1,2,3,4,5]
    >>> b = np.array([1,2,9,4,4])
    >>> distcorr(a, b)
    0.762676242417
    """
    X = np.atleast_1d(X)
    Y = np.atleast_1d(Y)
    if np.prod(X.shape) == len(X):
        X = X[:, None]
    if np.prod(Y.shape) == len(Y):
        Y = Y[:, None]
    X = np.atleast_2d(X)
    Y = np.atleast_2d(Y)
    n = X.shape[0]
    if Y.shape[0] != X.shape[0]:
        raise ValueError('Number of samples must match')
    a = squareform(pdist(X))
    b = squareform(pdist(Y))
    A = a - a.mean(axis=0)[None, :] - a.mean(axis=1)[:, None] + a.mean()
    B = b - b.mean(axis=0)[None, :] - b.mean(axis=1)[:, None] + b.mean()
    
    dcov2_xy = (A * B).sum()/float(n * n)
    dcov2_xx = (A * A).sum()/float(n * n)
    dcov2_yy = (B * B).sum()/float(n * n)
    dcor = np.sqrt(dcov2_xy)/np.sqrt(np.sqrt(dcov2_xx) * np.sqrt(dcov2_yy))
    return dcor

# -----------------------------------------------------------------------------

def dra_hst_sun(mjd):
    """
    Finds the difference in RA between HST and the Sun at a given MJD.

    Parameters
    ----------
    mjd : float
        The MJD of interest.

    Returns
    -------
    dra : float
        The difference in RA between HST and the Sun at the given MJD.
    """
    mjd = Time(mjd, format='mjd')

    # get sun position at this mjd
    sun_pos = get_sun(mjd)
    sun_ra = sun_pos.ra.deg
    sun_dec = sun_pos.dec.deg

    # get HST lat/long info
    hst = fits.getdata('HST_SubEarth.fits', 1)
    hst_mjd = hst.mjd
    hst_lat = hst.sublat
    hst_long =hst.sublng
    
    # get HST position at this mjd (using closest matching HST entry)
    first = hst_mjd[bisect(hst_mjd, mjd.value)]
    last = hst_mjd[bisect(hst_mjd, mjd.value) - 1]
    if abs(first-mjd.value) < abs(last-mjd.value):
        closest_mjd = first
    else:
        closest_mjd = last
    pos = float(hst_long[hst_mjd==closest_mjd])

    # turn position into hst ra: RA = GMST + long
    gst = mjd.sidereal_time('mean', 'greenwich').value
    hst_ra = pos + gst*15
    if hst_ra < 0:
        hst_ra += 360
    if hst_ra > 360:
        hst_ra -= 360

    # delta_ra [HST - Sun]
    dra = hst_ra - sun_ra
    if dra < -180:
        dra += 360
    if dra > 180:
        dra -= 360

    return dra

# -----------------------------------------------------------------------------

def find_correlations(t):
    """
    Searches for telemetry correlations with the parameter of interest.

    Parameters
    ----------
    t : str
        The table containing columns MJD and the values of the parameter of
        interest and the various telemetry parameters at those MJDs.

    Outputs
    -------
    correlations.txt
        A file summarizing the telemetry correlations.

    correlation_matrix.png
        A correlation matrix of the telemetry parameters that correlate most
        with the parameter of interest.
    """

    # Read in the data
    t = pd.read_csv(t)

    # Record the linear and non-linear correlations with the parameter
    cor_linear = []
    cor_dist = []
    for p in t.columns:
        # Find the correlation coefficients
        try:
            c = t[p].corr(t['param'])
            cor_linear.append(c.round(3))
            d = distcorr(t['param'], t[p])
            cor_dist.append(d.round(3))
        except TypeError:
            # transform unique strings to unique ints
            u, indices = np.unique(t[p], return_inverse=True)
            df = pd.DataFrame({'indices': indices}) 
            c = df['indices'].corr(t['param'])
            cor_linear.append(c.round(3))
            d = distcorr(t['param'], df['indices'])
            cor_dist.append(d.round(3))

    # Write out a summary table of the correlations
    t_out = pd.DataFrame({})
    t_out['param'] = t.columns
    t_out['R'] = cor_linear
    t_out['d'] = cor_dist
    t_out = t_out.fillna(0) # replace nans with zero
    t_out['R_abs'] = abs(t_out['R']) # sort by magnitude of correlation
    t_out = t_out.sort_values('R_abs', ascending=False).reset_index(drop=True)
    t_out = t_out.drop('R_abs', axis=1)
    t2_out = Table.from_pandas(t_out)
    t2_out.write('correlations.txt', format='ascii.fixed_width_two_line', 
                 overwrite=True)

    # Calculate the correlation matrix using the top 40 correlations with param
    names = list(t_out[0:41]['param'].values)
    t_top  = t[names].copy()
    correlations = t_top.corr()
    correlations = correlations.fillna(0) # replace nans with zero

    # Plot the correlation matrix
    fig = plt.figure()
    ax = fig.add_axes([0,0,2,2])
    cax = ax.matshow(correlations, vmin=-1, vmax=1, cmap='coolwarm')
    fig.colorbar(cax, label='Pearson Correlation Coefficient', pad=.04)
    ticks = np.arange(0, len(names), 1)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(names, rotation=90)
    ax.set_yticklabels(names)
    ax.set_xlim(-0.5, len(names)-1.5)
    ax.set_ylim(len(names)-1.5, -0.5)

    # Save the correlation matrix
    fig.savefig('correlation_matrix.png', bbox_inches='tight', dpi=200,
                overwrite=True)

# -----------------------------------------------------------------------------

def get_telemetry(t):
    """
    Finds the telemetry values at the MJDs of interest.

    Parameters
    ----------
    t : str
        The table containing columns MJD and the values of the parameter of
        interest at those MJDs.

    Outputs
    -------
    telemetry_table.txt
        A file containing all of the telemetry values at the MJDs of interest.
    """

    # Read in the user table containing the MJDs of interest
    t = pd.read_csv(t)

    # Get paths to all HST telemetry parameters
    hst_telemetry_files = glob.glob('/grp/hst/telemetry/I*')

    for i,f in enumerate(hst_telemetry_files):
        # Read in the telemetry parameter data
        param = os.path.basename(f)
        print('Finding matching {} entries'.format(param))

        if not param == 'I542_6month': # this is just a subset of I542
            data = pd.read_csv(f, names=['mjd', 'param'], header=None, 
                               delim_whitespace=True)

            # Some 2nd column parameters have spaces, i.e. 'Not Init'. To
            # get around this, read in the file as 3 columns.
            try:
                test = data['mjd'][0]
            except KeyError:
                data = pd.read_csv(f, names=['mjd', 'param', 'param_extra'], 
                                   header=None, delim_whitespace=True)

            # Find the matching parameter value for each MJD of interest
            def match(mjd):
                return data['param'][bisect(data['mjd'], mjd) - 1]

            matching_param = [match(mjd) for mjd in t['mjd']]

            # Add these values to the table
            t[param] = matching_param
        
        # tracker
        if i%10==0:
            print(str(i) + '/' + str(len(hst_telemetry_files)) + ' complete.')

    # Add difference in RA between HST and the Sun at MJDs to the table
    print('Finding difference in RA between HST and Sun at MJDs.')
    try:
        dras = [dra_hst_sun(mjd) for mjd in t['mjd']]
        t['dra_hst_sun'] = dras
    except IndexError:
        pass

    # Write out the final telemetry table
    t.to_csv('telemetry_table.txt', index=False)

# -----------------------------------------------------------------------------

def parse_args():
    """
    Parses command line arguments.
    
    Returns
    -------
    args : object
        Contains the table name that contains the MJD/parameter values.
    """

    t_help = ('The table containing columns MJD and the values of the '
              'parameter of interest at those MJDs.')

    parser = argparse.ArgumentParser()
    parser.add_argument('--t', dest='t', action='store', type=str, 
                        required=True, help=t_help)
    args = parser.parse_args()

    return args

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def wfc3_telemetry(t):
    """
    The main function of the wfc3_telemetry module. See module docstring for
    further details.

    Parameters
    ----------
    t : str
        The table containing columns MJD and the values of the parameter of
        interest at those MJDs.
    """

    # Create a table of the values of all of the telemetry parameters at the
    # MJDs of interest
    get_telemetry(t)
    print('Done finding the telemetry values at the given MJDs.')

    # Search for correlations between the telemetry parameters and the user-
    # defined parameter of interest
    find_correlations('telemetry_table.txt')
    print('Finished searching for telemetry correlations.')

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    
    args = parse_args()

    wfc3_telemetry(args.t)
