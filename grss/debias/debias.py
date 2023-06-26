import os
import numpy as np
import healpy as hp

__all__ = [ 'debias_lowres_path',
            'debias_hires_path',
            'debias_obs',
            'icrf2radec',
            'radec2icrf'
]

cwd = os.path.dirname(os.path.realpath(__file__))
debias_lowres_path = os.path.join(cwd,'lowres_data')
debias_hires_path = os.path.join(cwd,'hires_data')

def debias_obs(right_asc,dec,epoch,catalog,biasdf,nside=256):
    j2000_jd = 2451545.0
    # arcseconds to radians
    as2rad = 1/3600*np.pi/180
    # find pixel from RADEC
    idx = hp.ang2pix(nside, np.rad2deg(right_asc), np.rad2deg(dec), nest=False, lonlat=True)
    # find catalog data in pandas Data Frame
    colnames = [f'{catalog}_ra', f'{catalog}_dec', f'{catalog}_pm_ra', f'{catalog}_pm_dec']
    bias = biasdf[colnames].iloc[idx]
    # time from epoch in Julian years
    dt_jy = (epoch-j2000_jd)/365.25
    # bias correction
    ddec = (bias[colnames[1]]+dt_jy*bias[colnames[3]]/1000)*as2rad
    dec_deb = dec-ddec
    dra = (bias[colnames[0]]+dt_jy*bias[colnames[2]]/1000)*as2rad / np.cos(dec)
    ra_deb = right_asc-dra
    # Quadrant correction
    xyz = radec2icrf(ra_deb, dec_deb, deg=False)
    ra_deb, dec_deb = icrf2radec(xyz[0], xyz[1], xyz[2], deg=False)
    return ra_deb, dec_deb

def icrf2radec(pos_x, pos_y, pos_z, deg=True):
    pos = np.array([pos_x, pos_y, pos_z])
    dist = np.linalg.norm(pos,axis=0) if (pos.ndim>1) else np.linalg.norm(pos)
    phi = np.arctan2(pos_y/dist,pos_x/dist)
    delta = np.arcsin(pos_z/dist)
    if deg:
        right_asc = np.mod(np.rad2deg(phi)+360,360)
        dec = np.rad2deg(delta)
    else:
        right_asc = np.mod(phi+2*np.pi,2*np.pi)
        dec = delta
    return right_asc, dec

def radec2icrf(right_asc, dec, deg=True):
    if deg:
        alpha = np.deg2rad(right_asc)
        delta = np.deg2rad(dec)
    else:
        alpha = np.array(right_asc)
        delta = np.array(dec)
    cosd = np.cos(delta)
    pos_x = cosd*np.cos(alpha)
    pos_y = cosd*np.sin(alpha)
    pos_z = np.sin(delta)
    return np.array([pos_x, pos_y, pos_z])
