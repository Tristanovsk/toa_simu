import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

import chart_studio
chart_studio.tools.set_credentials_file(username='tristanovsk', api_key='e2aYFdkLAwuAnCZn8aWe')
import chart_studio.plotly as py
import plotly.express as px
import plotly.offline as po

import coxmunk.coxmunk as coxmunk
import RTxploitation
# from grs import acutils
# import grs.auxdata as grsdata
import toa_simu.auxdata as ad

from Py6S import *

publish= False
sza = 30
vza = 10
azi = 10
wind = 2
wind_azi = 0
stats = 'bh2006'
dirlut = os.path.abspath('/DATA/projet/VRTC/lut/norm_rad_toa')

wl = np.linspace(400, 2400, 200)
N = len(wl)

####################################
#    Sunglint reflectance
####################################
wri = ad.water_refractive_index()
r_g = [coxmunk.sunglint(sza, vza, azi, m=m).sunglint(wind, wind_azi, stats=stats) for m in wri.n_harmel(wl)]

####################################
#    Rayleigh and aerosol opt. thickness
####################################
rot = ad.rot().get_data(wl)

####################################
#     6S absorbing gases transmittance
####################################

s = SixS()
s.geometry.solar_z = sza
s.geometry.solar_a = 0
s.geometry.view_z = vza
s.geometry.view_a = azi
s.aero_profile = AeroProfile.PredefinedType(AeroProfile.Maritime)
wavelengths, results = SixSHelpers.Wavelengths.run_wavelengths(s, wl / 1000)
F0, trans_gas, irradiance = [], [], []
for i in range(N):
    res = results[i]
    F0 = np.append(F0, res.solar_spectrum)
    trans_gas = np.append(trans_gas, res.total_gaseous_transmittance)
    irradiance = np.append(irradiance,
                           res.diffuse_solar_irradiance + res.direct_solar_irradiance)

Es_toa = F0 * np.cos(np.radians(sza))

####################################
#   Open LUT file
####################################

models = ['rg0.10_sig0.46_nr1.51_ni-0.0200', 'rg0.10_sig0.46_nr1.45_ni-0.0001',
          'rg0.80_sig0.60_nr1.51_ni-0.0200', 'rg0.80_sig0.60_nr1.35_ni-0.0010']
ichoice=2
tiltes= ['TOA radiance', 'TOA normalized radiance','TOA normalized radiance corrected for gaseous absorption']
params= ['L_TOA', 'L_NTOA','L_NTOA_gascorr']
norms = [trans_gas * Es_toa, trans_gas, 1]

norm = norms[ichoice]
df_ = pd.DataFrame()

for model in models:

    lutfile = os.path.join(dirlut, 'lut_toa_rad_aero_' + model + '.nc')

    lut = xr.open_dataarray(lutfile, group='I')
    aot = xr.open_dataarray(lutfile, group='aot')
    Ln_toa_ = lut.isel(z=1).interp(sza=sza, vza=vza, azi=azi)
    Ln_toa = Ln_toa_.interp(wl=wl / 1000, method='quadratic')  # 'cubic')
    # Ln_toa = Ln_toa_.interp( wl=wl / 1000, method='nearest')
    aot_ = aot.interp(wl=wl / 1000, method='quadratic')

    for aot550 in [0.01, 0.1, 0.2, 0.3, 0.5, 0.8]:
        trans_atmo = np.exp(
            -(rot + aot_.interp(aot=aot550)) * (1 / np.cos(np.radians(sza)) + 1 / np.cos(np.radians(vza))))
        Ltoa = norm * Ln_toa.interp(aot=aot550)
        Ltoa_g = Ltoa + norm * trans_atmo * r_g
        df = pd.DataFrame({'wl': wl, 'Ltoa': Ltoa, 'Ltoa_g': Ltoa_g})
        df['aot550'] = aot550
        df['model'] = model
        df_ = pd.concat([df_, df])

####################################
#  Plot spectra
####################################
df__ = df_.melt(id_vars=['wl', 'aot550','model'])
by_ = 'aot550'

fig = px.scatter(df__, x="wl", y="value", color="variable",facet_col='model',facet_col_wrap=2,
                 # hover_name="site", hover_data=['date','satellite','Rrs_count'],
                 opacity=0.5,
                 # range_x=[-0.005,0.04],range_y=[-0.005,0.04],
                 animation_frame=by_,
                 title=tiltes[ichoice]+" with and without sunglint, sza= {:.1f}°, vza= {:.1f}° ".format(sza,vza),height =950)#,width=1200)
po.plot(fig,filename=params[ichoice]+'_test.html')
if publish:
    py.plot(fig,filename=params[ichoice]+'_test.html',auto_open=False)

####################################
#     Set SMAC patrameters for
#    absorbing gases correction
####################################
import lowtran
import lowtran.plots as lp

lowtran.nm2lt7(200, 2500, 20)
c1 = {'model': 6,
      'h1': 0,
      'angle': [0, 30, 60],
      'wlshort': 300,
      'wllong': 2600,
      'wlstep': 5,

      }

TR = lowtran.transmittance(c1)
lp.plottrans(TR, c1)
TR = lowtran.radiance(c1)
lp.plotradiance(TR, c1)
TR = lowtran.irradiance(c1)
lp.plotirrad(TR,c1)



s = SixS()
s.geometry.solar_z = sza
s.geometry.solar_a = 0
s.geometry.view_z = vza
s.geometry.view_a = azi
s.aero_profile = AeroProfile.PredefinedType(AeroProfile.Maritime)
parameter = 'apparent_radiance'
parameter = 'direct_solar_irradiance'
params = ['transmittance_no2.total', 'total_gaseous_transmittance', 'apparent_radiance', 'direct_solar_irradiance',
          'diffuse_solar_irradiance']

ss = []
for wl_ in wl / 1000:
    s.wavelength = Wavelength(wl_)
    s.run()
    ss.append(s.outputs)

wavelengths, values = SixSHelpers.Wavelengths.run_wavelengths(s, wl / 1000, n=8, output_name=params[0])
plt.figure()
SixSHelpers.Wavelengths.plot_wavelengths(wavelengths, values, 'Pixel Radiance (W/m^2)')

sensordata = grsdata.sensordata('S2A')
smac = acutils.smac(sensordata.smac_bands, sensordata.smac_dir)
smac.set_gas_param()
# smac.set_values(o3du=340, h2o=2)
smac.set_standard_values(l2h.pressure_msl)
l2h.aux.no2 = smac.uno2
