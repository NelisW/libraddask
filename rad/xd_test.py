__author__ = 'DGriffith'

from libraddask.rad import xd
import xarray
import matplotlib.pyplot as plt

spec_trans = xarray.DataArray([ 0.0, 0.1, 0.3, 0.5, 0.7, 1.0, 1.0, 0.8, 0.0],
                   [('wvl', [550., 600, 650, 700, 750, 800, 850, 950, 1000])],
                   name='trn',
                   attrs={'trn_units': '1', 'wvl_units': 'nm', 'extrap_val': 0.0})
spec_transB = xarray.DataArray([ 0.0, 0.5, 0.3, 0.25, 0.4, 0.45, 0.6, 0.7],
                   [('wvl', [551., 600, 660, 715, 755, 845, 851, 956])],
                   name='trn',
                   attrs={'trn_units': '', 'wvl_units': 'nm', 'extrap_val': 0.0})

print(spec_trans)
print('+++++++++++++++++++++++++++++++++++')
print(spec_transB)
print('+++++++++++++++++++++++++++++++++++')
xd_harmonise_interp((spec_trans, spec_transB))
