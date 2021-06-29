import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline

print('chi2_pro.py called!')
print('function to call: grid_val(ptx,pty,ptz)')
print('use return as: ')
print('X_val=grid_val(X_pts,Y_pts,Z_pts)[0]')
print('Y_val=grid_val(X_pts,Y_pts,Z_pts)[1]')
print('Z_val=grid_val(X_pts,Y_pts,Z_pts)[2]')


def grid_val(ptx,pty,ptz):
    
    my_df = pd.DataFrame.from_dict({'y':pty, 'x':ptx, 'fxy': ptz})
    
    def bivariate_interpolation(df, FXY= 'fxy', Y='y', X='x'): 
        df_copy = df.copy().sort_values(by=[Y, X], ascending=True)
        x = np.array(df_copy[X].unique(), dtype='float64')
        y = np.array(df_copy[Y].unique(), dtype='float64')
        z = np.array(df_copy[FXY].values, dtype='float64')
        Z = z.reshape(len(y), len(x))

        interp_spline = interpolate.RectBivariateSpline(y, x, Z, kx=1, ky=1)
        return interp_spline


    interpolador= bivariate_interpolation(my_df)

    def chi2_func(y,x):
        return interpolador(y, x, grid=False)

    FunVec=np.vectorize(chi2_func)

    y_value=np.unique(pty)
    x_value=np.unique(ptx)
    X_grid, Y_grid = np.meshgrid(x_value, y_value)
    Z_grid= FunVec(Y_grid,X_grid)
    
    return X_grid, Y_grid, Z_grid