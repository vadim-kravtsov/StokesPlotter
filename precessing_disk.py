''' Stokes parameters in binay sisyems -- BME + disk. (Kravtsov et al. 2020 + disk)
    Use the ``bokeh serve`` command to run the example by executing:
    bokeh serve precessing_disk.py at your command prompt. Then navigate to the URL
    http://localhost:5006/sliders in your browser.
'''

import numpy as np
from numpy import mean, pi, cos, sin, array, sqrt, arctan2, tan, abs
from bokeh.io import curdoc
from bokeh.layouts import row, column
from bokeh.plotting import figure
from bokeh.models import LinearColorMapper
import pickle
import pandas as pd

from bokeh.models import ColumnDataSource, BoxZoomTool, CrosshairTool, ResetTool
from bokeh.models import SaveTool, PanTool, FileInput, Div, DataTable, TableColumn
from bokeh.models.widgets import Slider, TextInput, Button
from sqlalchemy import null

# from plotter import bin_not_data, bin_t60_data

# Parameters of the model (angles in degrees)

e           = 0.0   # eccentricity 
i           = 0     # inclination of the orbit 
Omega       = 0     # position angle of the orbit on the sky
beta        = 0     # precession angle
l_p         = 0     # longitude of periastron
phi_0       = 0     # orbital phase of periastron
gamma_0     = 0     # superorbital phase of some special orientation of the disk axis
q0          = 0     # constant level of Q
u0          = 0     # constant level of U
f0_disk     = 0     # typical scattering fraction of the disk
f0_cloud    = 0.5     # typical scattering fraction of the cloud

# P_orb       = 5.59983     # orbital period
# P_sup       = 13*P_orb    # superorbital period

rad = np.pi/180.0

def read_data(filename):
    try:
        data = pd.read_csv(filename,
                           delimiter='\s+',
                           names=['phi', 'q', 'u', 'error'], skiprows=1)
        exitcode = 0
        return data, exitcode 
    except FileNotFoundError:
        exitcode = 1
        return None, exitcode

def errorbar(fig, x, y, xerr=None, yerr=None, color='red', 
             point_kwargs={}, error_kwargs={}):
             
    fig.circle(x, y, color=color, **point_kwargs)

#   x = np.array(x)
#   y = np.array(y)
#   xerr = np.array()
    y_err_x = []
    y_err_y = []
    for px, py, err in zip(x, y, yerr):
        y_err_x.append((px, px))
        y_err_y.append((py - err, py + err))
    fig.multi_line(y_err_x, y_err_y, color=color, **error_kwargs)


def calc_lambda_from_phi(phi, e, l_p, phi_0):
    M = 2*pi*(phi - phi_0)
    E = M
    for _ in range(100):
        E = M + e*sin(E)
    return 2. * arctan2(sqrt(1. + e) * tan(E/2.), sqrt(1. - e)) + l_p

    
def thomson_model(phi, phi_sup, e, i, Omega, beta, l_p, phi_0, gamma_0, q0, u0, f0_disk, f0_cloud):
    l = calc_lambda_from_phi(phi, e, l_p, phi_0)
    
    gamma = 2*pi*phi_sup + gamma_0
    cos_Psi = sin(l)*sin(beta)*sin(gamma) + cos(l)*sin(beta)*cos(gamma)
    cos_Sigma = cos(i)*cos(beta) + sin(i)*sin(beta)*cos(gamma)
    
    f0_disk_array = []

    for indx in range(len(l)):
        if beta <= pi/2.:
            if cos_Sigma[indx] > 0:
                # See the top of the disk
                if cos_Psi[indx] <= 0:
                    # The top of the disk is bright
                    f0_current = f0_disk
                else:
                    # The bottom of the disk is bright
                    f0_current = 0
            elif cos_Sigma[indx] < 0:
                # See the bottom of the disk
                if cos_Psi[indx] <= 0:
                    # The top of the disk is bright
                    f0_current = 0
                else:
                    # The bottom of the disk is bright
                    f0_current = f0_disk
            else:
                # See the disk from the edge
                f0_current = 0
        else:
            if cos_Sigma[indx] < 0:
                # See the top of the disk
                if cos_Psi[indx] >= 0:
                    # The top of the disk is bright
                    f0_current = f0_disk
                else:
                    # The bottom of the disk is bright
                    f0_current = 0
            elif cos_Sigma[indx] > 0:
                # See the bottom of the disk
                if cos_Psi[indx] >= 0:
                    # The top of the disk is bright
                    f0_current = 0
                else:
                    # The bottom of the disk is bright
                    f0_current = f0_disk
            else:
                # See the disk from the edge
                f0_current = 0
        f0_disk_array.append(f0_current)
 
    q_part = disk_and_cloud_q(l, l_p, i, e, f0_disk_array, f0_cloud, cos_Psi, cos_Sigma)
    u_part = disk_and_cloud_u(l, l_p, i, e, f0_disk_array, f0_cloud, cos_Psi, cos_Sigma)
        
    q_res = q_part*cos(2.*Omega) - u_part*sin(2.*Omega)
    u_res = q_part*sin(2.*Omega) + u_part*cos(2.*Omega)
        
    q_res += q0
    u_res += u0
        
    return (q_res, u_res)

    
def disk_and_cloud_q(l, l_p, i, e, f0_disk_array, f0_cloud, cos_Psi, cos_Sigma):
    f_sc_disk = f0_disk_array * (1 + e * cos(l - l_p))**2 * abs(cos_Psi)
    f_sc_cloud = f0_cloud * (1 + e * cos(l - l_p))**2
    return 3./16. * (sin(i)**2 - (1. + cos(i)**2) * cos(2.*l))*(f_sc_cloud + (1. - f_sc_cloud)*f_sc_disk)
    
    
def disk_and_cloud_u(l, l_p, i, e, f0_disk_array, f0_cloud, cos_Psi, cos_Sigma):
    f_sc_disk = f0_disk_array * (1 + e * cos(l - l_p))**2 * abs(cos_Psi)
    f_sc_cloud = f0_cloud * (1 + e * cos(l - l_p))**2
    return - 3./8. * cos(i)*sin(2.*l) * (f_sc_cloud + (1. - f_sc_cloud)*f_sc_disk)

            
def main(): 
# phi = original_dataframe.data['phi']
# phi_sup = np.zeros(len(phi))

# print(phi)
    npoints = 100
    rad = np.pi/180
    phi = np.linspace(0, 1, npoints)
    phi_sup = np.zeros(npoints)
    q, u = thomson_model(phi, phi_sup, e, i*rad, Omega*rad, beta*rad, l_p*rad, phi_0, gamma_0, q0, u0, f0_disk, f0_cloud)
# q += 0.4
# u += -4.5
    p = np.sqrt(q**2 + u**2)
    pa = 0.5*np.arctan2(u, q)*180/np.pi
    
    source = ColumnDataSource(data=dict(x=phi, q=q, u=u, p=p, pa=pa))


# Set up plot
    plotQ = figure(plot_height=300, plot_width=500, title="my sine wave",
              tools="crosshair,pan,reset,save,wheel_zoom, box_zoom",
              x_range=[0, 1], y_range=[-0.5, 0.5])


# plotQ.scatter('x', 'q', source=source, line_width=1, line_alpha=1)

    plotQ.ray(x=[0.25], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
    plotQ.ray(x=[0.50], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
    plotQ.ray(x=[0.75], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)

    plotU = figure(plot_height=300, plot_width=500, title="my sine wave",
               tools="crosshair,pan,reset,save,wheel_zoom, box_zoom",
               x_range=[0, 1], y_range=[-0.5, 0.5])

    plotU.ray(x=[0.25], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
    plotU.ray(x=[0.50], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
    plotU.ray(x=[0.75], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)

    
    plotP = figure(plot_height=300, plot_width=500, title="PD",
               tools="crosshair,pan,reset,save,wheel_zoom, box_zoom",
               x_range=[0, 1], y_range=[0, 10])

    plotP.ray(x=[0.25], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
    plotP.ray(x=[0.50], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
    plotP.ray(x=[0.75], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)


    plotPA = figure(plot_height=300, plot_width=500, title="PD",
               tools="crosshair,pan,reset,save,wheel_zoom, box_zoom",
               x_range=[0, 1], y_range=[-180, 180])

    plotPA.ray(x=[0.25], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
    plotPA.ray(x=[0.50], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
    plotPA.ray(x=[0.75], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
    
    
    plotQU = figure(plot_height=300, plot_width=500, title="QU",
              tools="crosshair,pan,reset,save,wheel_zoom, box_zoom",
              x_range=[-1, 1], y_range=[-1, 1])

# plotQU.scatter('q', 'u', source=source, line_width=1, line_alpha=1)
    plotQ.line('x', 'q', source=source, line_width=5, line_alpha=1)
    plotU.line('x', 'u', source=source, line_width=5, line_alpha=1)
    plotP.line('x', 'p', source=source, line_width=5, line_alpha=1)
    plotPA.line('x', 'pa', source=source, line_width=5, line_alpha=1)

    global original_dataframe
    original_dataframe = ColumnDataSource({'phi': [], 'q' : [], 'u' : [], 'error' : []})

    errors_div = Div(text="")    
    def file_input_handler(attr, old, new):
        filename = new
        data, exitcode = read_data(filename)
        if exitcode == 1: 
            errors_div.text = 'Error: data file should be in the same folder with the python script!'
            errors_div.style = {'color': 'red'}
        else:
            errors_div.text = ''
            errors_div.style = {'color': 'black'}
            original_dataframe.data = original_dataframe.from_df(data)
                
    
    def plot_data(event):
        global original_dataframe
        try:
            print(original_dataframe.data)
        
            min_y, max_y = min(original_dataframe.data['q']), max(original_dataframe.data['q'])
            plotQ.y_range.start = min_y - 0.1*(max_y - min_y)
            plotQ.y_range.end = max_y + 0.1*(max_y - min_y)
            errorbar(plotQ, original_dataframe.data['phi'], original_dataframe.data['q'], yerr = np.array(original_dataframe.data['error']), color='red')
            
            min_y, max_y = min(original_dataframe.data['u']), max(original_dataframe.data['u'])
            plotU.y_range.start = min_y - 0.1*(max_y - min_y)
            plotU.y_range.end = max_y + 0.1*(max_y - min_y)
            errorbar(plotU, original_dataframe.data['phi'], original_dataframe.data['u'], yerr = np.array(original_dataframe.data['error']), color='red')
            
            p_obs = np.sqrt(original_dataframe.data['q']**2 + original_dataframe.data['u']**2)
            pa = 0.5*np.arctan2(original_dataframe.data['u'], original_dataframe.data['q'])*180/np.pi
            
            min_y, max_y = min(p_obs), max(p_obs)
            plotP.y_range.start = min_y - 0.1*(max_y - min_y)
            plotP.y_range.end = max_y + 0.1*(max_y - min_y)
            errorbar(plotP, original_dataframe.data['phi'], p_obs, yerr = np.array(original_dataframe.data['error']), color='red')
            
        except ValueError:
            errors_div.text = 'Error: open the data file first!'
            errors_div.style = {'color': 'red'}

    file_input = FileInput(accept=".csv,.txt", multiple=False)
    file_input.on_change('filename', file_input_handler)
   

# Set up widgets
    e_slider          = Slider(title="e", value=e, start=0, end=1.0, step=0.01)
    i_slider          = Slider(title="i", value=i, start=0, end=180, step=1)
    omega_slider      = Slider(title="omega", value=Omega, start=-180, end=180, step=1)
    lp_slider         = Slider(title="l_p", value=l_p, start=0, end=360, step=1)
    q0_slider         = Slider(title="q0", value=q0, start=-2, end=2, step=0.0005)
    u0_slider         = Slider(title="u0", value=u0, start=-2, end=2, step=0.0005)
    beta_slider       = Slider(title="beta", value=beta, start=0, end=180, step=1)
    phi_0_slider      = Slider(title="phi_0", value=phi_0, start=0, end=1.0, step=0.01)
    gamma_0_slider    = Slider(title="gamma_0", value=gamma_0, start=0, end=360, step=1)
    f0_disk_slider    = Slider(title="f0_disk", value=f0_disk, start=0, end=10.0, step=0.01)
    f0_cloud_slider   = Slider(title="f0_cloud", value=f0_cloud, start=0, end=10.0, step=0.01)
# p_sup_frac_slider = Slider(title="p_sup_in_orb", value=P_sup/P_orb, start=1, end=100, step=1)
# multiplier_slider = Slider(title="p_sup_mult", value=1, start=1, end=101, step=100)

    data_plot_button = Button(label='Plot data')
    data_plot_button.on_click(plot_data)
# Set up callbacks
    plotQ.title.text = 'Stokes Q'
    plotU.title.text = 'Stokes U'
    plotP.title.text = 'PD'
    plotPA.title.text = 'PA'

    def update_data(attrname, old, new):

        # Get the current slider values
        e = e_slider.value
        i = i_slider.value*rad
        Omega = omega_slider.value*rad
        beta = beta_slider.value*rad
        l_p = lp_slider.value*rad
        phi_0 = phi_0_slider.value
        gamma_0 = gamma_0_slider.value*rad
        q0 = q0_slider.value
        u0 = u0_slider.value
        f0_disk = f0_disk_slider.value
        f0_cloud = f0_cloud_slider.value
       
        q, u = thomson_model(phi + phi_0, phi_sup, e, i, Omega, beta, l_p, phi_0, gamma_0, q0, u0, f0_disk, f0_cloud)
       
        p = np.sqrt(q**2 + u**2)
        pa = 0.5*np.arctan2(u, q)*180/np.pi
        source.data = dict(x=phi + phi_0, q=q, u=u, p=p, pa=pa)
        

    for w in [e_slider, i_slider ,omega_slider, lp_slider, q0_slider, u0_slider,
          beta_slider, phi_0_slider, gamma_0_slider, f0_disk_slider, f0_cloud_slider]:
          w.on_change('value', update_data)


    # Set up layouts and add to document

    logo = Div(text="""<b>Stokes Plotter</b>""",
               style={'text-align': 'left', 'font-size': 'large'})
    
    instructions = Div(text="""Please, open ".txt" or ".csv" file with the data in the following format: <i>phi, q, u, err</i>.""")
    
    
    inputs = column(logo, instructions, file_input, data_plot_button, e_slider, i_slider ,omega_slider, lp_slider, q0_slider, u0_slider,
                beta_slider, phi_0_slider, gamma_0_slider, f0_disk_slider, f0_cloud_slider)

    curdoc().add_root(row(column(plotQ, plotU, width=500), 
                          column(plotP, plotPA, width=500),
                          column(inputs, errors_div, width=300)))
    
    curdoc().title = "Stokes Parameters"


main()