''' Stokes parameters in binay sisyems -- BME + disk. (Kravtsov et al. 2020 + disk)
    Use the ``bokeh serve`` command to run the example by executing:
    bokeh serve precessing_disk.py at your command prompt. Then navigate to the URL
    http://localhost:5006/sliders in your browser.
'''

import numpy as np
from numpy import mean, pi, cos, sin, array, sqrt, arctan2, tan, abs
from bokeh.io import curdoc
from bokeh.layouts import row, column
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, TextInput
from bokeh.plotting import figure
from bokeh.models import LinearColorMapper
import pickle
import pandas as pd

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
f0_disk     = 1     # typical scattering fraction of the disk
f0_cloud    = 0     # typical scattering fraction of the cloud

P_orb       = 5.59983     # orbital period
P_sup       = 13*P_orb    # superorbital period


def calc_orbital_phase(jd, jd0=2441874.707, period=P_orb):
    jd = np.array(jd)
    dt = jd - jd0
    phi = dt/period - np.array(dt/period).astype(int)
    return phi

def calc_superorbital_phase(jd, jd0=2441874.707, period=P_sup):
    jd = np.array(jd)
    dt = jd - jd0
    phi = dt/period - np.array(dt/period).astype(int)
    return phi
      

def bin_not_data(filt='B', bin_size_minutes=30, cpa=np.pi/180.0*np.array([32.49, 39., 43.7]), callibrated=False):
    # db_not = pd.read_csv(f'polar_cyg/final/not.txt', delimiter='\s+', header=0)
    if callibrated:
        db_not_individual = pickle.load(open(f'polar_cyg/final/not_callibrated_{filt}.txt', 'rb'))
        cpa = np.pi/180.0*np.array([34.8, 32.39, 31.4])
    
        cpa = {'B': cpa[0], 'V': cpa[1], 'R': cpa[2]}
        q_inst = {'B' :-0.00171, 'V': -0.00197, 'R': -0.00446}
        u_inst = {'B' :-0.00287, 'V': -0.00028, 'R': 0.00102}
  
    else:
        db_not_individual = pd.read_csv(f'polar_cyg/final/not_{filt}_individual.txt', delimiter='\s+', header=0)
        cpa = {'B': cpa[0], 'V': cpa[1], 'R': cpa[2]}
        q_inst = {'B' :-0.0, 'V': .0, 'R': -0.}
        u_inst = {'B' :-0.0, 'V': .0, 'R': 0.}
    
    db_not_individual = db_not_individual[db_not_individual.qe < 1]
    db_not_individual = db_not_individual[db_not_individual.ue < 1]
    
    db_not_individual.q -= q_inst[filt]
    db_not_individual.u -= u_inst[filt]
    
    p = np.sqrt(db_not_individual.q**2 + db_not_individual.u**2)
    pa = 0.5*np.arctan2(db_not_individual.u, db_not_individual.q) + cpa[filt]
    
    q = p*np.cos(2*pa)
    u = p*np.sin(2*pa)
    
    db_not_individual.q = q
    db_not_individual.u = u
    
    night_numbers = np.unique(db_not_individual.NO)
 
    binned_dataset_not = pd.DataFrame({'jd': [], 
                                   'q': [],
                                   'u': [],
                                   'qe': [],
                                   'ue': []})
    
    
    bin_size_minutes = bin_size_minutes
    bin_size_days = bin_size_minutes/(60*24)
       
    for night_number in night_numbers:
        one_night_data = db_not_individual[db_not_individual.NO == night_number] 
     
        start, finish = min(one_night_data.jd), max(one_night_data.jd)
        duration_of_observtations = finish - start
        
        number_of_bins = int((duration_of_observtations/bin_size_days))
        
        for bin_number in range(number_of_bins):
            bin_start = start + bin_number*bin_size_days
            bin_finish = bin_start + bin_size_days
            
            bin_size = bin_finish - bin_start
            bin_center = bin_start + (bin_size)/2.0
         
            one_bin_data = one_night_data[one_night_data.jd >= bin_start]
            one_bin_data = one_bin_data[one_bin_data.jd < bin_finish]
            
            if len(one_bin_data > 0):
                q = np.average(one_bin_data.q, weights=1/one_bin_data.qe)
                u = np.average(one_bin_data.u, weights=1/one_bin_data.ue)
                
                q_err = np.std(one_bin_data.q)/np.sqrt(len(one_bin_data.q))
                u_err = np.std(one_bin_data.u)/np.sqrt(len(one_bin_data.u))
        
                phase = calc_orbital_phase(bin_center + 2400000)
                
                if q_err == 0: 
                    continue
                   
                if u_err == 0:
                    continue 
                    
                bin_df = pd.DataFrame({'jd': [bin_center + 2400000], 
                                   'q': [q],
                                   'u': [u],
                                   'qe': [q_err],
                                   'ue': [u_err]})
                
                binned_dataset_not = binned_dataset_not.append(bin_df, ignore_index=True)            
                
            else:
                continue
            
    return binned_dataset_not
    
   
def bin_t60_data(filt='B', bin_size_minutes=30):
    db_mean = pd.read_csv(f'polar_cyg/final/{filt}.txt', delimiter='\s+', header=0)
    
    db_koko = pd.read_csv(f'polar_cyg/final/{filt}_individual.txt', delimiter='\s+', header=0)
    db_koko = db_koko[db_koko.qe < 1]
    db_koko = db_koko[db_koko.ue < 1]
    
    cpa = np.pi/180.0*np.array([34.8, 32.39, 31.4])
    
    cpa = {'B': cpa[0], 'V': cpa[1], 'R': cpa[2]}
    q_inst = {'B' :-0.00171, 'V': -0.00197, 'R': -0.00446}
    u_inst = {'B' :-0.00287, 'V': -0.00028, 'R': 0.00102}
    
    db_koko.q -= q_inst[filt]
    db_koko.u -= u_inst[filt]
    
    p = np.sqrt(db_koko.q**2 + db_koko.u**2)
    pa = 0.5*np.arctan2(db_koko.u, db_koko.q) + cpa[filt]
    
    q = p*np.cos(2*pa)
    u = p*np.sin(2*pa)
    
    db_koko.q = q
    db_koko.u = u
    
    night_numbers = np.unique(db_koko.NO)
    
    binned_dataset = pd.DataFrame({'jd': [], 
                                   'q': [],
                                   'u': [],
                                   'qe': [],
                                   'ue': []})
    
    
    bin_size_minutes = bin_size_minutes
    bin_size_days = bin_size_minutes/(60*24)
       
    for night_number in night_numbers:
        one_night_data = db_koko[db_koko.NO == night_number] 
     
        start, finish = min(one_night_data.jd), max(one_night_data.jd)
        duration_of_observtations = finish - start
        
        number_of_bins = int((duration_of_observtations/bin_size_days))
        
        for bin_number in range(number_of_bins):
            bin_start = start + bin_number*bin_size_days
            bin_finish = bin_start + bin_size_days
            
            bin_size = bin_finish - bin_start
            bin_center = bin_start + (bin_size)/2.0
         
            one_bin_data = one_night_data[one_night_data.jd >= bin_start]
            one_bin_data = one_bin_data[one_bin_data.jd < bin_finish]
            
            if len(one_bin_data > 0):
                q = np.average(one_bin_data.q, weights=1/one_bin_data.qe)
                u = np.average(one_bin_data.u, weights=1/one_bin_data.ue)
                
                q_err = np.std(one_bin_data.q)/np.sqrt(len(one_bin_data.q))
                u_err = np.std(one_bin_data.u)/np.sqrt(len(one_bin_data.u))
        
                phase = calc_orbital_phase(bin_center + 2400000)
                
                if q_err == 0: 
                    continue
                    q_err = 0.03
                
                if u_err == 0:
                    continue 
                    u_err = 0.03
                    
                bin_df = pd.DataFrame({'jd': [bin_center + 2400000], 
                                   'q': [q],
                                   'u': [u],
                                   'qe': [q_err],
                                   'ue': [u_err]})
                
                binned_dataset = binned_dataset.append(bin_df, ignore_index=True)                
            else:
                continue
    return binned_dataset


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


# Set up data
N = 2000
jd = np.linspace(2441874.707, 2441874.707 + 2*P_orb, N)
phi = calc_orbital_phase(jd)
phi_sup = calc_superorbital_phase(jd)
rad = pi/180

'''
filt = 'B'

not_binned = bin_not_data(filt=filt, bin_size_minutes=30, callibrated=True)
t60_binned = bin_t60_data(filt=filt, bin_size_minutes=30)
  
    
set = pd.concat([not_binned, t60_binned])    
set = set.sort_values(by='jd')
        
jd = set.jd
phi = calc_orbital_phase(jd)
phi_sup = calc_superorbital_phase(jd)
'''
color_mapper = LinearColorMapper(palette='Turbo256', low=min(jd), high=max(jd))

q, u = thomson_model(phi, phi_sup, e, i*rad, Omega*rad, beta*rad, l_p*rad, phi_0, gamma_0, q0, u0, f0_disk, f0_cloud)
# q += 0.4
# u += -4.5
p = np.sqrt(q**2 + u**2)
pa = 0.5*np.arctan2(u, q)*180/np.pi + 180
   
source = ColumnDataSource(data=dict(x=phi, q=q, u=u, jd=jd, p=p, pa=pa))


# Set up plot
plotQ = figure(plot_height=500, plot_width=500, title="my sine wave",
              tools="crosshair,pan,reset,save,wheel_zoom, box_zoom",
              x_range=[0, 1], y_range=[-0.5, 0.5])


plotQ.scatter('x', 'q', source=source, line_width=1, line_alpha=1, color={'field': 'jd', 'transform': color_mapper})
plotQ.ray(x=[0.25], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
plotQ.ray(x=[0.50], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
plotQ.ray(x=[0.75], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)

plotU = figure(plot_height=500, plot_width=500, title="my sine wave",
               tools="crosshair,pan,reset,save,wheel_zoom, box_zoom",
               x_range=[0, 1], y_range=[-0.5, 0.5])

plotU.ray(x=[0.25], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
plotU.ray(x=[0.50], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
plotU.ray(x=[0.75], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)


plotQU = figure(plot_height=500, plot_width=500, title="QU",
              tools="crosshair,pan,reset,save,wheel_zoom, box_zoom",
              x_range=[-1, 1], y_range=[-1, 1])

plotQU.scatter('q', 'u', source=source, line_width=1, line_alpha=1,  color={'field': 'jd', 'transform': color_mapper})


'''
not_binned = bin_not_data(filt='R', bin_size_minutes=60, callibrated=True)
t60_binned = bin_t60_data(filt='R', bin_size_minutes=60)
  

full_set = pd.concat([not_binned, t60_binned])
full_set = full_set.sort_values(by='jd')
 
qu_ism = {'B':{'q':-0.09, 'u':-4.31, 'err':0.17},
          'V':{'q':-0.12, 'u':-3.91, 'err':0.14},
          'R':{'q':-0.19, 'u':-3.82, 'err':0.15}}

    
  
# full_set.q -= qu_ism['B']['q']
# full_set.u -= qu_ism['B']['u']
    
set1 = full_set[(full_set.jd > 2459676) & (full_set.jd < 2459694)]
set2 = full_set[(full_set.jd > 2459715) & (full_set.jd < 2459723)]
set3 = full_set[(full_set.jd > 2459754) & (full_set.jd < 2459794)]


data_set2 = ColumnDataSource(data=dict(x=calc_orbital_phase(set2.jd), q=set2.q, u=set2.u, jd=set2.jd))

plotQU.scatter('q', 'u', source=data_set2, line_alpha=1, line_width=1)
plotQ.scatter('x', 'q', source=data_set2, line_alpha=1, line_width=1)
plotU.scatter('x', 'u', source=data_set2, line_alpha=1, line_width=1)
'''

# plotP = figure(plot_height=500, plot_width=500, title="my sine wave",
#               tools="crosshair,pan,reset,save,wheel_zoom, box_zoom",
#               x_range=[0, 1], y_range=[0, 2.5])


# plotP.scatter('x', 'p', source=source, line_width=1, line_alpha=1, color={'field': 'jd', 'transform': color_mapper})
# plotP.ray(x=[0.25], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
# plotP.ray(x=[0.50], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
# plotP.ray(x=[0.75], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)


# plotPA = figure(plot_height=500, plot_width=500, title="my sine wave",
#               tools="crosshair,pan,reset,save,wheel_zoom, box_zoom",
#               x_range=[0, 1], y_range = [120, 180])


# plotPA.scatter('x', 'pa', source=source, line_width=1, line_alpha=1, color={'field': 'jd', 'transform': color_mapper})
# plotPA.ray(x=[0.25], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
# plotPA.ray(x=[0.50], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)
# plotPA.ray(x=[0.75], y=-5, length=10, angle=pi/2, line_width=1, alpha=0.3)


plotU.scatter('x', 'u', source=source, line_width=1, line_alpha=1,  color={'field': 'jd', 'transform': color_mapper})


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
f0_disk_slider    = Slider(title="f0_disk", value=f0_disk, start=0, end=1.0, step=0.01)
f0_cloud_slider   = Slider(title="f0_cloud", value=f0_cloud, start=0, end=1.0, step=0.01)
p_sup_frac_slider = Slider(title="p_sup_in_orb", value=P_sup/P_orb, start=1, end=100, step=1)
multiplier_slider = Slider(title="p_sup_mult", value=1, start=1, end=101, step=100)

# Set up callbacks
plotQ.title.text = 'Stokes Q'
plotU.title.text = 'Stokes U'

rad = np.pi/180


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
    p_sup_frac = p_sup_frac_slider.value
    multiplier = multiplier_slider.value
    phi_sup = calc_superorbital_phase(jd, period=p_sup_frac*P_orb*multiplier)
    
    q, u = thomson_model(phi + phi_0, phi_sup, e, i, Omega, beta, l_p, phi_0, gamma_0, q0, u0, f0_disk, f0_cloud)
    # q += 0.4
    # u += -4.56
    # q += np.random.normal(loc=0, scale=0.01, size=len(q))
    # u += np.random.normal(loc=0, scale=0.01, size=len(u))
    p = np.sqrt(q**2 + u**2)
    pa = 0.5*np.arctan2(u, q)*180/np.pi + 180
    source.data = dict(x=phi + phi_0, q=q, u=u, jd=jd, p=p, pa=pa)


for w in [e_slider, i_slider ,omega_slider, lp_slider, q0_slider, u0_slider,
          beta_slider, phi_0_slider, gamma_0_slider, f0_disk_slider, f0_cloud_slider, p_sup_frac_slider, multiplier_slider]:
    w.on_change('value', update_data)


# Set up layouts and add to document
inputs = column(e_slider, i_slider ,omega_slider, lp_slider, q0_slider, u0_slider,
                beta_slider, phi_0_slider, gamma_0_slider, f0_disk_slider, f0_cloud_slider, p_sup_frac_slider, multiplier_slider)

# curdoc().add_root(column(row(plotQ, plotU, inputs,  width=1000), row(plotP, plotPA, width=1000)))
curdoc().add_root(column(row(plotQ, plotU, inputs,  width=1000), row(plotQU, width=500)))

#curdoc().add_root(row(plotP, plotAng,  width=1300))

curdoc().title = "Stokes Parameters"
