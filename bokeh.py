''' Stokes parameters in binay sisyems. (Simmons, Boyle, 1984)
    Use the ``bokeh serve`` command to run the example by executing:
    bokeh serve sliders.py at your command prompt. Then navigate to the URL
    http://localhost:5006/sliders in your browser.
'''

import numpy as np
from numpy import mean, pi, cos, sin, array
from bokeh.io import curdoc
from bokeh.layouts import row, column
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, TextInput
from bokeh.plotting import figure


e0, i0, omega0, lp0, tau0, q0_0, u0_0 = 0.5, 40, 0, 40, 0.02, 0.52, -1.21

rad = pi/180
i = 0
e = 0.5
omega = 0
lp = 40
tau = 0.01

filt = 'v'
# fiQ, q_obs, q_std = np.loadtxt(f'q_{filt}.txt', unpack=True)
# fiU, u_obs, u_std = np.loadtxt(f'u_{filt}.txt', unpack=True)

fi, q_obs, q_std, u_obs, u_std = np.loadtxt(f'{filt}_all.csv', unpack=True, skiprows=1, delimiter=',')
fiQ, fiU = fi, fi

q0 = mean(q_obs)
u0 = mean(u_obs)
#q0 = 0
#u0 = 0


def Q_berdyugin(l, tau, i, e, lp):
    return -tau/(1-e**2)**2*(1+e*cos(l-lp))**2*((1+cos(i)**2)*cos(2*l)-sin(i)**2)


def U_berdyugin(l, tau, i, e, lp):
    return -2*tau/(1-e**2)**2*(1+e*cos(l-lp))**2*cos(i)*(sin(2*l))


def calc_q_berdyugin(l, q0, omega, i, e, lp, tau):
    l, omega, i, e, lp, tau = l, omega, i, e, lp, tau
    return Q_berdyugin(l, tau, i, e, lp)*cos(omega)-U_berdyugin(l, tau, i, e, lp)*sin(omega)+q0


def calc_u_berdyugin(l, u0, omega, i, e, lp, tau):
    l, omega, i, e, lp, tau = l, omega, i, e, lp, tau
    return U_berdyugin(l, tau, i, e, lp)*sin(omega)+Q_berdyugin(l, tau, i, e, lp)*cos(omega)+u0


def Q(l, tau, i, e, lp):
    S1 = e**2/4*cos(2*lp)
    S2 = e*cos(l+lp)
    S3 = (1+e**2/2)*cos(2*l)
    S4 = e*cos(3*l-lp)
    S5 = e**2/4*cos(2*(2*l-lp))
    C1 = (1+e**2/2)
    C2 = 2*e*cos(l-lp)*e**2/2*cos(2*(l-lp))
    return -tau/(1-e**2)**2*((1+cos(i)**2)*(S1+S2+S3+S4+S5)) - (sin(i)**2*(C1+C2))  # Simmons, Boyle 1984


def U(l, tau, i, e, lp):
    S1 = e**2/4*sin(2*lp)
    S2 = e*sin(l+lp)
    S3 = (1+e**2/2)*sin(2*l)
    S4 = e*sin(3*l-lp)
    S5 = e**2/4*sin(2*(2*l - lp))
    return -tau/(1-e**2/2)**2*2*cos(i)*(S1+S2+S3+S4+S5)


def calc_q(l, q0, omega, i, e, lp, tau):
    l, omega, i, e, lp, tau = l, omega, i, e, lp, tau
    return Q(l, tau, i, e, lp)*cos(omega)-U(l, tau, i, e, lp)*sin(omega)+q0


def calc_u(l, u0, omega, i, e, lp, tau):
    l, omega, i, e, lp, tau = l, omega, i, e, lp, tau
    return U(l, tau, i, e, lp)*sin(omega)+Q(l, tau, i, e, lp)*cos(omega)+u0


# Set up data
N = 200
l = np.linspace(0, 4*pi, N)
x = l
phase_reduction = 0.275*2*pi
#q_obs = q_obs - mean(q0)
#u_obs = u_obs - mean(u0)

q = calc_q_berdyugin(l, q0, e=e, i=i, omega=omega, lp=lp, tau=tau)
u = calc_u_berdyugin(l, u0, e=e, i=i, omega=omega, lp=lp, tau=tau)
p = np.sqrt(q**2+u**2)
ang = 1/2*np.arctan2(u, q)
source = ColumnDataSource(data=dict(x=x, q=q, u=u, p=p, ang=ang))


# Set up plot
plot = figure(plot_height=500, plot_width=500, title="my sine wave",
              tools="crosshair,pan,reset,save,wheel_zoom, box_zoom",
              x_range=[0, 4*pi], y_range=[-0.07, 0.17])
              #x_range=[0, 4*pi], y_range=[-1.5, 1.5])


plot.line('x', 'q', source=source, line_width=3, line_alpha=0.6)
plot.circle(fiQ*2*pi - phase_reduction, q_obs, size=7, color='#FF3333')
plot.circle(fiQ*2*pi+2*pi - phase_reduction, q_obs, size=7, color='#FF3333')

plot2 = figure(plot_height=500, plot_width=500, title="my sine wave",
               tools="crosshair,pan,reset,save,wheel_zoom, box_zoom",
               #x_range=[0, 4*pi], y_range=[mean(u_obs)-5*mean(u_std), mean(u_obs)+5*mean(u_std)])
               x_range=[0, 4*pi], y_range=[-1.35, -1.05])

plot2.line('x', 'u', source=source, line_width=3, line_alpha=0.6)
plot2.circle(fiU*2*pi - phase_reduction, u_obs, size=7, color='#FF3333')
plot2.circle(fiU*2*pi+2*pi - phase_reduction, u_obs, size=7, color='#FF3333')

plotP = figure(plot_height=500, plot_width=500, title="my sine wave",
               tools="crosshair,pan,reset,save,wheel_zoom, box_zoom",
               x_range=[0, 4*pi], y_range=[mean(p)-5*mean(p), mean(p)+5*mean(p)])
               #x_range=[0, 4*pi], y_range=[-1.5, 1.5])

plotP.line('x', 'p', source=source, line_width=3, line_alpha=0.6)
plotP.circle(fiU*2*pi - phase_reduction, np.sqrt(q_obs**2+u_obs**2), size=7)
plotP.circle(fiU*2*pi+2*pi - phase_reduction, np.sqrt(q_obs**2+u_obs**2), size=7)


plotAng = figure(plot_height=500, plot_width=500, title="my sine wave",
               tools="crosshair,pan,reset,save,wheel_zoom, box_zoom",
               x_range=[0, 4*pi], y_range=[-pi, pi])
               #x_range=[0, 4*pi], y_range=[-1.5, 1.5])

plotAng.line('x', 'ang', source=source, line_width=3, line_alpha=0.6)
plotAng.circle(fiU*2*pi - phase_reduction, (np.arctan2(u_obs, q_obs)/pi)%1*pi-pi/2, size=7)
plotAng.circle(fiU*2*pi+2*pi - phase_reduction, (np.arctan2(u_obs, q_obs)/pi)%1*pi-pi/2, size=7)


# Set up widgets
e_slider = Slider(title="e", value=e0, start=0, end=1.0, step=0.01)
i_slider = Slider(title="i", value=i0, start=-90, end=90, step=0.1)
omega_slider = Slider(title="omega", value=omega0, start=-180, end=180, step=0.1)
lp_slider = Slider(title="lambda_p", value=lp0, start=0, end=360, step=0.1)
tau_slider = Slider(title="tau", value=tau0, start=0, end=0.05, step=0.0001)
q0_slider = Slider(title="q0", value=q0_0, start=-2, end=2, step=0.0005)
u0_slider = Slider(title="u0", value=u0_0, start=-2, end=2, step=0.0005)


# Set up callbacks
plot.title.text = 'Stokes Q'
plot2.title.text = 'Stokes U'
plotP.title.text = 'Polarization'
plotAng.title.text = 'Angle'


def update_data(attrname, old, new):

    # Get the current slider values
    e = e_slider.value
    i = i_slider.value*rad
    omega = omega_slider.value*rad
    lp = lp_slider.value*rad
    tau = tau_slider.value
    q0 = q0_slider.value
    u0 = u0_slider.value

    # Generate the new curve
    #q = calc_q_berdyugin(l, q0, e=e, i=i, omega=omega, lp=lp, tau=tau)
    #u = calc_u_berdyugin(l, u0, e=e, i=i, omega=omega, lp=lp, tau=tau)
    q = calc_q(l, q0, e=e, i=i, omega=omega, lp=lp, tau=tau)
    u = calc_u(l, u0, e=e, i=i, omega=omega, lp=lp, tau=tau)
    p = np.sqrt(q**2+u**2)
    ang = 1/2*np.arctan2(u, q)
    source.data = dict(x=l, q=q, u=u, p=p, ang=ang)


for w in [e_slider, i_slider, omega_slider, lp_slider, tau_slider, q0_slider, u0_slider]:
    w.on_change('value', update_data)


# Set up layouts and add to document
inputs = column(e_slider, i_slider, omega_slider, lp_slider, tau_slider, q0_slider, u0_slider)

curdoc().add_root(row(plot, plot2, inputs,  width=1300))
#curdoc().add_root(row(plotP, plotAng,  width=1300))

curdoc().title = "Stokes Parameters"
