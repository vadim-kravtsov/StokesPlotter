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

rad = pi/180
i = 0
e = 0
omega = 0
lp = 0
tau = 0

filt = 'b'
fiQ, q_obs, q_std = np.loadtxt(f'q_{filt}.txt', unpack=True)
fiU, u_obs, u_std = np.loadtxt(f'u_{filt}.txt', unpack=True)
q0 = mean(q_obs)
u0 = mean(u_obs)


def Q(l, tau, i, e, lp):
    S1 = e**2/4*cos(2*lp)
    S2 = e*cos(l+lp)
    S3 = (1+e**2/2)*cos(2*l)
    S4 = 3*e*cos(3*l-lp)
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


def calc_q(l, omega, i, e, lp, tau):
    l, omega, i, e, lp, tau = l, omega, i, e, lp, tau
    return q0+Q(l, tau, i, e, lp)*cos(omega)-U(l, tau, i, e, lp)*sin(omega)


def calc_u(l, omega, i, e, lp, tau):
    l, omega, i, e, lp, tau = l, omega, i, e, lp, tau
    return u0+Q(l, tau, i, e, lp)*sin(omega)+U(l, tau, i, e, lp)*cos(omega)


# Set up data
N = 200
l = np.linspace(0, 4*pi, N)
x = l
q = calc_q(l, e=e, i=i, omega=omega, lp=lp, tau=tau)
u = calc_q(l, e=e, i=i, omega=omega, lp=lp, tau=tau)
source = ColumnDataSource(data=dict(x=x, q=q, u=u))


# Set up plot
plot = figure(plot_height=400, plot_width=400, title="my sine wave",
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[0, 4*pi], y_range=[min(q_obs)-3*mean(q_std), min(q_obs)+3*mean(q_std)])


plot.line('x', 'q', source=source, line_width=3, line_alpha=0.6)
plot.circle(fiQ*2*pi, q_obs, size=3, color='#FF3333')
plot.circle(fiQ*2*pi+2*pi, q_obs, size=3, color='#FF3333')

plot2 = figure(plot_height=400, plot_width=400, title="my sine wave",
               tools="crosshair,pan,reset,save,wheel_zoom",
               x_range=[0, 4*pi], y_range=[min(u_obs)-3*mean(u_std), min(u_obs)+3*mean(u_std)])

plot2.line('x', 'u', source=source, line_width=3, line_alpha=0.6)
plot2.circle(fiU*2*pi, u_obs, size=3, color='#FF3333')
plot2.circle(fiU*2*pi+2*pi, u_obs, size=3, color='#FF3333')


# Set up widgets
e_slider = Slider(title="e", value=0.0, start=0, end=1.0, step=0.01)
i_slider = Slider(title="i", value=0.0, start=0, end=360, step=1)
omega_slider = Slider(title="omega", value=0.0, start=0, end=360, step=1)
lp_slider = Slider(title="lambda_p", value=0.0, start=0, end=360, step=1)
tau_slider = Slider(title="tau", value=0.0, start=0, end=0.1, step=0.001)


# Set up callbacks
plot.title.text = 'Stokes Q'
plot2.title.text = 'Stokes U'


def update_data(attrname, old, new):

    # Get the current slider values
    e = e_slider.value
    i = i_slider.value*rad
    omega = omega_slider.value*rad
    lp = lp_slider.value*rad
    tau = tau_slider.value
    
    # Generate the new curve
    q = calc_q(l, e=e, i=i, omega=omega, lp=lp, tau=tau)
    u = calc_u(l, e=e, i=i, omega=omega, lp=lp, tau=tau)

    source.data = dict(x=l, q=q, u=u)


for w in [e_slider, i_slider, omega_slider, lp_slider, tau_slider]:
    w.on_change('value', update_data)


# Set up layouts and add to document
inputs = column(e_slider, i_slider, omega_slider, lp_slider, tau_slider)

curdoc().add_root(row(plot, plot2, inputs,  width=1200))
curdoc().title = "Stokes Parameters"
