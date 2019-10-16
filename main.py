''' Present an interactive function explorer with slider widgets.
Scrub the sliders to change the properties of the ``sin`` curve, or
type into the title text box to update the title of the plot.
Use the ``bokeh serve`` command to run the example by executing:
    bokeh serve sliders.py
at your command prompt. Then navigate to the URL
    http://localhost:5006/sliders
in your browser.
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
q0 = 0
u0 = 0

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


def calc_q(l, omega=0, i=25, e=.573, lp=40.5, tau=0.1, q0=1):
    l, q0, omega, i, e, lp, tau = l, q0, omega, i, e, lp, tau
    return q0+Q(l, tau, i, e, lp)*cos(omega)-U(l, tau, i, e, lp)*sin(omega)


def calc_u(l, omega=0, i=25, e=.573, lp=40.5, tau=0.1, u0=1):
    l, u0, omega, i, e, lp, tau = l, u0, omega, i, e, lp, tau
    return u0+Q(l, tau, i, e, lp)*sin(omega)+U(l, tau, i, e, lp)*cos(omega)


# Set up data
N = 200
l = np.linspace(0, 4*pi, N)
x = l
q = calc_q(l, e=e, i=i, omega=omega, lp=lp)
u = calc_q(l, e=e, i=i, omega=omega, lp=lp)
source = ColumnDataSource(data=dict(x=x, q=q, u=u))


# Set up plot
plot = figure(plot_height=400, plot_width=400, title="my sine wave",
              tools="crosshair,pan,reset,save,wheel_zoom",
              x_range=[0, 4*pi], y_range=[-1, 3])

plot.line('x', 'q', source=source, line_width=3, line_alpha=0.6)


plot2 = figure(plot_height=400, plot_width=400, title="my sine wave",
               tools="crosshair,pan,reset,save,wheel_zoom",
               x_range=[0, 4*pi], y_range=[-1, 3])

plot2.line('x', 'u', source=source, line_width=3, line_alpha=0.6)


# Set up widgets
e_slider = Slider(title="e", value=0.0, start=0, end=1.0, step=0.01)
i_slider = Slider(title="i", value=0.0, start=0, end=180, step=1)
omega_slider = Slider(title="omega", value=0.0, start=0, end=180, step=1)
lp_slider = Slider(title="lambda_p", value=0.0, start=0, end=180, step=1)


# Set up callbacks
plot.title.text = 'Stokes Q'
plot2.title.text = 'Stokes U'


def update_data(attrname, old, new):

    # Get the current slider values
    e = e_slider.value
    i = i_slider.value*rad
    omega = omega_slider.value*rad
    lp = lp_slider.value*rad

    # Generate the new curve
    q = calc_q(l, e=e, i=i, omega=omega, lp=lp)
    u = calc_u(l, e=e, i=i, omega=omega, lp=lp)

    source.data = dict(x=l, q=q, u=u)


for w in [e_slider, i_slider, omega_slider, lp_slider]:
    w.on_change('value', update_data)


# Set up layouts and add to document
inputs = column(e_slider, i_slider, omega_slider, lp_slider)

curdoc().add_root(row(plot, plot2, inputs,  width=1200))
curdoc().title = "Stokes Parameters"
