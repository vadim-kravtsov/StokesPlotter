# StokesPlotter
<b>StokesPlotter</b> is a program which plots the polarization, produced by the Thomson scattering of stellar radiation off the orbiting cloud of electrons in a binary system as a function of the orbital phase. The polarizaiton in that case depends on the orbital parameters of the binary star, as well as on the orientation of this orbit on the sky (see e.g. [Brown, McLean, 1977](https://ui.adsabs.harvard.edu/link_gateway/1977A%26A....57..141B/ADS_PDF) or [Kravtsov et al., 2020](https://doi.org/10.1051/0004-6361/202038745)).  

To run the program, install Bokeh and run the local server: 

```
bokeh serve precessing_disk.py
```

Then open the generated link in your favourite browser.
You will see the following interface, in which you can play with the orbital parameters to see how the polarization properties depend on them. 

![screenshot](https://github.com/vadim-kravtsov/StokesPlotter/blob/master/screenshot.png?raw=true)
