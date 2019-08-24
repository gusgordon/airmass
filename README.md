# airmass
Python module for computing atmospheric airmass and column density, using [Reed Meyer's airmass program](http://reed.gigacorp.net/vitdownld.html#airmass). I have converted his C program to a Python extension.

For a detailed description of the model used, see [Reed's description](https://github.com/gusgordon/airmass/blob/master/extensions/airmassc.c).

## Installation
```
pip install https://github.com/gusgordon/airmass/zipball/master
```

## Usage
```
import airmass

# Specifying pressure and temperature
# relative_humidity_percent is an optional parameter. If specified, water vapor will be considered.
airmass.compute(zenith_degrees, altitude_meters, day_of_year, latitude, light_wavelength_angstroms, pressure_mb, temperature_c)

# Using an altitude
# Temperature and pressure are taken from experimental values at the given altitude
airmass.from_altitude(zenith_degrees, altitude_meters, day_of_year, latitude, light_wavelength_angstroms)
```
See the documentation of each function for more details.

## Solar
There is additionally a `solar` module, which is a work in progress. The goal of this is for accurate computation of solar intensity at various altitudes, so that high-altitude solar power output can be modeled ([for use in solar planes](https://github.com/gusgordon/atmosat)).
