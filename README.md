# airmass
Python module for computing atmospheric airmass and column density, using [Reed Meyer's airmass program](http://reed.gigacorp.net/vitdownld.html#airmass). I have converted his C program to a Python extension. This is a very accurate calculation - to quote Reed, it "was developed to find out to just what degree the standard approximate formulas for airmass are inaccurate".

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
There is additionally a `solar` module, which uses extinction coefficients at each wavelength to compute the solar irradiance under different conditions (namely at different altitudes and times) on Earth. Extinction from dry air, ozone, and water vapor are included. The goal of this is for accurate computation of solar intensity at various altitudes, so that high-altitude solar power output can be modeled ([for use in solar planes](https://github.com/gusgordon/atmosat)).

### Solar usage
```
airmass.solar_intensity_time(
    altitude,
    day_of_year,
    latitude,
    hour_of_day,
    include_diffuse_sky=True
)
```
