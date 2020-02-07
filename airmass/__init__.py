from airmassc import c_compute_airmass as _compute_airmass
from .airmass import compute, from_altitude, water_vapor_col_density
from .solar import (
    get_solar_spectrum_modtran, solar_intensity, yearly_solar_intensity,
    get_ozone_ppm, get_cross_section, get_air_mass_extinction,
    get_ozone_particle_col_density, get_cloud_attentuation,
    get_solar_irradiance,
)
