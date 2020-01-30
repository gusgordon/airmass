import numpy as np
import pandas as pd
from . import airmass as airmass_lib


def solar_intensity(
    altitude,
    day_of_year,
    latitude,
    wavelength=5000,
    hour_of_day=None,
    include_albedo=False,
):
    """
    Experimental - outputs solar intensities for a given wavelength of light.

    Todo:
     - Don't use airmass, but use column density with the extinction coefficients of atmospheric
       components (dry air, ozone, dust, and water vapor). Then compute intensity using
       Reed's exctinction formula.
    """
    albedo_zero_zenith = 0.25
    solar_irradiance_space = 1367.0  # W/m^2
    earth_radius = 6375000  # m

    # Declination angle from http://www.pveducation.org/pvcdrom/properties-of-sunlight/declination-angle
    declination = -23.45 * np.cos(360 / 365 * (day_of_year + 10) * np.pi / 180)

    # Loop through day's hours
    if hour_of_day is None:
        hour_angles = np.arange(0, 360)
    else:
        hour_angles = [int(np.rint((hour_of_day * 360) / 24)) + 180]

    intensities = []
    for hour_angle in hour_angles:
        p1 = np.sin(latitude * np.pi / 180) * np.sin(declination * np.pi / 180)
        p2 = (
            np.cos(latitude * np.pi / 180)
            * np.cos(declination * np.pi / 180)
            * np.cos(hour_angle * np.pi / 180)
        )
        zenith = np.arccos(p1 + p2) * 180 / np.pi

        if zenith > 89.9:
            intensities.append(0)
            continue

        airmass = airmass_lib.get_airmass_data(
            zenith, altitude, day_of_year, latitude, wavelength
        )
        airmass_sea_level_zenith = airmass_lib.get_airmass_data(
            0, 0, day_of_year, latitude, wavelength
        )
        relative_airmass = airmass / airmass_sea_level_zenith

        # Modified from http://www.pveducation.org/pvcdrom/2-properties-sunlight/air-mass
        intensity = (
            (1 + min(relative_airmass, 1) / 10)
            * 1353
            * 0.7 ** (relative_airmass ** 0.678)
            * np.cos(zenith * np.pi / 180)
        )

        # Albedo model: http://slideplayer.com/slide/10671870/
        if include_albedo:
            albedo = (
                albedo_zero_zenith
                + 4.9115e-9 * zenith ** 4
                + 6.0372e-8 * zenith ** 3
                - 2.1793e-5 * zenith ** 2
                + 1.3798e-3 * zenith
            )
            local_ff = (
                earth_radius / (earth_radius + altitude)
            ) ** 2  # Local form factor to Earth
            albedo_intensity = (
                solar_irradiance_space
                * albedo
                * local_ff
                * np.cos(zenith * np.pi / 180)
            )
            intensity += albedo_intensity

        intensities.append(intensity)

    return np.mean(np.concatenate([intensities]))


def yearly_solar_intensity(altitude, latitude=20, wavelength=5000, cloud_cover=None):
    mis = []
    for day in np.arange(0, 359, 30):
        mi = solar_intensity(
            altitude, day_of_year=day, latitude=latitude, wavelength=wavelength
        )
        mis.append(mi)
    mis_mean = np.mean(mis)
    if cloud_cover:
        attn = 0.135  # Clouds only let 13.5% of light through on average
        mis_mean = cloud_cover * mis_mean * attn + (1 - cloud_cover) * mis_mean * 1
    return mis_mean, mis_mean * 60 * 60 * 24 / (1000 * 3600)
