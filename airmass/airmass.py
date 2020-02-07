import os
import numpy as np
import pandas as pd
from airmassc import c_compute_airmass as _compute_airmass

dir_path = os.path.dirname(os.path.realpath(__file__))

def compute(
    zenith_degrees, altitude, day_of_year, latitude, light_wavelength, P_Pa, T_K
):
    """Computes airmass and column density on Earth for a
    given wavelength of light.

    Uses a CPython version of the airmass program from Reed Meyer:
    http://reed.gigacorp.net/vitdownld.html#airmass

    Parameters
    ----------
    zenith_degrees : float
        The apparent zenith angle (not the true angle in the
        absence of atmosphere, but angle including refraction) at which
        you want the airmass.
    altitude : float
        The altitude of the observer above sea level, in meters.
        Consider using the ``from_altitude`` function below if you are using
        standard atmospheric conditions, and would like temperature and pressure
        automatically applied for your supplied altitude.
    day_of_year : float
        The date, in days since the beginning of the year
        (0 = midnight Jan. 1; the default is 80, corresponding roughly to
        the vernal equinox); don't worry about fractions of a day, since the
        seasonal effect on airmass is relatively small.
    latitude : float
        The observer's latitude on Earth.
    light_wavelength : float
        The observing wavelength in Angstroms. Don't worry too much about an
        exact figure here, as the wavelength dependence is very small.
    P_Pa : float
        Specifies the local pressure in Pascals.
        Consider using the ``from_altitude`` function below if you are using
        standard atmospheric conditions, and would like temperature and pressure
        automatically applied for your supplied altitude.
    T_K : float
        Specifies the local temperature in Kelvin.
        Consider using the ``from_altitude`` function below if you are using
        standard atmospheric conditions, and would like temperature and pressure
        automatically applied for your supplied altitude.

    Returns
    -------
    dict
        Contains airmass and column densities for dry air.
    """

    res = _compute_airmass(
        zenith_degrees, altitude, day_of_year, latitude, light_wavelength, P_Pa, T_K
    )

    if res[0] == "zenith_degrees not valid":
        raise ValueError(
            f"Zenith angle of {zenith_degrees} degrees is not valid. The zenith angle must be between 0 and 90 degrees."
        )
    elif res[0] == "regular":
        return {
            "type": "standard_airmass",
            "column_density": res[1],
            "zenith_column_density": res[2],
            "airmass": res[3],
        }
    else:
        raise Exception("Internal exception occurred. Try another set of values.")


def water_vapor_col_density(
    zenith_degrees,
    altitude,
    day_of_year,
    latitude,
    light_wavelength,
    relative_humidity_sea_level,
):
    """Computes a rough approximation for the column density of water vapor.

    Reed's program has two issues with the humidity calculations:
      1. It uses a relative humidity at altitude, whereas most relative humidity measurements
      are at sea level. This version allow a relative humidity at sea level
      to be specified.
      2. Hangs on this calculation for altitudes greater than roughly
      4000 m, so this functions attempts to "patch" that by computing the column
      density at sea level then extrapolating to higher altitudes.

    Since the original water vapor model is simple, this is relatively simple
    to approximate roughly.

    Parameters
    ----------
    zenith_degrees : float
        The apparent zenith angle (not the true angle in the
        absence of atmosphere, but angle including refraction) at which
        you want the airmass.
    altitude : float
        The altitude of the observer above sea level, in meters.
        Consider using the ``from_altitude`` function below if you are using
        standard atmospheric conditions, and would like temperature and pressure
        automatically applied for your supplied altitude.
    day_of_year : float
        The date, in days since the beginning of the year
        (0 = midnight Jan. 1; the default is 80, corresponding roughly to
        the vernal equinox); don't worry about fractions of a day, since the
        seasonal effect on airmass is relatively small.
    latitude : float
        The observer's latitude on Earth.
    light_wavelength : float
        The observing wavelength in Angstroms. Don't worry too much about an
        exact figure here, as the wavelength dependence is very small.
    P_Pa : float
        Specifies the local pressure in Pascals.
        Consider using the ``from_altitude`` function below if you are using
        standard atmospheric conditions, and would like temperature and pressure
        automatically applied for your supplied altitude.
    T_K : float
        Specifies the local temperature in Kelvin.
        Consider using the ``from_altitude`` function below if you are using
        standard atmospheric conditions, and would like temperature and pressure
        automatically applied for your supplied altitude.

    Returns
    -------
    float
        Water vapor column density.
    """

    troposphere_top_m = 17000
    lower_stratosphere_top_m = 30000

    if altitude < troposphere_top_m:
        rel_humid = relative_humidity_sea_level
    elif altitude < lower_stratosphere_top_m:
        rel_humid = relative_humidity_sea_level * 0.5
    else:
        return 0.0

    water_vapor_col_density_0_alt = (
        _compute_airmass(
            zenith_degrees,
            0,
            day_of_year,
            latitude,
            light_wavelength,
            101325,
            288.15,
            rel_humid,
        )[1]
        * -0.5
    )

    mean_humid_alt_falloff = (troposphere_top_m + lower_stratosphere_top_m) / 2

    water_vapor_col_density = max(
        (1 - altitude / mean_humid_alt_falloff) * water_vapor_col_density_0_alt, 0
    )

    return water_vapor_col_density


def from_altitude(zenith_degrees, altitude, day_of_year, latitude, light_wavelength):
    """Computes airmass and column density on Earth for a
    given altitude. Uses the 1976 atmospheric data to compute temperature and
    pressure for the given altitude. If you would like to specify temperature and
    pressure yourself, use the compute function above.

    Parameters
    ----------
    zenith_degrees : float
        The apparent zenith angle (not the true angle in the
        absence of atmosphere, but angle including refraction) at which
        you want the airmass.
    altitude : float
        The altitude of the observer above sea level, in meters.
    day_of_year : float
        The date, in days since the beginning of the year
        (0 = midnight Jan. 1; the default is 80, corresponding roughly to
        the vernal equinox); don't worry about fractions of a day, since the
        seasonal effect on airmass is relatively small.
    latitude : float
        The observer's latitude on Earth.
    light_wavelength : float
        The observing wavelength in Angstroms. Don't worry too much about an
        exact figure here, as the wavelength dependence is very small.

    Returns
    -------
    dict
        Contains airmass and column densities for dry air.
    """

    # Atmospheric data from http://www.digitaldutch.com/atmoscalc/tableatmosphere.htm
    atm_data = pd.read_csv(os.path.join(dir_path, "data/st_atm_1976.csv")).set_index("Altitude [m]")
    atmosphere = (
        atm_data.reindex(np.arange(atm_data.index[0], atm_data.index[-1] + 1))
        .interpolate()
        .values
    )

    T_K, P_Pa, _, _, _ = atmosphere[int(altitude)]

    return compute(
        zenith_degrees,
        altitude,
        day_of_year,
        latitude,
        light_wavelength,
        P_Pa,
        T_K,
    )
