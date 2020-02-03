import os
import numpy as np
import pandas as pd
from airmassc import c_compute_airmass as _compute_airmass

# Atmospheric data from http://www.digitaldutch.com/atmoscalc/tableatmosphere.htm
dir_path = os.path.dirname(os.path.realpath(__file__))
atm_data = pd.read_csv(os.path.join(dir_path, "st_atm_1976.csv")).set_index("Altitude [m]")
atmosphere = (
    atm_data.reindex(np.arange(atm_data.index[0], atm_data.index[-1] + 1))
    .interpolate()
    .values
)


def compute(
    zenith_degrees,
    altitude,
    day_of_year,
    latitude,
    light_wavelength,
    P_Pa,
    T_K,
    relative_humidity=None,
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
    relative_humidity : float
        Additionally compute an estimate of the water vapor column density
        and airmass, using a crude model for the water vapor profile based
        on the local relative humidity at ground level in percent (100 = saturated).

    Returns
    -------
    dict
        Contains airmass and column densities, for dry air by default and for water
        vapor if relative_humidity is specified.
    """
    if relative_humidity is None:
        res = _compute_airmass(
            zenith_degrees, altitude, day_of_year, latitude, light_wavelength, P_Pa, T_c
        )
    else:
        res = _compute_airmass(
            zenith_degrees,
            altitude,
            day_of_year,
            latitude,
            light_wavelength,
            P_Pa,
            T_K,
            relhumid=relative_humidity,
        )

    if res[0] == "zenith_degrees not valid":
        raise ValueError(
            f"Zenith angle of {zenith_degrees} degrees is not valid. The zenith angle must be between 0 and 90 degrees."
        )
    elif res[0] == "water vapor":
        return {
            "type": "water_vapor_airmass",
            "water_vapor_column_density": res[1],
            "water_vapor_airmass": res[2],
        }
    elif res[0] == "regular":
        return {
            "type": "standard_airmass",
            "column_density": res[1],
            "zenith_column_density": res[2],
            "airmass": res[3],
        }
    else:
        raise Exception("Internal exception occurred. Try another set of values.")


def from_altitude(
    zenith_degrees, altitude, day_of_year, latitude, light_wavelength, relative_humidity=None
):
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
    relative_humidity : float
        Additionally compute an estimate of the water vapor column density
        and airmass, using a crude model for the water vapor profile based
        on the local relative humidity at ground level in percent (100 = saturated).

    Returns
    -------
    dict
        Contains airmass and column densities, for dry air by default and for water
        vapor if relative_humidity is specified.
    """
    T_K, P_Pa, _, _, _ = atmosphere[int(altitude)]

    return compute(
        zenith_degrees,
        altitude,
        day_of_year,
        latitude,
        light_wavelength,
        P_Pa,
        T_K,
        relative_humidity=relative_humidity,
    )
