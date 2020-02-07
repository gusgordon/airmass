import numpy as np
import pandas as pd
import airmass as airmass_lib

N_Avogadro = 6.02214076e23

atmosphere_composition_by_mass = {
    "N2": 0.75519454,
    "O2": 0.23140128,
    "Ar": 0.01288169,
    "CO2": 0.0004771,
    "SO2": 0.00000221,
    "CH4": 0.00000111,
}

molar_mass_g_per_mol = {
    "N2": 28.0134,
    "O2": 31.9988,
    "Ar": 39.948,
    "CO2": 44.0095,
    "SO2": 64.0638,
    "CH4": 16.04246,
}

air_molar_mass = 0.0289647  # kg/mol


def get_solar_spectrum_modtran():
    """
    Load the solar spectrum from MODTRAN 6.

    Returns
    -------
    pd.DataFrame
        Solar spectrum in increments of 0.1 nm, containing
        irradiance and power contribution for the wavelength.
    """
    # Load in MODTRAN extraterrestrial specta
    # Source: https://www.nrel.gov/grid/solar-resource/spectra.html
    solar_spectrum_all = (
        pd.read_csv("data/MODTRAN-W_M-2_nm-1.csv")
        .drop("(CM-1)", axis=1)
        .set_index("nm")
    )

    # For each wavelength, take the mean across all samples
    solar_spectrum_mean = solar_spectrum_all.mean(axis=1)

    # Reduce the precision to something more reasonable
    reduced_solar_spectrum = solar_spectrum_mean.groupby(
        solar_spectrum_mean.index.to_series().round(1)
    ).mean()

    # Also compute the power contribution for each wavelength
    # This gives a high estimate, as there is imprecision on the high end of the spectrum
    # For cases involving solar cells, the high end of the spectrum is cut out anyway,
    # so it doesn't matter here.
    reduced_solar_spectrum = pd.DataFrame(
        {"irradiance": reduced_solar_spectrum}
    ).reset_index()
    reduced_solar_spectrum["wavelength_width"] = reduced_solar_spectrum.nm.diff()
    reduced_solar_spectrum.loc[0, "wavelength_width"] = 0.1
    reduced_solar_spectrum["power_contribution"] = (
        reduced_solar_spectrum.irradiance * reduced_solar_spectrum.wavelength_width
    )

    return reduced_solar_spectrum.set_index("nm")


def get_ozone_ppm():
    """
    Load the ozone concentration at different altitudes.

    Returns
    -------
    pd.DataFrame
        Solar spectrum in increments of 0.1 nm, containing
        irradiance and power contribution for the wavelength.
    """

    # Load in NASA ozone vs. altitude data
    # Source: https://ozonewatch.gsfc.nasa.gov/facts/SH.html
    ozone_ppm_vs_altitude = pd.read_csv(
        "data/ozone_vs_altitude.csv", index_col="altitude_km"
    ).ozone_ppm.sort_index()

    ozone_frac_vs_altitude = ozone_ppm_vs_altitude * 1e-6
    ozone_frac_vs_altitude.index = ozone_frac_vs_altitude.index * 1e3
    ozone_frac_vs_altitude.index.name = "altitude_m"

    new_index = np.concatenate(
        [
            ozone_frac_vs_altitude.index,
            np.arange(0, ozone_frac_vs_altitude.index[-1] + 1),
        ]
    )

    # Interpolate ozone for all altitudes
    ozone_frac_vs_altitude = (
        (ozone_frac_vs_altitude.reindex(new_index))
        .sort_index()
        .interpolate(method="cubic")
    )

    # Fill low altitude values as 0 concentration
    ozone_frac_vs_altitude = ozone_frac_vs_altitude.fillna(0)

    # Round the altitude to the nearest meter so it's easier to query
    ozone_frac_vs_altitude.index = ozone_frac_vs_altitude.index.to_series().round()

    return ozone_frac_vs_altitude


def get_cross_section(compound):
    """
    Load the extinction cross section (absorbtion + Rayleigh) for a
    given compound. Other sources of extinction are negligible.
    Returns cross sections in units of cm^2.

    This data was generated using the generate_extinction.ipynb
    notebook in the etc/ directory. The basis of this notebook
    is https://github.com/sukritranjan/ranjanwordsworthsasselov2016/blob/bd5de2cb081b129bfd6b3b5bc173ca4e9ce7702b/extract_gas_cross_sections.py.

    Supported compounds:
     - N2
     - Ar (absorbtion cross section is ignored, need to add this, but it's negligible)
     - O2
     - CO2
     - H2O
     - CH4
     - O3
     - SO2

    Returns
    -------
    pd.DataFrame
        Total extintion cross section in increments of 0.1 nm.
    """

    cross_section = pd.read_csv(
        f"data/composite_xc_extended_{compound.lower()}",
        delim_whitespace=True,
        index_col="wav_nm",
    )

    return cross_section


def get_air_mass_extinction():
    """
    Load the total extinction cross section (absorbtion + Rayleigh)
    for air in cm^2/g.

    Included compounds:
     - N2
     - Ar (absorbtion cross section is ignored, need to add this, but it's negligible)
     - O2
     - CO2
     - CH4
     - SO2

    Returns
    -------
    pd.Series
        Total extintion cross section of air in increments of 0.1 nm.
    """

    xs = {}
    for molecule, contribution in atmosphere_composition_by_mass.items():
        xs[molecule] = get_cross_section(molecule).total_cm2.mul(
            N_Avogadro / molar_mass_g_per_mol[molecule]
        )

        xs[molecule] = xs[molecule].mul(contribution)

        xs[molecule] = (
            xs[molecule].groupby(xs[molecule].index.to_series().round(1)).mean()
        )

    # Interpolate missing values linearly. This mostly fills in the CO2 spectrum
    xs = pd.DataFrame(xs)[200:900].interpolate()

    return xs.sum(axis=1)


def get_ozone_particle_col_density(zenith_degrees, altitude):
    """Computes the ozone particle column density in particles/m^2.

    Parameters
    ----------
    zenith_degrees : float
        The apparent zenith angle (not the true angle in the
        absence of atmosphere, but angle including refraction) at which
        you want the airmass.
    altitude : float
        The altitude of the observer above sea level, in meters.

    Returns
    -------
    float
        Ozone particle column density in particles/m^2.
    """

    atm_data = pd.read_csv("data/st_atm_1976.csv").set_index("Altitude [m]")
    atmosphere = atm_data.reindex(
        np.arange(atm_data.index[0], atm_data.index[-1] + 1)
    ).interpolate()

    n_air_particles = atmosphere["Density [kg/m3]"] / air_molar_mass * N_Avogadro

    ozone = pd.DataFrame({"ppm": get_ozone_ppm(), "n_air_particles": n_air_particles})
    ozone["rel_ppm"] = ozone.ppm.div(np.cos(zenith_degrees * np.pi / 180))
    ozone = ozone[:50000]
    ozone["ozone_particles"] = ozone.n_air_particles * ozone.rel_ppm / 1e6

    return ozone[altitude:].ozone_particles.sum()


def get_cloud_attentuation(altitude, cloud_level_m=9000):
    """Computes the cloud attenuation using a simple model of a single
    cloud at the given altitude.
    
    Note this is effectively the worst case scenario for cloud cover.

    Parameters
    ----------
    altitude : float
        The altitude of the observer above sea level, in meters.

    Returns
    -------
    float
        Cloud attenuation.
    """

    if altitude > cloud_level_m:
        return 0
    else:
        return 0.85


def get_solar_irradiance(
    zenith_degrees, altitude, day_of_year, latitude, include_diffuse_sky=True
):
    """Get the solar irradiance in watts/m^2 at a certain altitude and given certain conditions.=
    Utilizes exctinction coefficients across 200-900 nm wavelengths for dry air, water vapor, ozone,
    and clouds (using a simple model for clouds).
    
    Wavelengths outside of 200-900 nm are extrapolated, however this is roughly valid given that
    Rayleigh scattering is the primary contributor to extinction and it is fairly consistent at
    the larger wavelengths.
    
    A way to improve this function would be to incorporate extinction coefficients at larger
    wavelengths.
    
    Outputs power in watts/m^2.

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
    include_diffuse_sky : boolean, default True
        Approximately add in the contribution from diffuse sky irradiance.
        This is a rough model, but regardless, the contribution is fairly small.

    Returns
    -------
    float
        Solar irradiance in watts/m^2.
    """

    # Dry air
    cd = {}
    for wavelength in np.arange(200, 900.1, 0.1):
        if np.round(wavelength, 1) % 10 == 0:
            cd[wavelength] = airmass_lib.from_altitude(
                zenith_degrees, altitude, day_of_year, latitude, wavelength * 10
            )["column_density"]
        else:
            cd[wavelength] = np.nan

    wavelength_col_densities = pd.Series(cd).interpolate()

    wavelength_col_densities.index = wavelength_col_densities.index.to_series().round(1)

    air_extinction_coefs = wavelength_col_densities.mul(get_air_mass_extinction())

    # Ozone
    o3_cross_section = get_cross_section("O3").total_cm2
    o3_cross_section = o3_cross_section.reindex(
        np.concatenate([o3_cross_section.index, np.arange(200.0, 900.1, 0.1)])
    ).sort_index()
    o3_cross_section = o3_cross_section.interpolate()
    o3_cross_section = o3_cross_section.groupby(
        o3_cross_section.index.to_series().round(1)
    ).mean()

    ozone_extinction_coefs = o3_cross_section[200:900] * get_ozone_particle_col_density(
        zenith_degrees, altitude
    )

    # Water vapor
    wv_col_density = airmass_lib.water_vapor_col_density(
        zenith_degrees, altitude, day_of_year, latitude, 5500, 0.5
    )

    h2o_cross_section = get_cross_section("h2o").total_cm2[200:900]
    h2o_cross_section = h2o_cross_section.reindex(
        np.concatenate([h2o_cross_section.index, np.arange(200.0, 900.1, 0.1)])
    ).sort_index()
    h2o_cross_section = h2o_cross_section.interpolate()
    h2o_cross_section = h2o_cross_section.groupby(
        h2o_cross_section.index.to_series().round(1)
    ).mean()
    h2o_extinction_coefs = h2o_cross_section.mul(N_Avogadro / 18.01528)

    # Clouds
    cloud_extinction = get_cloud_attentuation(altitude)

    # Sum the extinctions to get the total extinction
    total_extinction = (
        air_extinction_coefs
        + ozone_extinction_coefs
        + h2o_extinction_coefs
        + cloud_extinction
    )

    solar_spectrum_power = get_solar_spectrum_modtran().power_contribution

    # Get the total power *in the 200-900 nm range*
    visible_power_spectrum = solar_spectrum_power * np.power(np.e, -total_extinction)
    visible_power_total = visible_power_spectrum.sum()

    # Since we only had exctinction coefficients in the visible spectrum, let's multiply this back
    # out as an approximation for the total spectrum. This should be roughly accurate, since most
    # radiation is from Rayleigh scattering which is roughly constant for larger wavelengths.
    frac_visible_power = (
        solar_spectrum_power[200:900].sum() / solar_spectrum_power.sum()
    )
    total_power_approx = visible_power_total / frac_visible_power

    # We use a very simple model for the diffuse sky irradiance
    # since the diffuse sky irradiance is fairly insignificant.
    # We assume that at sea level, there is an 8% contribution from diffuse
    # sky, and this drops off linearly up until 20,000 m, above which the
    # contribution is assumed to be zero.
    diffuse_sky_factor = max(0.08 * min(1, (20000 - altitude) / 20000), 0)

    # We also multiply by the cosine of the zenith angle to account for the projection
    # of sunlight onto a flat surface
    total_power_incedent = (
        total_power_approx
        * (1 + diffuse_sky_factor)
        * np.cos(zenith_degrees * np.pi / 180)
    )

    return total_power_incedent

def solar_intensity_time(
    altitude,
    day_of_year,
    latitude,
    hour_of_day,
    include_diffuse_sky=True,
):
    """Get the solar irradiance in watts/m^2 at a certain time.
    Outputs power in watts/m^2. See doc for get_solar_irradiance above
    for more details.

    Parameters
    ----------
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
    hour_of_day : float
        The current hour, for example 7.5 for 7:30 AM.
    include_diffuse_sky : boolean, default True
        Approximately add in the contribution from diffuse sky irradiance.
        This is a rough model, but regardless, the contribution is fairly small.

    Returns
    -------
    float
        Solar irradiance in watts/m^2.
    """

    # Declination angle from http://www.pveducation.org/pvcdrom/properties-of-sunlight/declination-angle
    declination = -23.45 * np.cos(360 / 365 * (day_of_year + 10) * np.pi / 180)

    hour_angle = int(np.rint((hour_of_day * 360) / 24)) + 180

    p1 = np.sin(latitude * np.pi / 180) * np.sin(declination * np.pi / 180)
    p2 = (
        np.cos(latitude * np.pi / 180)
        * np.cos(declination * np.pi / 180)
        * np.cos(hour_angle * np.pi / 180)
    )
    zenith = np.arccos(p1 + p2) * 180 / np.pi

    if zenith > 89.9:
        return 0.0

    return get_solar_irradiance(
        zenith, altitude, day_of_year, latitude,
        include_diffuse_sky=include_diffuse_sky
    )

    return intensity
