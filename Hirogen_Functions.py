import sys
import os
import math
import spectres

import numpy as np
import scipy.constants as constants
import mysql.connector

from astropy.io import fits
from astropy import units as u
from astropy.convolution import convolve, Box1DKernel
from scipy.signal import medfilt, savgol_filter
from scipy.ndimage import gaussian_filter1d
from dust_extinction.parameter_averages import F99
from numpy.fft import *
from mysql.connector import errorcode

np.set_printoptions(threshold=sys.maxsize)
c = constants.value('speed of light in vacuum') / 1000

Debug = False
NERSC = False

# Imports NERSC specific modules or common modules on different paths
if NERSC:
    # DESI Specific Packages

    import fitsio
    # import healpy as hp
    # When last I checked there was an issue with importing healpy on NERSC due to a circular import issue

    import desispec.io

    from glob import glob
    from collections import defaultdict

    # ADM Note that we use the commissioning targeting mask, as we're working with mini-SV data from commissioning.
    from desitarget.cmx.cmx_targetmask import cmx_mask


def user_config(user):
    if user == 'Peter':
        database_name = 'CoronalLineGalaxies'
        database_user = 'root'
        database_password = 'Th0thArchive'
        # Likely not the best idea to have the password essentially public but its not like its actually sensitive info

        # Used only in SDSS mode
        # main_spectra_path = '/Volumes/GoogleDrive/My Drive/TDE_Project/Coronal_Line_Galaxies/Existing_Sample_Data/SDSS_Spectra_DR16'
        main_spectra_path = f"/Volumes/GoogleDrive/My Drive/TDE_Project/Coronal_Line_Galaxies"
        # main_spectra_path = f"/Volumes/GoogleDrive/My Drive/TDE_Project/Changing_Look_QSOs/Followup/dr16/eboss/spectro/redux/v5_13_0/spectra/lite"

    elif user == 'Joe':
        # Add your paths + details here when you have them set up!
        database_name = 'CoronalLineGalaxies'  # Not the specific database table name
        database_user = 'root'
        database_password = 'C5bhRxH4TKzZ9bhSgAfd'
        main_spectra_path = '.'  # This should be the folder than contains the dr16 folder from the download script

    else:
        print('User not recognised - please configure your details and paths in the Hirogen_Functions.py file.')
        sys.exit()

    return database_name, database_user, database_password, main_spectra_path


def main_config():
    """ This function defines criteria used by all parts of Hirogen to ensure consistency"""

    # These four parameters control what regions of the spectra are 'cut out' for analysis around each line
    lower_wave_region = 125
    upper_wave_region = 125

    lower_v_shift = -3000
    upper_v_shift = 3000

    # Any object with a candidate score equal to or exceeding this value will be explored by the detailed line fitter
    detailed_candidate_threshold = 5

    # Any object with a candidate score equal to or exceeding this value will be treated as an ECLE candidate
    candidate_threshold = 7

    # Maximum possible candidate score
    candidate_score_max = 11

    # To qualify as a line detection for candidate selection the overall feature pEQW must exceed this value:
    line_detection_peqw_threshold = -2

    # To qualify as a line detection for candidate selection the maximum flux point of the feature must exceed
    # this value
    line_detection_max_threshold = 1.05

    line_detection_max_above_av_continua_threshold = 1.05

    # To qualify as a line detection for candidate selection the maximum flux point of the feature must occur
    # within this +_ kms^-1 of the zero velocity point
    line_detection_peak_tolerance = 350

    # The region inside which the peak max is considered - a peak outside this region will not be considered to be the
    # true line peak. Helps filter out where lines are closely spaced with potentially unrelated features (kms^-1)
    line_peak_location_region = 750

    # If the minimum of the line peak region falls below this value the line cannot score - designed to weed out very
    # noisy lines

    # line_peak_region_min_threshold = 0.93
    line_peak_region_min_threshold = -0.3  # This is for DESI

    # Thresholds for identifying strong lines. First is overall line strength, second is for detecting sharp peaks
    strong_eqw_threshold = -7.5
    strong_peak_max_threshold = 2

    return lower_wave_region, upper_wave_region, lower_v_shift, upper_v_shift, candidate_threshold, \
           candidate_score_max, line_detection_peqw_threshold, line_detection_max_threshold, \
           line_detection_max_above_av_continua_threshold, line_detection_peak_tolerance, line_peak_location_region, \
           line_peak_region_min_threshold, strong_eqw_threshold, strong_peak_max_threshold, \
           detailed_candidate_threshold


""" Makes the main_config function accessible to the other functions in this file by default"""
Main_Config_Data = main_config()
R_ext = Main_Config_Data[0]

""" Line Lists"""


def lines_for_analysis():
    """Serves as the internal line list in and information storage location.

    All lines included have the EQWs measured. 'Air' as these are line locations at STP"""

    # The [OII] line is actually a doublet of 3727.320 and 3729.225

    spectral_lines_for_analysis_air = {
        'Halpha': [6562.79, 'indigo', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        'Hbeta': [4861.35, 'indigo', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        'Hgamma': [4340.472, 'indigo', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        'Hdelta': [4101.734, 'indigo', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],

        '[NII]6548': [6548, 'deeppink', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        '[NII]6584': [6584, 'deeppink', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],

        '[FeVII]3759': [3759, 'orangered', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        '[FeVII]5160': [5160, 'orangered', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        '[FeVII]5722': [5722, 'orangered', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        '[FeVII]6088': [6088, 'orangered', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        '[FeX]6376': [6376, 'orangered', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        '[FeXI]7894': [7894, 'orangered', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        '[FeXIV]5304': [5304, 'orangered', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],

        '[OI]6300': [6300, 'mediumaquamarine', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        '[OI]6363': [6363, 'mediumaquamarine', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        '[OII]3728': [3728.2725, 'mediumaquamarine', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        '[OIII]4959': [4959, 'mediumaquamarine', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        '[OIII]5007': [5007, 'mediumaquamarine', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],

        'HeI4478': [4478, 'royalblue', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        'HeII4686': [4686, 'royalblue', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],

        'NaID': [5892.935, 'orange', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],

        '[SII]6716': [6716, 'goldenrod', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],
        '[SII]6731': [6731, 'goldenrod', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],

        '[SII]6716_6731': [6723.5, 'goldenrod', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],

        '[SIII]6313': [6313, 'goldenrod', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],

        '[NeIII]3896': [3896, 'blueviolet', [], [], [], [], [], [], [], [], [], [], [], [], [], [], []],

    }

    return spectral_lines_for_analysis_air


def lines_for_scoring():
    """The lines of that factor into the candidate score determination"""
    spectral_lines_for_scoring_air = [
        'Halpha',
        '[FeVII]6088',
        '[FeX]6376',
        '[FeXI]7894',
        '[FeXIV]5304',
        '[OIII]5007',
        'HeII4686'
    ]

    return spectral_lines_for_scoring_air


def coronal_lines():
    """Fe lines"""
    spectral_coronal_lines_air = [
        '[FeVII]6088',
        '[FeX]6376',
        '[FeXI]7894',
        '[FeXIV]5304'
    ]

    return spectral_coronal_lines_air


def primary_lines():
    """The lines of primary interest and those which are included in the main snapshot plot"""

    lines_primary = [

        '[FeVII]6088',
        '[FeX]6376',
        '[FeXI]7894',
        '[FeXIV]5304',

        'Halpha',
        'Hdelta',

        'HeII4686',

        '[OIII]5007',

    ]

    return lines_primary


def primary_line_labels():
    """Line labels for the primary lines of interest formatted in a way to properly display them on plots"""
    line_labels_primary = [
        r'[FeVII] $\lambda$ 6088 $\AA$',
        r'[FeX] $\lambda$ 6376 $\AA$',
        r'[FeXI] $\lambda$ 7894 $\AA$',
        r'[FeXIV] $\lambda$ 5304 $\AA$',

        r'H$\alpha$',
        r'H$\delta$',

        r'HeII $\lambda$ 4686 $\AA$',

        r'[OIII] $\lambda$ 5007 $\AA$',
    ]
    return line_labels_primary


def primary_line_labels_short():
    """Shortened Line labels for the primary lines for when space is an issue"""
    short_line_labels_primary = [
        r'[FeVII]',
        r'[FeX]',
        r'[FeXI]',
        r'[FeXIV]',

        r'H$\alpha$',
        r'H$\delta$',

        r'HeII',

        r'[OIII]',
    ]
    return short_line_labels_primary


def secondary_lines():
    """The lines of secondary interest: Balmer series, He features, the [NII] lines around Halpha, [OIII] and [SII]"""
    lines_secondary = [

        'Halpha',
        'Hbeta',
        'Hgamma',
        'Hdelta',

        'HeI4478',
        'HeII4686',

        '[NII]6548',
        '[NII]6584',

        '[OI]6300',
        '[OI]6363',
        '[OII]3728',
        '[OIII]4959',
        '[OIII]5007',

        '[SII]6716',
        '[SII]6731',
        '[SII]6716_6731',
        '[SIII]6313'
    ]

    return lines_secondary


def secondary_line_labels():
    """Long form line labels for the secondary lines"""
    line_labels_secondary = [
        r'H$\alpha$',
        r'H$\beta$',
        r'H$\gamma$',
        r'H$\delta$',

        r'HeI $\lambda$ 4478 $\AA$',
        r'HeII $\lambda$ 4686 $\AA$',

        r'[NII] $\lambda$ 6548 $\AA$',
        r'[NII] $\lambda$ 6584 $\AA$',

        r'[OI] $\lambda$ 6300 $\AA$',
        r'[OI] $\lambda$ 6363 $\AA$',
        r'[OII] $\lambda$ 3728 $\AA$',
        r'[OIII] $\lambda$ 4959 $\AA$',
        r'[OIII] $\lambda$ 5007 $\AA$',

        r'[SII] $\lambda$ 6716 $\AA$',
        r'[SII] $\lambda$ 6731 $\AA$',
        r'[SII] $\lambda$$\lambda$ 6716, 6731 $\AA$',
        r'[SIII] $\lambda$ 6313 $\AA$'

    ]
    return line_labels_secondary


def tertiary_lines():
    """Lines of tertiary interest: Weaker Fe lines, NaID, [NeIII]"""
    lines_tertiary = [

        '[FeVII]3759',
        '[FeVII]5160',
        '[FeVII]5722',

        'NaID',

        '[NeIII]3896',

    ]

    return lines_tertiary


def tertiary_line_labels():
    """Long form line labels for the tertiary lines"""
    line_labels_tertiary = [
        r'[FeVII] $\lambda$ 3759 $\AA$',
        r'[FeVII] $\lambda$ 5160 $\AA$',
        r'[FeVII] $\lambda$ 5722 $\AA$',

        'NaID',

        r'[NeIII]$\lambda$ 3896 $\AA$',
    ]
    return line_labels_tertiary


def comparison_lines():
    """The lines of primary interest and those which are included in the main snapshot plot"""

    lines_comparison = [

        '[FeVII]3759',
        '[FeVII]5160',
        '[FeVII]5722',

        '[FeVII]6088',
        '[FeX]6376',
        '[FeXI]7894',
        '[FeXIV]5304',

        'Halpha',
        'Hdelta',

        'HeII4686',

        '[OIII]5007',

    ]

    return lines_comparison


def comparison_plot_line_labels_short():
    """Shortened Line labels for the comparison plot lines for when space is an issue"""
    short_line_labels_comparison = [
        r'[FeVII]',
        r'[FeVII]',
        r'[FeVII]',
        r'[FeVII]',
        r'[FeX]',
        r'[FeXI]',
        r'[FeXIV]',

        r'H$\alpha$',
        r'H$\delta$',

        r'HeII',

        r'[OIII]',
    ]
    return short_line_labels_comparison


def comparison_lines_extended():
    """Extended list of lines to include in the followup comparison plot """

    lines_comparison = [

        '[FeVII]3759',
        '[FeVII]5160',
        '[FeVII]5722',

        '[FeVII]6088',
        '[FeX]6376',
        '[FeXI]7894',
        '[FeXIV]5304',

        'Halpha',
        'Hbeta',
        'Hdelta',

        'HeI4478',
        'HeII4686',
    ]

    return lines_comparison


def comparison_plot_line_labels_short_extended():
    """Shortened Line labels for the followup comparison plot lines for when space is an issue"""
    short_line_labels_comparison = [
        r'[FeVII]',
        r'[FeVII]',
        r'[FeVII]',
        r'[FeVII]',
        r'[FeX]',
        r'[FeXI]',
        r'[FeXIV]',

        r'H$\alpha$',
        r'H$\beta$',
        r'H$\delta$',

        r'HeI',
        r'HeII',
    ]
    return short_line_labels_comparison


def presentation_zoom_lines():
    """The lines of primary interest and those which are included in the main snapshot plot"""

    lines_primary = [
        '[FeVII]6088',
        '[FeX]6376',
        '[FeXI]7894',
        '[FeXIV]5304',
    ]

    return lines_primary


def presentation_zoom_line_labels():
    """Line labels for the presentation lines of interest formatted in a way to properly display them on plots"""
    line_labels_presentation = [
        r'[FeVII] $\lambda$ 6088 $\AA$',
        r'[FeX] $\lambda$ 6376 $\AA$',
        r'[FeXI] $\lambda$ 7894 $\AA$',
        r'[FeXIV] $\lambda$ 5304 $\AA$',
    ]
    return line_labels_presentation


def interesting_lines():
    # More general line list - currently used in the template comparison plots

    lines_of_interest = [
        'Halpha',
        'Hbeta',
        'Hgamma',
        'Hdelta',

        'HeI4478',
        'HeII4686',

        '[OI]6300',
        '[OII]3728',
        '[OIII]4959',
        '[OIII]5007',

        '[FeVII]3759',
        '[FeVII]5160',
        '[FeVII]5722',
        '[FeVII]6088',
        '[FeX]6376',
        '[FeXI]7894',
        '[FeXIV]5304',

        '[SII]6716_6731'
    ]

    return lines_of_interest


def interesting_line_labels():
    line_labels_interesting = [

        r'H$\alpha$',
        r'H$\beta$',
        r'H$\gamma$',
        r'H$\delta$',

        r'HeI $\lambda$ 4478 $\AA$',
        r'HeII $\lambda$ 4686 $\AA$',

        r'[OI] $\lambda$ 6300 $\AA$',
        r'[OII] $\lambda$ 3728 $\AA$',
        r'[OIII] $\lambda$ 4959 $\AA$',
        r'[OIII] $\lambda$ 5007 $\AA$',

        r'[FeVII] $\lambda$ 3759 $\AA$',
        r'[FeVII] $\lambda$ 5160 $\AA$',
        r'[FeVII] $\lambda$ 5722 $\AA$',
        r'[FeVII] $\lambda$ 6088 $\AA$',
        r'[FeX] $\lambda$ 6376 $\AA$',
        r'[FeXI] $\lambda$ 7894 $\AA$',
        r'[FeXIV] $\lambda$ 5304 $\AA$',

        r'[SII] $\lambda$$\lambda$ 6716, 6731 $\AA$'
    ]
    return line_labels_interesting


def interesting_line_labels_short():
    line_labels_interesting = [

        r'H$\alpha$',
        r'H$\beta$',
        r'H$\gamma$',
        r'H$\delta$',

        r'HeI4478',
        r'HeII4686',

        r'[OI]6300',
        r'[OII]3728',
        r'[OIII]4959',
        r'[OIII]5007',

        r'[FeVII]3759',
        r'[FeVII]5160',
        r'[FeVII]5722',
        r'[FeVII]6088',
        r'[FeX]6376',
        r'[FeXI]7894',
        r'[FeXIV]5304',

        r'[SII]6716_6731'
    ]
    return line_labels_interesting


"""Spectroscopic analysis functions"""


def region_cutter(shift, wave, flux, low_cut, high_cut, mode='Shift', error=None):
    """
    Cuts down a spectrum, flux, wavelength and determined shift, based on either a shift or wavelength region
    Defaults to shift for compatibility
    Expects wave and flux to have a unit but does NOT preserve these
    """

    wave = wave.value
    flux = flux.value

    cut_shift = []
    cut_wave = []
    cut_flux = []
    cut_error = []

    if mode.upper() in {'SHIFT', 'VELOCITY'}:
        for xx, item in enumerate(shift):

            if low_cut <= shift[xx] <= high_cut:
                cut_shift.append(shift[xx])
                cut_wave.append(wave[xx])
                cut_flux.append(flux[xx])
                if error:
                    cut_error.append(error[xx])

            elif shift[xx] > high_cut:
                break

    elif mode.upper() in {'WAVELENGTH', 'LAMBDA'}:

        for xx, item in enumerate(wave):

            if low_cut <= wave[xx] <= high_cut:
                cut_shift.append(shift[xx])
                cut_wave.append(wave[xx])
                cut_flux.append(flux[xx])
                if error:
                    cut_error.append(error[xx])

            elif wave[xx] > high_cut:
                break

    else:
        print('Mode selection not recognised.\nPlease check and try again.')
        sys.exit()

    return cut_shift, cut_wave, cut_flux, cut_error


def continua_maker(spec_region, flux, line_name, line_loc, object_name):
    """This version of the continua generator and corrector uses a wavelength region rather than a relative velocity"""

    if line_name.lower() in {'halpha', '[nii]6548', '[nii]6584'}:

        # These are the values Decker used for her paper with Or.
        # They will include parts of broad features so I am tweaking.
        """
        blue_start = 6492.8
        blue_end = 6532.8

        red_start = 6592.8
        red_end = 6632.8
        """
        blue_start = 6465
        blue_end = 6470

        red_start = 6620
        red_end = 6635

    elif line_name.lower() in {'hbeta', '[oiii]5007'}:

        blue_start = line_loc - 35
        blue_end = line_loc - 30

        red_start = line_loc + 30
        red_end = line_loc + 35

    elif line_name.lower() in {'[fexiv]5304', 'heii4686'}:

        blue_start = line_loc - 30
        blue_end = line_loc - 25

        red_start = line_loc + 25
        red_end = line_loc + 30

    elif line_name.lower() in {'[oiii]4959'}:

        blue_start = line_loc - 20
        blue_end = line_loc - 15

        red_start = line_loc + 15
        red_end = line_loc + 20

    else:
        blue_start = line_loc - 50
        blue_end = line_loc - 30

        red_start = line_loc + 30
        red_end = line_loc + 50

    continuum_regions = [blue_start, blue_end, red_start, red_end]

    blue_middle = (blue_start + blue_end) / 2
    red_middle = (red_start + red_end) / 2

    blue_flux = []
    red_flux = []

    for xx, wavelength_point in enumerate(spec_region):
        if blue_start < spec_region[xx] < blue_end:
            blue_flux.append(flux[xx])
        if red_start < spec_region[xx] < red_end:
            red_flux.append(flux[xx])

    try:
        m = (np.nanmean(red_flux) - np.nanmean(blue_flux)) / (red_middle - blue_middle)
        c = np.nanmean(blue_flux) - (m * blue_middle)

    except FloatingPointError:
        print("Something has gone wrong with the selection of part of the continuum - likely no valid points"
              "\nSkipping object for now")
        print(f"{object_name}\t{line_name}")

        print(spec_region)

        return [], []

    continuum = []

    for xx, velocity in enumerate(spec_region):
        continuum.append((m * spec_region[xx]) + c)

    try:
        scaled_flux = np.array(flux) / np.array(continuum)
    except:
        print("Continuum Scaling Failure")
        print(flux)
        print(continuum)

        scaled_flux = np.array(flux)

    return continuum, scaled_flux, continuum_regions


def eqw_finder(flux, continuum, xaxis, object_name, start="A", stop="B", xstep="C"):
    """Equivalent Width Finder - Configured for specifically for Hirogen - An older version that is not currently in use

    Uses the 'area below the continuum' integration method

    Will work for both absorption and emission features.
    Always returns the CALCULATED value. Emission = +ve, Absorption = -ve

    Requires: Flux, Continuum, Wavelength/Velocity(xaxis) Start, Stop (optional, if none given will run for the
    full number of points),
    xstep (optional, works this out itself if the axis is in wavelength - MUST be given if a velocity)
    Returns: [EQW, EQW Points, EQW_Region_Flux]"""

    if len(flux) != len(continuum) or len(flux) != len(xaxis):
        print("It appears one of the data lists passed to the EQW Finder isn't the correct length."
              "\nPlease check and try again\n\nLengths\nFlux:{}\nContinuum:{}\nX axis:{}"
              .format(str(len(flux)), str(len(continuum)), str(len(xaxis))))
        sys.exit()

    if start == "A":
        start = xaxis[0]

    if stop == "B":
        stop = xaxis[-1]

    if xstep == "C":
        xstep = xaxis[1] - xaxis[0]

    eqw_points = []
    flux_points = []
    velocity_points = []

    for xx, item in enumerate(xaxis):

        if start < xaxis[xx] < stop:
            flux_point = flux[xx] / continuum[xx]
            eqw_point = flux_point * xstep
            eqw_points.append(eqw_point)
            flux_points.append(flux_point)
            velocity_points.append(xaxis[xx])

    # eqw = np.nansum(eqw_points) * -1  # -1 makes sure absorption features are -ve and emission features are +ve
    eqw = np.nansum(eqw_points)

    if Debug:

        if eqw == np.nan:
            print("Nan EQW Detected")
            print(object_name)
            print(eqw)
            print(eqw_points)

        print(flux_points)
        print(eqw_points)
        print(continuum)
        print(eqw)

    return eqw, eqw_points, flux_points, velocity_points


def eqw_measurement(
        flux,
        true_continuum,
        scaled_flux,
        xaxis,
        line_loc,
        measurement_length=12
):
    """Newer version of the EQW measurement - still uses the under the continuum summation
    Defaults to measuring 12 angstroms either side of the given line location though this can be changed as needed"""

    prepared_subtracted_flux = []
    prepared_scaled_flux = []
    prepared_xaxis = []

    for ii, item in enumerate(xaxis):
        if (line_loc - measurement_length) < xaxis[ii] < (line_loc + measurement_length):
            prepared_subtracted_flux.append(flux[ii] - true_continuum[ii])
            prepared_scaled_flux.append(scaled_flux[ii])
            prepared_xaxis.append(xaxis[ii])

    continuum = np.ones(len(prepared_scaled_flux))  # Assumes the input flux has been scaled to 1
    sum = np.sum((continuum - prepared_scaled_flux) / continuum)  # Now uses the standard definition of EQW
    line_flux = np.sum(prepared_subtracted_flux)

    delta = xaxis[1] - xaxis[0]
    eqw = delta * sum

    return eqw, line_flux


def lin_con_candidate_line_identifier(line_name, line_peqw, peqw_threshold, continuum_removed_flux, shift_points,
                                      peak_threshold, above_continuum_threshold, peak_max_region, peak_region,
                                      feature_region, peak_region_minima_threshold):
    """Line scoring identification
    Checks the given line against the scoring criteria and determines if the line passes the cut"""

    feature_region_shift = []
    feature_region_flux = []

    peak_region_shift = []
    peak_region_flux = []

    central_region_shift = []
    central_region_flux = []

    for ii, item in enumerate(shift_points):

        if -1 * feature_region < shift_points[ii] < feature_region:
            feature_region_shift.append(shift_points[ii])
            feature_region_flux.append(continuum_removed_flux[ii])

            if -1 * peak_region < shift_points[ii] < peak_region:
                peak_region_shift.append(shift_points[ii])
                peak_region_flux.append(continuum_removed_flux[ii])

                if -1 * peak_max_region < shift_points[ii] < peak_max_region:
                    central_region_shift.append(shift_points[ii])
                    central_region_flux.append(continuum_removed_flux[ii])

    try:
        peak_region_max = np.nanmax(peak_region_flux)
        peak_region_max_loc = peak_region_shift[np.nanargmax(peak_region_flux)]
        peak_region_mean = np.nanmean(peak_region_flux)
    except ValueError:
        print(f"Peak region max cannot be located: likely spectral gap.\n{line_name}")
        peak_region_max = np.nan
        peak_region_max_loc = np.nan
        peak_region_mean = np.nan

    try:
        feature_region_max = np.nanmax(feature_region_flux)
        feature_region_max_loc = feature_region_shift[np.nanargmax(feature_region_flux)]
        feature_region_mean = np.nanmean(feature_region_flux)
    except ValueError:
        print(f"Feature region max cannot be located: likely spectral gap.\n{line_name}")
        feature_region_max = np.nan
        feature_region_max_loc = np.nan
        feature_region_mean = np.nan

    if Debug:
        print(line_name)
        print(line_peqw)
        print(peak_region_max)
        print(peak_region_max_loc)
        print(peak_region_mean)
        print(feature_region_max)
        print(feature_region_max_loc)
        print(feature_region_mean)
        print(feature_region_mean * above_continuum_threshold)
        print("\n")

    if line_peqw <= peqw_threshold and peak_region_max > peak_threshold and \
            peak_region_max > feature_region_mean * above_continuum_threshold and \
            (-1 * peak_max_region) < peak_region_max_loc < peak_max_region and \
            np.isnan(np.sum(peak_region_flux)) == False and \
            np.nanmin(central_region_flux) >= peak_region_minima_threshold:

        if line_name == '[FeXIV]5304' or line_name == '[FeVII]6088' or line_name == '[FeX]6376' or \
                line_name == '[FeXI]7894':
            return 2, peak_region_mean, peak_region_max, peak_region_max_loc, feature_region_mean, feature_region_max, \
                   feature_region_max_loc, 'Pass', 'Pass: +2 Score'
        else:
            return 1, peak_region_mean, peak_region_max, peak_region_max_loc, feature_region_mean, feature_region_max, \
                   feature_region_max_loc, 'Pass', 'Pass: +1 Score'

    else:
        return 0, peak_region_mean, peak_region_max, peak_region_max_loc, feature_region_mean, feature_region_max, \
               feature_region_max_loc, 'Fail', 'Fail'


def lin_con_candidate_line_identifier_verbose(
        line_name, line_peqw, peqw_threshold, continuum_removed_flux, shift_points,
        peak_threshold, above_continuum_threshold, peak_max_region, peak_region, feature_region,
        peak_region_minima_threshold):
    """Line scoring identification but with more informative information on which criteria a line failed the tests"""

    feature_region_shift = []
    feature_region_flux = []

    peak_region_shift = []
    peak_region_flux = []

    central_region_shift = []
    central_region_flux = []

    for ii, item in enumerate(shift_points):

        if -1 * feature_region < shift_points[ii] < feature_region:
            feature_region_shift.append(shift_points[ii])
            feature_region_flux.append(continuum_removed_flux[ii])

            if -1 * peak_region < shift_points[ii] < peak_region:
                peak_region_shift.append(shift_points[ii])
                peak_region_flux.append(continuum_removed_flux[ii])

                if -1 * peak_max_region < shift_points[ii] < peak_max_region:
                    central_region_shift.append(shift_points[ii])
                    central_region_flux.append(continuum_removed_flux[ii])
    try:
        peak_region_max = np.nanmax(peak_region_flux)
        peak_region_max_loc = peak_region_shift[np.nanargmax(peak_region_flux)]
        peak_region_mean = np.nanmean(peak_region_flux)
    except ValueError:
        print(f"Peak region max cannot be located: likely spectral gap.\n{line_name}")
        peak_region_max = np.nan
        peak_region_max_loc = np.nan
        peak_region_mean = np.nan

    try:
        feature_region_max = np.nanmax(feature_region_flux)
        feature_region_max_loc = feature_region_shift[np.nanargmax(feature_region_flux)]
        feature_region_mean = np.nanmean(feature_region_flux)
    except ValueError:
        print(f"Feature region max cannot be located: likely spectral gap.\n{line_name}")
        feature_region_max = np.nan
        feature_region_max_loc = np.nan
        feature_region_mean = np.nan

    if np.isnan(np.sum(peak_region_flux)) == False and len(peak_region_flux):
        print("No nans in peak detection region")

        if np.nanmin(central_region_flux) >= peak_region_minima_threshold:
            print(np.nanmin(central_region_flux), peak_region_minima_threshold)
            print("Central region minimum above threshold")

            if line_peqw <= peqw_threshold:
                print(line_peqw, peqw_threshold)
                print("Line pEQW exceeds threshold")

                if peak_region_max > peak_threshold:
                    print(peak_region_max, peak_threshold)
                    print("Region max exceeds threshold")

                    if peak_region_max > feature_region_mean * above_continuum_threshold:
                        print(peak_region_max, (peak_region_mean * above_continuum_threshold))

                        print("Region max exceeds mean * continuum threshold")

                        if (-1 * peak_max_region) < peak_region_max_loc < peak_max_region:
                            print(shift_points[np.nanargmax(continuum_removed_flux)])
                            print("Max occurs within allowed velocity region")

                            if line_name == '[FeXIV]5304' or line_name == '[FeVII]6088' or line_name == '[FeX]6376' or \
                                    line_name == '[FeXI]7894':
                                print(f"{line_name}: Pass")
                                return 2, peak_region_mean, peak_region_max, peak_region_max_loc, feature_region_mean, \
                                       feature_region_max, feature_region_max_loc, 'Pass', 'Pass: +2 Score'
                            else:
                                print(f"{line_name}: Pass")
                                return 1, peak_region_mean, peak_region_max, peak_region_max_loc, feature_region_mean, \
                                       feature_region_max, feature_region_max_loc, 'Pass', 'Pass: +1 Score'
                        else:
                            print(f"{line_name}:Region Max Occurs Outside Allowed Velocity Region Fail")
                            return 0, peak_region_mean, peak_region_max, peak_region_max_loc, feature_region_mean, \
                                   feature_region_max, feature_region_max_loc, 'Fail: Peak Location', \
                                   'Fail: Peak Outside Allowed Region'
                    else:
                        print(f"{line_name}: Region Max Above Mean Threshold Fail")
                        return 0, peak_region_mean, peak_region_max, peak_region_max_loc, feature_region_mean, \
                               feature_region_max, feature_region_max_loc, 'Fail: Peak strength above mean', \
                               'Fail: Peak Not Sufficiently Above Region Mean'
                else:
                    print(f"{line_name}: Peak Threshold Fail")
                    return 0, peak_region_mean, peak_region_max, peak_region_max_loc, feature_region_mean, \
                           feature_region_max, feature_region_max_loc, 'Fail: Absolute peak strength', \
                           'Fail: Peak Does Not Reach Min Height Threshold'

            else:
                print(f"{line_name}: pEQW cut Fail")
                return 0, peak_region_mean, peak_region_max, peak_region_max_loc, feature_region_mean, \
                       feature_region_max, feature_region_max_loc, \
                       'Fail: Feature pEQW', 'Fail: Feature Does Not Meet pEQW Cut'

        else:
            print(f"{line_name}: Peak region minimum falls below threshold")
            return 0, peak_region_mean, peak_region_max, peak_region_max_loc, feature_region_mean, feature_region_max, \
                   feature_region_max_loc, 'Fail: Central region minimum below threshold', \
                   'Fail: Central Region Minimum Below Threshold'

    else:
        print(f"{line_name}: nan detected in peak detection region")
        return 0, peak_region_mean, peak_region_max, peak_region_max_loc, feature_region_mean, feature_region_max, \
               feature_region_max_loc, 'Fail: Peak region nan detected', \
               'Fail: nan Detected In Peak Identification Region or peak region not found'


"""Spectral readin functions"""


########################################################################
# Observed Wavelength to Rest Wavelength Converter
########################################################################

def rest_wavelength_converter(observer_frame_wave, z):
    i = 0

    RestWave = []

    for item in observer_frame_wave:
        if observer_frame_wave[i] != np.nan:
            RestSingle = (observer_frame_wave[i] / (1 + z))
            RestWave.append(RestSingle)
        else:
            RestWave.append(np.nan)
        i = i + 1

    RestWave = np.array(RestWave)

    return RestWave


def spec_reader(filepath, verbose=False):
    """Generic fits spectrum reader - does not include a flux error"""
    hdul = fits.open(filepath)
    header = hdul[0].header

    if verbose:
        hdul.info()
        print(hdul)
        print(header)

    data = hdul[1].data
    hdul.close()

    return header, data


def spec_reader_ascii(filepath, verbose=False):
    """Generic ascii spectrum reader"""

    asciiSpectrum = np.genfromtxt(filepath, unpack=True, comments='#')
    ObservedWave = np.array(asciiSpectrum[0])
    Flux = np.array(asciiSpectrum[1])

    # This will generate a fake flux error column if one isn't present - not particularly desirable
    try:
        FluxError = np.array(asciiSpectrum[2])
    except IndexError:
        FluxError = np.array([np.float64(0)] * len(Flux))

    return ObservedWave, Flux, FluxError


def spec_reader_fits(filepath, verbose=False):
    """Generic fits spectrum reader"""
    """Assumes that the data is all held within tables rather than the header"""

    hdul = fits.open(filepath, comments='#')

    header = hdul[0].header
    spectrum_data = hdul[1].data

    if verbose:
        print(header)

    """Assumes the spectrum data is in the format:
    Wavelength, Flux, FluxError"""

    ObservedWave = np.array(spectrum_data[0][0])
    Flux = np.array(spectrum_data[0][1])

    # This will generate a fake flux error column if one isn't present - not particularly desirable
    try:
        FluxError = np.array(spectrum_data[0][2])
    except IndexError:
        FluxError = np.array([np.float64(0)] * len(Flux))

    return ObservedWave, Flux, FluxError


def sdss_spectra_file_path_generator(
        general_path,
        plate_list,
        mjd_list,
        fiber_id_list,
        survey_list,
        run2d_list,
        override_path_flag_list=[],
        override_path_list=[]
):
    """Generates the filepaths for the SDSS spectra assuming they follow the standard file structure from the direct
    download"""

    file_paths = []

    for ii, item in enumerate(plate_list):

        if override_path_flag_list[ii] == 0:

            if survey_list[ii] in ['sdss', 'segue1', 'segue2']:

                file_paths.append(
                    f"{general_path}/dr16/sdss/spectro/redux/{run2d_list[ii]}/spectra/{plate_list[ii]}/"
                    f"spec-{plate_list[ii]}-{mjd_list[ii]}-{fiber_id_list[ii]}.fits"
                )

            elif survey_list[ii] in ['eboss', 'boss']:

                file_paths.append(
                    f"{general_path}/dr16/eboss/spectro/redux/{run2d_list[ii]}/spectra/full/{plate_list[ii]}/"
                    f"spec-{plate_list[ii]}-{mjd_list[ii]}-{fiber_id_list[ii]}.fits"
                )

            else:
                print("Not sure how to handle this survey - exiting for now")
                print(f"{ii}\t{survey_list[ii]}")
                sys.exit()

        else:
            file_paths.append(
                f"{override_path_list[ii]}"
            )

    return file_paths


def sdss_peaks_spectra_file_path_generator(general_path, plate_list, mjd_list, fiber_id_list, survey_list, run2d_list,
                                           override_path_flag_list=[], override_path_list=[], Mode='peaks'):
    file_paths = []

    if Mode.lower() == 'peaks':

        for ii, item in enumerate(plate_list):

            if override_path_flag_list[ii] == 0:

                if survey_list[ii] in ['sdss', 'segue1', 'segue2']:

                    file_paths.append(
                        f"{general_path}/dr16/sdss/spectro/redux/{run2d_list[ii]}/spectra/{plate_list[ii]}/"
                        f"spec-{plate_list[ii]}-{mjd_list[ii]}-{fiber_id_list[ii]}-peaks.fits"
                    )

                elif survey_list[ii] in ['eboss', 'boss']:

                    file_paths.append(
                        f"{general_path}/dr16/eboss/spectro/redux/{run2d_list[ii]}/spectra/full/{plate_list[ii]}/"
                        f"spec-{plate_list[ii]}-{mjd_list[ii]}-{fiber_id_list[ii]}-peaks.fits"
                    )

                else:
                    print("Not sure how to handle this survey - exiting for now")
                    print(f"{ii}\t{survey_list[ii]}")
                    sys.exit()

            else:
                file_paths.append(
                    f"{override_path_list[ii]}"
                )

    elif Mode.lower() == 'fakes':

        for ii, item in enumerate(plate_list):

            if override_path_flag_list[ii] == 0:

                if survey_list[ii] in ['sdss', 'segue1', 'segue2']:

                    file_paths.append(
                        f"{general_path}/dr16/sdss/spectro/redux/{run2d_list[ii]}/spectra/{plate_list[ii]}/"
                        f"spec-{plate_list[ii]}-{mjd_list[ii]}-{fiber_id_list[ii]}-fake.fits"
                    )

                elif survey_list[ii] in ['eboss', 'boss']:

                    file_paths.append(
                        f"{general_path}/dr16/eboss/spectro/redux/{run2d_list[ii]}/spectra/full/{plate_list[ii]}/"
                        f"spec-{plate_list[ii]}-{mjd_list[ii]}-{fiber_id_list[ii]}-fake.fits"
                    )

                else:
                    print("Not sure how to handle this survey - exiting for now")
                    print(f"{ii}\t{survey_list[ii]}")
                    sys.exit()

            else:
                file_paths.append(
                    f"{override_path_list[ii]}"
                )

    else:
        print('mode selection not recognised.\nPlease check and try again.')
        sys.exit()

    return file_paths


def ascii_spectrum_reader(
        filepath,
        z,
        extinction,
        smoothing=True, smoothing_boxcar=5,
        median_filter=False, med_filter_kernel=3,
        sav_gol=False, sav_gol_window=5, sav_gol_poly=2,
        flag_check=False,
        z_correction_flag=0, extinction_correction_flag=0
):
    """ Very similar in function to SDSS spectrum reader so comments are only included for the differences
    Makes a few assumptions about the file and sets the unknown SDSS values to defaults """

    spectrum = spec_reader_ascii(filepath=filepath)
    wave = spectrum[0]
    flux = spectrum[1] * 1e17
    flux_err = spectrum[2] * 1e17

    # The multiplication factors are to convert into the default SDSS data format where these are baked into the data
    # This allows for easier direct comparisons between original and follow up spectra - these may need to be removed
    # if the spectrum already has these factors included e.g. DESI spectra

    lamb_observed = wave * u.AA
    spec_res = (lamb_observed[1] - lamb_observed[0])  # Assumes that the file has a fixed wavelength resolution

    flux_in = flux * u.Unit('erg cm-2 s-1 AA-1')
    error = flux_err

    all_flag_percentage = -999
    counted_flag_percentage = -999

    flux = flux_in

    ###########
    # Smoothing
    ###########

    if smoothing:  # If active runs a boxcar smooth, defaults to a width of 5 but can be set via the function call
        smoothed_flux = convolve(flux, Box1DKernel(smoothing_boxcar))
        flux = smoothed_flux

    ###########
    # Median Filter
    ###########

    if median_filter:  # If active runs a median filter smooth, defaults to a kernel size of 10 but can be set via the
        # function call
        filtered_flux = medfilt(flux, kernel_size=med_filter_kernel)
        flux = filtered_flux

    ###########
    # Savitzky-Golay Filter
    ###########

    if sav_gol:
        # If active runs a Savitzky-Golay filter on the spectral data, defaults to a window size of 5 and 2nd order
        # polynomial though can be configured via the call. All other parameters are the scipy defaults
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html

        filtered_flux = savgol_filter(x=flux, window_length=sav_gol_window, polyorder=sav_gol_poly)
        flux = filtered_flux

    ###########
    # Extinction Correction
    ###########

    # SDSS Extinction comes from the mean value of the per band extinctions converted using the constants given in the
    # 2002 Early Data release paper - this uses the Schlegel, Finkbeiner and Davis dust maps

    # New version of the extinction correction using the "dust_extinction" package

    if extinction_correction_flag == 0:
        extinction_model = F99(Rv=3.1)
        flux_extinction_corrected = flux / extinction_model.extinguish(
            lamb_observed, Ebv=extinction
        )
        flux = flux_extinction_corrected

    ###########
    # Redshift Correction
    ###########

    if z_correction_flag == 0:
        lamb_rest = rest_wavelength_converter(observer_frame_wave=lamb_observed.value, z=z) * u.AA
    else:
        lamb_rest = lamb_observed

    # These flags lists are empty as an ascii file won't have them

    flags = []
    masked_flags = []

    return lamb_rest, flux, error, spec_res, flags, masked_flags, all_flag_percentage, counted_flag_percentage


def fits_spectrum_reader(
        filepath,
        z,
        extinction,
        smoothing=True, smoothing_boxcar=5,
        median_filter=False, med_filter_kernel=3,
        sav_gol=False, sav_gol_window=5, sav_gol_poly=2,
        flag_check=False,
        z_correction_flag=0, extinction_correction_flag=0
):
    """ Very similar in function to SDSS spectrum reader so comments are only included for the differences
    Makes a few assumptions about the file and sets the unknown SDSS values to defaults """

    spectrum = spec_reader_fits(filepath=filepath)
    wave = spectrum[0]
    flux = spectrum[1]  # * 1e17 # This factor appears to be included in the MMT spectra at least
    flux_err = spectrum[2]  # * 1e17

    # The multiplication factors are to convert into the default SDSS data format where these are baked into the data
    # This allows for easier direct comparisons between original and follow up spectra

    lamb_observed = wave * u.AA
    spec_res = (lamb_observed[1] - lamb_observed[0])  # Assumes that the file has a fixed wavelength resolution

    flux_in = flux * u.Unit('erg cm-2 s-1 AA-1')
    error = flux_err

    all_flag_percentage = -999
    counted_flag_percentage = -999

    flux = flux_in

    ###########
    # Smoothing
    ###########

    if smoothing:  # If active runs a boxcar smooth, defaults to a width of 5 but can be set via the function call
        smoothed_flux = convolve(flux, Box1DKernel(smoothing_boxcar))
        flux = smoothed_flux

    ###########
    # Median Filter
    ###########

    if median_filter:  # If active runs a median filter smooth, defaults to a kernel size of 10 but can be set via the
        # function call
        filtered_flux = medfilt(flux, kernel_size=med_filter_kernel)
        flux = filtered_flux

    ###########
    # Savitzky-Golay Filter
    ###########

    if sav_gol:
        # If active runs a Savitzky-Golay filter on the spectral data, defaults to a window size of 5 and 2nd order
        # polynomial though can be configured via the call. All other parameters are the scipy defaults
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html

        filtered_flux = savgol_filter(x=flux, window_length=sav_gol_window, polyorder=sav_gol_poly)
        flux = filtered_flux

    ###########
    # Extinction Correction
    ###########

    # SDSS Extinction comes from the mean value of the per band extinctions converted using the constants given in the
    # 2002 Early Data release paper - this uses the Schlegel, Finkbeiner and Davis dust maps

    # New version of the extinction correction using the "dust_extinction" package
    if extinction_correction_flag == 0:
        extinction_model = F99(Rv=3.1)
        flux_extinction_corrected = flux / extinction_model.extinguish(
            lamb_observed, Ebv=extinction
        )
        flux = flux_extinction_corrected

    ###########
    # Redshift Correction
    ###########

    if z_correction_flag == 0:
        lamb_rest = rest_wavelength_converter(observer_frame_wave=lamb_observed.value, z=z) * u.AA
    else:
        lamb_rest = lamb_observed

    # These flags lists are empty as an ascii file won't have them

    flags = []
    masked_flags = []

    return lamb_rest, flux, error, spec_res, flags, masked_flags, all_flag_percentage, counted_flag_percentage


def sdss_spectrum_reader(
        filepath,
        z,
        extinction,
        smoothing=True, smoothing_boxcar=5,
        median_filter=False, med_filter_kernel=3,
        header_print=False, flag_check=False,
        sav_gol=False, sav_gol_window=5, sav_gol_poly=2,
        z_correction_flag=0, extinction_correction_flag=0
):
    """Reads in and does initial processing on SDSS spectra"""

    spectrum = spec_reader(filepath=filepath)
    spec_header = spectrum[0]

    if header_print:
        print(spec_header)

    data = spectrum[1]

    lamb_observed = 10 ** data['loglam'] * u.AA  # Observed wavelength
    spec_res = lamb_observed[1] - lamb_observed[0]  # Spectral resolution
    # This is VERY rough currently given the non fixed SDSS spectrum resolution in angstroms

    flux_in = data['flux'] * u.Unit('erg cm-2 s-1 AA-1')  # Flux
    error = data['ivar']  # Flux error

    flags = data['and_mask']  # Bitmask and mask map

    ignored_flag_values = np.array([0, 4, 16, 131076])  # Pixels with these flags are not removed

    bright_sky = np.argwhere(flags == 8388608)  # The bitmask for bright sky lines specifically
    all_flags = np.argwhere(flags > 0)

    # masked_flags = np.argwhere(np.isin(flags, ignored_flag_values, invert=True))
    masked_flags = np.argwhere(np.isin(flags, bright_sky, invert=False))
    # Swapping this out so only bright sky lines are masked for now

    if Debug:
        print(f"The bright sky flagged pixels are: {bright_sky}")
        print(f"The pixels masked for any reason are: {masked_flags}")

    all_flag_percentage = (len(all_flags) / len(lamb_observed)) * 100
    counted_flag_percentage = (len(masked_flags) / len(lamb_observed)) * 100
    # Calculates the percentage of pixels flagged relative to the full spectrum

    for ii, item in enumerate(flux_in):  # This loops sets the flagged pixels to be nan
        if ii in masked_flags:
            flux_in[ii] = np.nan

    flux = flux_in

    ###########
    # Smoothing
    ###########

    if smoothing:  # If active runs a boxcar smooth, defaults to a width of 5 but can be set via the function call
        smoothed_flux = convolve(flux, Box1DKernel(smoothing_boxcar))
        flux = smoothed_flux

    ###########
    # Median Filter
    ###########

    if median_filter:  # If active runs a median filter smooth, defaults to a kernel size of 3 but can be set via the
        # function call
        filtered_flux = medfilt(flux, kernel_size=med_filter_kernel)
        flux = filtered_flux

    ###########
    # Savitzky-Golay Filter
    ###########

    if sav_gol:
        # If active runs a Savitzky-Golay filter on the spectral data, defaults to a window size of 5 and 2nd order
        # polynomial though can be configured via the call. All other parameters are the scipy defaults
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html

        filtered_flux = savgol_filter(x=flux, window_length=sav_gol_window, polyorder=sav_gol_poly)
        flux = filtered_flux

    ###########
    # Extinction Correction
    ###########

    # SDSS Extinction comes from the mean value of the per band extinctions converted using the constants given in the
    # 2002 Early Data release paper - this uses the Schlegel, Finkbeiner and Davis dust maps

    # New version of the extinction correction using the "dust_extinction" package
    if extinction_correction_flag == 0:
        extinction_model = F99(Rv=3.1)
        flux_extinction_corrected = flux / extinction_model.extinguish(
            lamb_observed, Ebv=extinction
        )
        flux = flux_extinction_corrected

    ###########
    # Redshift Correction
    ###########

    # Applies a spectral redshift correction
    if z_correction_flag == 0:
        lamb_rest = rest_wavelength_converter(observer_frame_wave=lamb_observed.value, z=z) * u.AA
    else:
        lamb_rest = lamb_observed

    if flag_check:  # If active displays the flagged pixels and the associated wavelengths

        for ii, item in enumerate(flux):
            print(format(flags[ii], '031b'))
            print(lamb_rest[ii])

    return lamb_rest, flux, error, spec_res, flags, masked_flags, all_flag_percentage, counted_flag_percentage


"""DESI specific functions"""


def DESI_check_env():
    print("Now running a quick check to confirm that the environment variables have been configured correctly.\n")

    for env in ('DESI_SPECTRO_REDUX', 'SPECPROD'):
        if env in os.environ:
            print('${}={}'.format(env, os.getenv(env)))
        else:
            print('Required environment variable {} not set!'.format(env))

    reduxdir = desispec.io.specprod_root()
    if not os.path.exists(reduxdir):
        print("ERROR: {} doesn't exist; check $DESI_SPECTRO_REDUX/$SPECPROD".format(reduxdir))
    else:
        print('OK: {} exists\n'.format(reduxdir))


def desi_datafile_reader(file):
    """Processes the Hirogen analysis output file making its contents available to the current file"""

    data = np.genfromtxt(file, unpack=True, comments='#', dtype=str, delimiter="\t", skip_header=1)

    # ID
    col1 = data[0].astype(int)

    # Date
    col2 = data[1].astype(int)

    # Tile
    col3 = data[2].astype(int)

    # Petal
    col4 = data[3].astype(int)

    # Fiber
    col5 = data[4].astype(int)

    # FiberStatus
    col6 = data[5].astype(int)

    # Perfile Index
    col7 = data[6].astype(int)

    # RA + DEC
    col8 = (np.array(data[7])).astype(float)
    col9 = (np.array(data[8])).astype(float)

    # Redshift + Err
    col10 = (np.array(data[9])).astype(float)
    col11 = (np.array(data[10])).astype(float)

    # EBV
    col12 = (np.array(data[11])).astype(float)

    # line pEQWs

    # Balmer
    col13 = (np.array(data[12])).astype(float)
    col14 = (np.array(data[13])).astype(float)
    col15 = (np.array(data[14])).astype(float)
    col16 = (np.array(data[15])).astype(float)

    # NII
    col17 = (np.array(data[16])).astype(float)
    col18 = (np.array(data[17])).astype(float)

    # Fe
    col19 = (np.array(data[18])).astype(float)
    col20 = (np.array(data[19])).astype(float)
    col21 = (np.array(data[20])).astype(float)
    col22 = (np.array(data[21])).astype(float)
    col23 = (np.array(data[22])).astype(float)
    col24 = (np.array(data[23])).astype(float)
    col25 = (np.array(data[24])).astype(float)

    # O
    col26 = (np.array(data[25])).astype(float)
    col27 = (np.array(data[26])).astype(float)
    col28 = (np.array(data[27])).astype(float)
    col29 = (np.array(data[28])).astype(float)

    # He
    col30 = (np.array(data[29])).astype(float)
    col31 = (np.array(data[30])).astype(float)

    # NaID
    col32 = (np.array(data[31])).astype(float)

    # S
    col33 = (np.array(data[32])).astype(float)
    col34 = (np.array(data[33])).astype(float)
    col35 = (np.array(data[34])).astype(float)

    # line Fluxes

    # Balmer
    col36 = (np.array(data[35])).astype(float)
    col37 = (np.array(data[36])).astype(float)
    col38 = (np.array(data[37])).astype(float)
    col39 = (np.array(data[38])).astype(float)

    # NII
    col40 = (np.array(data[39])).astype(float)
    col41 = (np.array(data[40])).astype(float)

    # Fe
    col42 = (np.array(data[41])).astype(float)
    col43 = (np.array(data[42])).astype(float)
    col44 = (np.array(data[43])).astype(float)
    col45 = (np.array(data[44])).astype(float)
    col46 = (np.array(data[45])).astype(float)
    col47 = (np.array(data[46])).astype(float)
    col48 = (np.array(data[47])).astype(float)

    # O
    col49 = (np.array(data[48])).astype(float)
    col50 = (np.array(data[49])).astype(float)
    col51 = (np.array(data[50])).astype(float)
    col52 = (np.array(data[51])).astype(float)

    # He
    col53 = (np.array(data[52])).astype(float)
    col54 = (np.array(data[53])).astype(float)

    # NaID
    col55 = (np.array(data[54])).astype(float)

    # S
    col56 = (np.array(data[55])).astype(float)
    col57 = (np.array(data[56])).astype(float)
    col58 = (np.array(data[57])).astype(float)

    # Ratios
    col59 = (np.array(data[58])).astype(float)
    col60 = (np.array(data[59])).astype(float)
    col61 = (np.array(data[60])).astype(float)
    col62 = (np.array(data[61])).astype(float)
    col63 = (np.array(data[62])).astype(float)
    col64 = (np.array(data[63])).astype(float)

    # Lick Indices
    col65 = (np.array(data[64])).astype(float)
    col66 = (np.array(data[65])).astype(float)
    col67 = (np.array(data[66])).astype(float)
    col68 = (np.array(data[67])).astype(float)
    col69 = (np.array(data[68])).astype(float)
    col70 = (np.array(data[69])).astype(float)

    # Flags
    col71 = (np.array(data[70])).astype(int)
    col72 = (np.array(data[71])).astype(int)
    col73 = (np.array(data[72])).astype(int)
    col74 = (np.array(data[73])).astype(int)
    col75 = (np.array(data[74])).astype(int)
    col76 = (np.array(data[75])).astype(int)
    col77 = (np.array(data[76])).astype(int)

    # Pixel / Spectrum Quality Information
    col78 = (np.array(data[77])).astype(float)
    col79 = (np.array(data[78])).astype(float)

    # ECLE Score
    col80 = (np.array(data[79])).astype(int)

    # Flux Data - Uses FiberFluxes but FLUX_IVAR
    # G
    col81 = (np.array(data[80])).astype(float)
    col82 = (np.array(data[81])).astype(float)
    # R
    col83 = (np.array(data[82])).astype(float)
    col84 = (np.array(data[83])).astype(float)
    # Z
    col85 = (np.array(data[84])).astype(float)
    col86 = (np.array(data[85])).astype(float)

    # WISE
    # W1
    col87 = (np.array(data[86])).astype(float)
    col88 = (np.array(data[87]))
    # W2
    col89 = (np.array(data[88])).astype(float)
    col90 = (np.array(data[89]))

    for ii, item in enumerate(col88):
        if col88[ii] == '--':
            col88[ii] = -995
        if col90[ii] == '--':
            col90[ii] = -995

    col88 = col88.astype(float)
    col90 = col90.astype(float)

    # Spectrum ID Key
    col91 = data[90]

    # Settings Config Store
    col92 = data[91]

    return \
        col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, \
        col17, col18, col19, col20, col21, col22, col23, col24, col25, col26, col27, col28, col29, col30, col31, \
        col32, col33, col34, col35, col36, col37, col38, col39, col40, col41, col42, col43, col44, col45, col46, \
        col47, col48, col49, col50, col51, col52, col53, col54, col55, col56, col57, col58, col59, col60, col61, \
        col62, col63, col64, col65, col66, col67, col68, col69, col70, col71, col72, col73, col74, col75, col76, \
        col77, col78, col79, col80, col81, col82, col83, col84, col85, col86, col87, col88, col89, col90, col91, col92


def desi_datafile_pre_Hirogen_analysis_reader(file):
    """Processes the pre-Hirogen analysis metafile making its contents available to the current file"""

    data = np.genfromtxt(file, unpack=True, comments='#', dtype=str, delimiter="\t", skip_header=1)

    # ID
    col1 = data[0].astype(int)

    # Date
    col2 = data[1].astype(int)

    # Tile
    col3 = data[2].astype(int)

    # PETAL
    col4 = data[3].astype(int)

    # FIBER
    col5 = data[4].astype(int)

    # FIBER_STATUS
    col6 = data[5].astype(int)

    # SPECTRUM_NUMBER
    col7 = data[6].astype(int)

    # RA
    col8 = data[7].astype(float)

    # DEC
    col9 = data[8].astype(float)

    # REDSHIFT
    col10 = data[9].astype(float)

    # REDSHIFT_ERR
    col11 = data[10].astype(float)

    # EBV
    col12 = data[11].astype(float)

    # GFiberFlux
    col13 = data[12].astype(float)

    # GFiberIvar
    col14 = data[13].astype(float)

    # RFiberFlux
    col15 = data[14].astype(float)

    # RFiberIvar
    col16 = data[15].astype(float)

    # ZFiberFlux
    col17 = data[16].astype(float)

    # ZFiberIvar
    col18 = data[17].astype(float)

    # W1Flux
    col19 = data[18].astype(float)

    # W1Ivar
    col20 = data[19].astype(float)

    # W2Flux
    col21 = data[20].astype(float)

    # W2Ivar
    col22 = data[21].astype(float)

    return col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, \
           col17, col18, col19, col20, col21, col22


"""Tools"""


# Largely isolated components of larger functions

def spectral_flux_smooth(flux, smoothing_boxcar=5):
    """Runs a boxcar smooth, defaults to a width of 5 but can be set via the function call"""

    smoothed_flux = convolve(flux, Box1DKernel(smoothing_boxcar))

    return smoothed_flux


def spectral_median_filter(flux, med_filter_kernel=3):
    """Runs a median filter, defaults to a kernel of size 3 but can be set via the function call"""

    filtered_flux = medfilt(flux, kernel_size=med_filter_kernel)

    return filtered_flux


def spectral_extinction_correction(flux, wave, extinction):
    extinction_model = F99(Rv=3.1)
    flux_extinction_corrected = flux / extinction_model.extinguish(wave, Ebv=extinction)

    return flux_extinction_corrected


# Fast Fourier Transform Smoothing Function

def fft_filter_signal(signal, threshold=1e8):
    fourier = rfft(signal)
    frequencies = rfftfreq(signal.size, d=20e-3 / signal.size)
    fourier[frequencies > threshold] = 0
    # without setting n = len(signal) one flux point is removed
    return irfft(fourier, n=len(signal))


# Gaussian Smoothing Function

def gaussian_smoothing(flux, sigma_value=2):
    smoothed_flux = gaussian_filter1d(input=flux, sigma=sigma_value)
    return smoothed_flux


# Savitzky-Golay filter Function - with defaults

def sav_gol_filter(flux, sav_gol_window=5, sav_gol_poly=2):
    smoothed_flux = savgol_filter(x=flux, window_length=sav_gol_window, polyorder=sav_gol_poly)
    return smoothed_flux


# Spectral rebinning with spectres

def spectres_rebin(flux, wave, desired_res):
    if len(wave) != len(flux):
        print('The lengths of the wave and flux array do not match.\nPlease check and try again')
        print(f'Wave:{len(wave)}\tFlux:{len(flux)}')
        sys.exit()

    rounded_wave_start = round_down_to_nearest(wave[0], desired_res)
    rounded_wave_end = round_up_to_nearest(wave[-1], desired_res)

    lin_space_value = int(((rounded_wave_end - rounded_wave_start) / desired_res) + 1)
    rebinned_wave = np.linspace(rounded_wave_start, rounded_wave_end, lin_space_value)

    rebinned_flux = spectres.spectres(
        new_wavs=rebinned_wave,
        spec_wavs=wave,
        spec_fluxes=flux,
        fill=0
    )

    return rebinned_flux, rebinned_wave


#############################
# Line Fitting Functions
#############################

#############################################
# Gaussian Fitting Function
#############################################

def single_gaussian(x, height, center, width):
    result = height * np.exp(-(x - center) ** 2 / (2 * width ** 2))
    return result


#############################################
# Gaussian EQW Function
#############################################

# Due to the uncertainty about this equation - it's not used
# Gaussian EQWs are calculated the same way as the narrow features

def gaussian_eqw(amp, sigma):
    eqw_gaussian = (amp * sigma) / 0.3989  # Not 100% that this equation is right

    return eqw_gaussian


#############################################
# STDDEV to FWHM Conversion
#############################################

def stddev_to_fwhm(list):
    result = np.ndarray.tolist(np.array(list) * (2 * ((2 * np.log(2)) ** 0.5)))

    return result


def stddev_to_fwhm_non_list(value):
    result = value * (2 * ((2 * np.log(2)) ** 0.5))

    return result

#############################################
# FWHM to STDDEV Conversion
#############################################

def fwhm_to_stddev(list):
    result = np.ndarray.tolist(np.array(list) / (2 * ((2 * np.log(2)) ** 0.5)))

    return result


def fwhm_to_stddev_non_list(value):
    result = value / (2 * ((2 * np.log(2)) ** 0.5))

    return result

#############################################
# FWHM to Velocity
#############################################

"""
Line center should be given in the same units as fwhm is measured in
Velocity will be returned in kms^-1
"""

def fwhm_to_velocity(fwhm, line_center):
    velocity = c * (fwhm / line_center)
    return velocity

#############################################
# Velocity to FWHM
#############################################

"""
Line center should be given in the same units as fwhm is measured in
Velocity will be returned in kms^-1
"""

def velocity_to_fwhm(velocity, line_center):
    fwhm = (velocity / c ) * line_center
    return fwhm



########################################################################
# Megaparsec --> cm and visa versa conversion
########################################################################

def mpc_to_cm(mpc):
    cm = mpc * 3.08567758128e+24
    return cm


def cm_to_mpc(cm):
    mpc = cm / 3.08567758128e+24
    return mpc


########################################################################
# Flux to luminosity and visa versa
########################################################################

def flux_to_luminosity(flux, distance):
    luminosity = flux * ((4 * math.pi) * distance ** 2)
    # Doesn't automatically check units - distance should most likely be in cm
    return luminosity


def luminosity_to_flux(luminosity, distance):
    flux = luminosity / ((4 * math.pi) * distance ** 2)
    # Doesn't automatically check units - distance should most likely be in cm
    return flux


#############################################
# Flux to luminosity conversion
#############################################

def flux_to_log_line_luminosity(flux, distance):
    # Distance is expected to be in Mpc
    distance_cm = mpc_to_cm(distance)
    line_lum = flux_to_luminosity(flux, distance_cm)
    log_line_lum = np.log10(line_lum)
    # The log version is returned
    return log_line_lum


#############################################
# Luminosity to flux conversion
#############################################

def log_line_lum_to_flux(log_line_lum, distance):
    # Distance is expected to be in Mpc and the luminosity in the log10 version
    distance_cm = mpc_to_cm(distance)
    line_lum = 10 ** log_line_lum
    flux = luminosity_to_flux(line_lum, distance_cm)
    # The log version is returned
    return flux


########################
# Lick Index Calculation
#
# Python translated version of the IDL script by Decker French
#######################

# Measure D4000 and H-delta A

def lick_index_calculation(wave, flux, err, verbose=False):
    s = {
        "d4000_n": 0.0,
        "d4000_n_err": -1.0,
        "lick_hd_a": 0.0,
        "lick_hd_a_err": -1.0,
        "lick_hg_a": 0.0,
        "lick_hg_a_err": -1.0
    }

    # D4000_N (Balogh)
    wl_b1 = 3850.0
    wl_b2 = 3950.0
    wl_b = []
    flux_b = []
    err_b = []

    wl_r1 = 4000.0
    wl_r2 = 4100.0
    wl_r = []
    flux_r = []
    err_r = []

    if np.nanmin(wave) < wl_b1 and np.nanmax(wave) > wl_r2:

        for ii, item in enumerate(wave):
            if wl_b1 < wave[ii] < wl_b2:
                wl_b.append(wave[ii])
                flux_b.append(flux[ii])
                err_b.append(err[ii])
            elif wl_r1 < wave[ii] < wl_r2:
                wl_r.append(wave[ii])
                flux_r.append(flux[ii])
                err_r.append(err[ii])

        wl_b = np.array(wl_b)
        flux_b = np.array(flux_b)
        err_b = np.array(err_b)

        wl_r = np.array(wl_r)
        flux_r = np.array(flux_r)
        err_r = np.array(err_r)

        fl_b = np.trapz(x=wl_b, y=flux_b * wl_b * wl_b)
        fl_r = np.trapz(x=wl_r, y=flux_r * wl_r * wl_r)

        fl_b_err = np.sqrt(np.trapz(x=wl_b, y=(err_b * wl_b * wl_b) ** 2))
        fl_r_err = np.sqrt(np.trapz(x=wl_r, y=(err_r * wl_r * wl_r) ** 2))

        d4000_n = (wl_b2 - wl_b1) / (wl_r2 - wl_r1) * fl_r / fl_b
        d4000_n_err = (wl_b2 - wl_b1) / (wl_r2 - wl_r1) * np.sqrt((fl_r_err / fl_b) ** 2 +
                                                                  (fl_b_err * fl_r / fl_b ** 2) ** 2)
        s.update({'d4000_n': d4000_n, 'd4000_n_err': d4000_n_err})
    else:
        if verbose:
            print("Insufficient wavelength coverage for D4000")
        s.update({'d4000_n': -999, 'd4000_n_err': -999})

    # Loop over H-delta and H-gamma
    index = 0
    while index <= 1:

        if index == 0:
            # Lick H-delta A
            wli = [4083.500, 4122.250]
            wlb = [4041.600, 4079.750]
            wlr = [4128.500, 4161.000]

        elif index == 1:
            # Lick H-gamma A
            wli = [4319.750, 4363.500]
            wlb = [4283.500, 4319.750]
            wlr = [4367.250, 4419.750]

        else:
            print("Script Failure\nLick index loop overrunning\nExiting")
            sys.exit()

        if np.nanmin(wave) < wlb[0] and np.nanmax(wave) > wlr[1]:

            wave_blue_cont = []
            flux_blue_cont = []
            err_blue_cont = []

            wave_red_cont = []
            flux_red_cont = []
            err_red_cont = []

            wave_lick = []
            flux_lick = []
            err_lick = []

            for ii, item in enumerate(wave):

                if wlb[0] < wave[ii] < wlb[1]:
                    wave_blue_cont.append(wave[ii])
                    flux_blue_cont.append(flux[ii])
                    err_blue_cont.append(err[ii])

                if wlb[0] < wave[ii] < wlb[1]:
                    wave_red_cont.append(wave[ii])
                    flux_red_cont.append(flux[ii])
                    err_red_cont.append(err[ii])

                if wli[0] < wave[ii] < wli[1]:
                    wave_lick.append(wave[ii])
                    flux_lick.append(flux[ii])
                    err_lick.append(err[ii])

            wave_blue_cont = np.array(wave_blue_cont)
            flux_blue_cont = np.array(flux_blue_cont)
            err_blue_cont = np.array(err_blue_cont)

            wave_red_cont = np.array(wave_red_cont)
            flux_red_cont = np.array(flux_red_cont)
            err_red_cont = np.array(err_red_cont)

            wave_lick = np.array(wave_lick)
            flux_lick = np.array(flux_lick)
            err_lick = np.array(err_lick)

            f_blue_cont = np.trapz(x=wave_blue_cont, y=flux_blue_cont) / (wlb[1] - wlb[0])
            err_blue_cont = np.sqrt(np.trapz(x=wave_blue_cont, y=err_blue_cont ** 2) / (wlb[1] - wlb[0]))
            wl_blue = (wlb[0] + wlb[1]) / 2.0

            f_red_cont = np.trapz(x=wave_red_cont, y=flux_red_cont) / (wlr[1] - wlr[0])
            err_red_cont = np.sqrt(np.trapz(x=wave_red_cont, y=err_red_cont ** 2) / (wlr[1] - wlr[0]))
            wl_red = (wlr[0] + wlr[1]) / 2.0

            # Linear fit across continuum windows
            m = (f_red_cont - f_blue_cont) / (wl_red - wl_blue)
            m_err = np.sqrt(err_blue_cont ** 2 + err_red_cont ** 2) / (wl_red - wl_blue)
            b = f_blue_cont
            b_err = err_blue_cont

            f_cont = m * (wave_lick - wl_blue) + b
            err_cont = np.sqrt((m_err * (wave_lick - wl_blue)) ** 2 + b_err ** 2)

            index_value = np.trapz(x=wave_lick, y=(1 - np.array(flux_lick) / f_cont))

            nerri = (np.array(flux_lick) / f_cont) * np.sqrt((err_lick / np.array(flux_lick)) ** 2
                                                             + (err_cont / f_cont) ** 2)

            index_err = np.sqrt(np.trapz(x=wave_lick, y=nerri ** 2))

            if index == 0:
                s.update({'lick_hd_a': index_value, 'lick_hd_a_err': index_err})

            if index == 1:
                s.update({'lick_hg_a': index_value, 'lick_hg_a_err': index_err})

        else:
            if index == 0:
                if verbose:
                    print("Wavelength coverage insufficient for Lick Hdelta index calculation. \nExiting")

                s.update({'lick_hd_a': -999, 'lick_hd_a_err': -999})
            elif index == 1:
                if verbose:
                    print("Wavelength coverage insufficient for Lick Hgamma index calculation. \nExiting")

                s.update({'lick_hg_a': -999, 'lick_hg_a_err': -999})
            else:
                print("Script Failure\nLick index loop wavelength check overrunning\nExiting")
                sys.exit()

        index += 1

    return s


########################################################################
# Rounding to nearest multiple
########################################################################

def round_up_to_nearest(x, a):
    try:
        rounded = math.ceil(x / a) * a
    except ValueError:
        print('Rounding up error encountered')
        rounded = -999

    return rounded


def round_down_to_nearest(x, a):
    try:
        rounded = math.floor(x / a) * a
    except ValueError:
        print('Rounding down error encountered')
        rounded = -999

    return rounded


########################################################################
# Find index of value in array closest to given value
########################################################################
def find_nearest_index_only(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


# Actual value as well as the index

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()

    return array[idx], idx


########################
# Directory creation function
########################

def directory_creation(Path):
    try:
        os.makedirs(Path)
    except OSError:
        pass


########################################################################
# mySQL Functions
########################################################################


def database_connection(user, password, database):
    print(f"Now attempting to connect to the database: {database}")

    try:
        connection = mysql.connector.connect(user=user, password=password,
                                             database=database)
    except mysql.connector.Error as err:

        connection = 'Failed'

        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print()
            print("Something is wrong with your user name or password!")
            print()
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print()
            print("Database does not exist!")
            print()
        else:
            print(err)

    print("Connection Complete\n")

    return connection


def commit(connection):
    connection.commit()
    # print("\nmySQL Database update has been completed.")


def commit_and_close(connection, cursor):
    connection.commit()
    cursor.close()
    connection.close()

    # print("\nmySQL Database update has been completed.\nThe connection is now closed.")


def no_commit_and_close(connection, cursor):
    cursor.close()
    connection.close()

    # print("\nThe connection to the mySQL Database is now closed.\nNo update has been saved.")
