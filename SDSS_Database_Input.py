import sys
import time
import Hirogen_Functions

import numpy as np
import astropy.units as u

from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord

startT = time.time()

np.seterr(invalid='raise')

##########################
# Configuration
##########################

User = 'Joe'
Input_Data_File = "confirmed_spectra.csv"
TableID = "SDSS_Confirmed_Spectra"

# Input_Data_File = "MyResult_2021615_RA_Above_20_Dec_Above_10_25000.csv"
# Input_Data_File = "MyResult_2021615_RA_Above_20_Dec_Below_10_25000.csv"

# Input_Data_File = "MyResult_202176_RA_Above_20_Dec_Above_10_50000.csv"
# Input_Data_File = "MyResult_2021621_RA_Above_20_Dec_Below_10_50000.csv"

# TableID = "SDSS_Test_Sample"

# TableID = "SDSS_Test_QSO_Sample"
# Input_Data_File = "All_QSOs_Less_Than_Point2.csv"
# Input_Data_File = "All_QSOs_Less_Than_Point4.csv"

# Input_Data_File = "SDSS_CL_QSOs_Updated.csv"
# TableID = "SDSS_CL_QSOs"

User_Config = Hirogen_Functions.user_config(user=User)

Database = User_Config[0]
Database_User = User_Config[1]
Database_Password = User_Config[2]

# SDSS Extinction Conversion Constants
u_extinction_constant = 5.155
g_extinction_constant = 3.793
r_extinction_constant = 2.751
i_extinction_constant = 2.086
z_extinction_constant = 1.479


# These are from the 'Early Data Release Paper 2002' and are the conversion factors from E(B-V) to per-band extinction


def sdss_dr16_photom_datafile_reader(file):
    """Processes the casjobs dr16 photometry query output file into a format capable of being added to the database"""

    spec_data = np.genfromtxt(file, unpack=True, comments='#', dtype=str, delimiter=",", skip_header=1)

    # Photometric ID
    col0 = spec_data[0].astype(int)

    # RA + DEC
    col1 = (np.array(spec_data[1])).astype(float)
    col2 = (np.array(spec_data[2])).astype(float)

    # Parent ID
    col3 = spec_data[3].astype(int)

    # Spectroscopic ID

    """
    col5 = []
    for ii, item in enumerate(spec_data[4]):
        try:
            col5.append(spec_data[4][ii].astype(int))
            print(len(spec_data[4][ii]))
        except OverflowError:
            print(ii)
            print(spec_data[4][ii])
            print(len(spec_data[4][ii]))
    """

    col4 = np.array(spec_data[4])

    # survey
    col5 = spec_data[5]

    # run2d
    col6 = spec_data[6]

    # Field
    col7 = spec_data[7].astype(int)

    # Plate
    col8 = spec_data[8]

    # MJD
    col9 = spec_data[9].astype(int)

    # fiberID
    col10 = spec_data[10]

    # spec z
    col11 = (np.array(spec_data[11])).astype(float)

    # spec z err
    col12 = (np.array(spec_data[12])).astype(float)

    # spec z warning
    col13 = spec_data[13].astype(int)

    # spec Median SNR
    col14 = (np.array(spec_data[14])).astype(float)

    # SDSS 'type'
    col15 = spec_data[15].astype(int)

    # SDSS Spectral Classification
    col16 = spec_data[16]

    # SDSS Spectral Sub-Classification
    col17 = spec_data[17]

    # SDSS Spectral Human Comments
    col18 = spec_data[18]

    # Extinction u
    col19 = (np.array(spec_data[19])).astype(float)
    # Extinction g
    col20 = (np.array(spec_data[20])).astype(float)
    # Extinction r
    col21 = (np.array(spec_data[21])).astype(float)
    # Extinction i
    col22 = (np.array(spec_data[22])).astype(float)
    # Extinction z
    col23 = (np.array(spec_data[23])).astype(float)

    # Clean flag
    col24 = spec_data[24].astype(int)

    # Phot flags, phot_mjd, phot_skyVersion, phot_run, phot_rerun
    col25 = (np.array(spec_data[25])).astype(int)
    col26 = (np.array(spec_data[26])).astype(int)
    col27 = (np.array(spec_data[27])).astype(int)
    col28 = (np.array(spec_data[28])).astype(int)
    col29 = (np.array(spec_data[29])).astype(int)

    # Petrosian u mag data: Mag, Err, Radius, Radius Err
    col30 = (np.array(spec_data[30])).astype(float)
    col31 = (np.array(spec_data[31])).astype(float)
    col32 = (np.array(spec_data[32])).astype(float)
    col33 = (np.array(spec_data[33])).astype(float)

    # Petrosian g mag data: Mag, Err, Radius, Radius Err
    col34 = (np.array(spec_data[34])).astype(float)
    col35 = (np.array(spec_data[35])).astype(float)
    col36 = (np.array(spec_data[36])).astype(float)
    col37 = (np.array(spec_data[37])).astype(float)

    # Petrosian r mag data: Mag, Err, Radius, Radius Err
    col38 = (np.array(spec_data[38])).astype(float)
    col39 = (np.array(spec_data[39])).astype(float)
    col40 = (np.array(spec_data[40])).astype(float)
    col41 = (np.array(spec_data[41])).astype(float)

    # Petrosian i mag data: Mag, Err, Radius, Radius Err
    col42 = (np.array(spec_data[42])).astype(float)
    col43 = (np.array(spec_data[43])).astype(float)
    col44 = (np.array(spec_data[44])).astype(float)
    col45 = (np.array(spec_data[45])).astype(float)

    # Petrosian z mag data: Mag, Err, Radius, Radius Err
    col46 = (np.array(spec_data[46])).astype(float)
    col47 = (np.array(spec_data[47])).astype(float)
    col48 = (np.array(spec_data[48])).astype(float)
    col49 = (np.array(spec_data[49])).astype(float)

    # Model u mag data: Mag, Err
    col50 = (np.array(spec_data[50])).astype(float)
    col51 = (np.array(spec_data[51])).astype(float)

    # Model g mag data: Mag, Err
    col52 = (np.array(spec_data[52])).astype(float)
    col53 = (np.array(spec_data[53])).astype(float)

    # Model r mag data: Mag, Err
    col54 = (np.array(spec_data[54])).astype(float)
    col55 = (np.array(spec_data[55])).astype(float)

    # Model i mag data: Mag, Err
    col56 = (np.array(spec_data[56])).astype(float)
    col57 = (np.array(spec_data[57])).astype(float)

    # Model z mag data: Mag, Err
    col58 = (np.array(spec_data[58])).astype(float)
    col59 = (np.array(spec_data[59])).astype(float)

    return \
        col0, col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, \
        col17, col18, col19, col20, col21, col22, col23, col24, col25, col26, col27, col28, col29, col30, col31, \
        col32, col33, col34, col35, col36, col37, col38, col39, col40, col41, col42, col43, col44, col45, col46, \
        col47, col48, col49, col50, col51, col52, col53, col54, col55, col56, col57, col58, col59


debug = False
long_mode_distance_calc_mode = False
# This is a slower but more robust method of dealing with the distance calculation - useful for catching mathematical
# errors e.g. nans if the initial pull down from SDSS or DESI hasn't caught them all

##########################
# Cosmology Configuration
##########################

cosmology = FlatLambdaCDM(H0=73 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.27)
# Matches the LOSS rates cosmology - 70 kms might be better though as the LOSS papers are rather old at this point

##########################
# Other Configuration
##########################
np.set_printoptions(threshold=np.inf)
SDSS_Explore_url = 'http://skyserver.sdss.org/dr16/en/tools/explore/Summary.aspx?id='
SDSS_Navigate_url = 'http://skyserver.sdss.org/dr16/en/tools/chart/navi.aspx?ra='

Data_Input = sdss_dr16_photom_datafile_reader(f"{Input_Data_File}")

if debug:
    print(Data_Input)

DR16_PhotometricID = Data_Input[0]
RA = Data_Input[1]
DEC = Data_Input[2]

Co_Ords = SkyCoord(np.array(RA), np.array(DEC), unit='deg')
Co_Ord_Strings = Co_Ords.to_string('hmsdms', precision=4)

RA_HMS = []
DEC_DMS = []
SDSS_Shortnames = []

for ii, item in enumerate(Co_Ord_Strings):
    RA_HMS.append(Co_Ord_Strings[ii][0:14])
    DEC_DMS.append(Co_Ord_Strings[ii][15:])
    SDSS_Shortnames.append(f'SDSS J{Co_Ord_Strings[ii][0:2]}{Co_Ord_Strings[ii][3:5]}{Co_Ord_Strings[ii][15]}'
                           f'{Co_Ord_Strings[ii][16:18]}{Co_Ord_Strings[ii][19:21]}')

DR16_ParentID = Data_Input[3]
DR16_SpectroscopicID = Data_Input[4]

Survey = Data_Input[5]
Run2D = Data_Input[6]

Field = Data_Input[7]
Plate = Data_Input[8]
MJD = Data_Input[9]
FiberID = Data_Input[10]

#######
# Adds the leading zeros to the plate and fiber ids so they match those of the spectral files
#######
for ii, item in enumerate(Plate):
    if len(str(Plate[ii])) < 4:
        Plate[ii] = Plate[ii].zfill(4)
    if len(str(FiberID[ii])) < 4:
        FiberID[ii] = FiberID[ii].zfill(4)

Spec_z = Data_Input[11]
Spec_z_err = Data_Input[12]
Spec_z_warning = Data_Input[13]

if long_mode_distance_calc_mode:

    Distance = []
    MaxDistance = []
    MinDistance = []

    for ii, item in enumerate(Spec_z):
        try:
            Distance_Value = cosmology.luminosity_distance(Spec_z[ii]).value
        except FloatingPointError:
            print('Distance Calculation Failure')
            print(Spec_z[ii])
            print(Spec_z_err[ii])
            print(f'{SDSS_Explore_url}{DR16_PhotometricID[ii]}')
            Distance_Value = np.nan

        try:
            Max_Distance_Value = cosmology.luminosity_distance(Spec_z[ii] + Spec_z_err[ii]).value
        except FloatingPointError:
            print('Max Distance Calculation Failure')
            print(Spec_z[ii])
            print(Spec_z_err[ii])
            print(f'{SDSS_Explore_url}{DR16_PhotometricID[ii]}')
            Max_Distance_Value = np.nan

        try:
            Min_Distance_Value = cosmology.luminosity_distance(Spec_z[ii] - Spec_z_err[ii]).value
        except FloatingPointError:
            print('Min Distance Calculation Failure')
            print(Spec_z[ii])
            print(Spec_z_err[ii])
            print(f'{SDSS_Explore_url}{DR16_PhotometricID[ii]}')
            Min_Distance_Value = np.nan

        Distance.append(Distance_Value)
        MaxDistance.append(Max_Distance_Value)
        MinDistance.append(Min_Distance_Value)

    Distance = np.array(Distance)
    MaxDistance = np.array(MaxDistance)
    MinDistance = np.array(MinDistance)

else:

    Distance = cosmology.luminosity_distance(Spec_z).value
    MaxDistance = cosmology.luminosity_distance(Spec_z + Spec_z_err).value
    MinDistance = cosmology.luminosity_distance(Spec_z - Spec_z_err).value

Max_Distance_Uncertainty = MaxDistance - Distance
Min_Distance_Uncertainty = Distance - MinDistance

Average_Distance_Uncertainty = (Max_Distance_Uncertainty + Min_Distance_Uncertainty) / 2

Distance_Modulus = []
Distance_Modulus_Err = []

for ii, item in enumerate(DR16_SpectroscopicID):

    try:
        Distance_Modulus.append(5 * np.log10(Distance[ii] * 1e6) - 5)
    except FloatingPointError:
        Distance_Modulus.append(0)

    try:

        Max_Distance_Modulus = 5 * np.log10(MaxDistance[ii] * 1e6) - 5
        Min_Distance_Modulus = 5 * np.log10(MinDistance[ii] * 1e6) - 5

        Distance_Modulus_Err.append((Max_Distance_Modulus - Min_Distance_Modulus) / 2)

    except FloatingPointError:
        Distance_Modulus_Err.append(0)

Median_SNR = Data_Input[14]

SDSS_Type = Data_Input[15]
Spec_Classification = Data_Input[16]
Spec_Sub_Classification = Data_Input[17]

Spec_Human_Comments = Data_Input[18]

Extinction_u = Data_Input[19]
Extinction_g = Data_Input[20]
Extinction_r = Data_Input[21]
Extinction_i = Data_Input[22]
Extinction_z = Data_Input[23]

E_BminusV_u = Extinction_u / u_extinction_constant
E_BminusV_g = Extinction_g / g_extinction_constant
E_BminusV_r = Extinction_r / r_extinction_constant
E_BminusV_i = Extinction_i / i_extinction_constant
E_BminusV_z = Extinction_z / z_extinction_constant

Mean_E_BminusV = np.mean([E_BminusV_u, E_BminusV_g, E_BminusV_r, E_BminusV_i, E_BminusV_z], axis=0)

Clean_Flag = Data_Input[24]

Phot_Flags = Data_Input[25]
Phot_MJD = Data_Input[26]
Phot_skyVersion = Data_Input[27]
Phot_Run = Data_Input[28]
Photo_Rerun = Data_Input[29]

# Note: Per the standard for the general photometric database, extinction correction is now applied only to the
# absolute magnitude values. Colour values are however corrected

Petrosian_u = Data_Input[30]
Petrosian_u_err = Data_Input[31]
Petrosian_u_rad = Data_Input[32]
Petrosian_u_rad_err = Data_Input[33]
Petrosian_u_AB = Petrosian_u - Distance_Modulus - Extinction_u
Petrosian_u_AB_err = (Petrosian_u_err ** 2 + (0.434 * Average_Distance_Uncertainty / Distance) ** 2) ** 0.5

Petrosian_g = Data_Input[34]
Petrosian_g_err = Data_Input[35]
Petrosian_g_rad = Data_Input[36]
Petrosian_g_rad_err = Data_Input[37]
Petrosian_g_AB = Petrosian_g - Distance_Modulus - Extinction_g
Petrosian_g_AB_err = (Petrosian_g_err ** 2 + (0.434 * Average_Distance_Uncertainty / Distance) ** 2) ** 0.5

Petrosian_r = Data_Input[38]
Petrosian_r_err = Data_Input[39]
Petrosian_r_rad = Data_Input[40]
Petrosian_r_rad_err = Data_Input[41]
Petrosian_r_AB = Petrosian_r - Distance_Modulus - Extinction_r
Petrosian_r_AB_err = (Petrosian_r_err ** 2 + (0.434 * Average_Distance_Uncertainty / Distance) ** 2) ** 0.5

Petrosian_i = Data_Input[42]
Petrosian_i_err = Data_Input[43]
Petrosian_i_rad = Data_Input[44]
Petrosian_i_rad_err = Data_Input[45]
Petrosian_i_AB = Petrosian_i - Distance_Modulus - Extinction_i
Petrosian_i_AB_err = (Petrosian_i_err ** 2 + (0.434 * Average_Distance_Uncertainty / Distance) ** 2) ** 0.5

Petrosian_z = Data_Input[46]
Petrosian_z_err = Data_Input[47]
Petrosian_z_rad = Data_Input[48]
Petrosian_z_rad_err = Data_Input[49]
Petrosian_z_AB = Petrosian_z - Distance_Modulus - Extinction_z
Petrosian_z_AB_err = (Petrosian_z_err ** 2 + (0.434 * Average_Distance_Uncertainty / Distance) ** 2) ** 0.5

Model_u = Data_Input[50]
Model_u_err = Data_Input[51]
Model_u_AB = Model_u - Distance_Modulus - Extinction_u
Model_u_AB_err = (Model_u_err ** 2 + (0.434 * Average_Distance_Uncertainty / Distance) ** 2) ** 0.5

Model_g = Data_Input[52]
Model_g_err = Data_Input[53]
Model_g_AB = Model_g - Distance_Modulus - Extinction_g
Model_g_AB_err = (Model_g_err ** 2 + (0.434 * Average_Distance_Uncertainty / Distance) ** 2) ** 0.5

Model_r = Data_Input[54]
Model_r_err = Data_Input[55]
Model_r_AB = Model_r - Distance_Modulus - Extinction_r
Model_r_AB_err = (Model_r_err ** 2 + (0.434 * Average_Distance_Uncertainty / Distance) ** 2) ** 0.5

Model_i = Data_Input[56]
Model_i_err = Data_Input[57]
Model_i_AB = Model_i - Distance_Modulus - Extinction_i
Model_i_AB_err = (Model_i_err ** 2 + (0.434 * Average_Distance_Uncertainty / Distance) ** 2) ** 0.5

Model_z = Data_Input[58]
Model_z_err = Data_Input[59]
Model_z_AB = Model_z - Distance_Modulus - Extinction_z
Model_z_AB_err = (Model_z_err ** 2 + (0.434 * Average_Distance_Uncertainty / Distance) ** 2) ** 0.5

Petrosian_u_minus_r = Petrosian_u_AB - Petrosian_r_AB
Petrosian_u_minus_r_observation_err = (Petrosian_u_AB_err ** 2 + Petrosian_r_AB_err ** 2) ** 0.5

Petrosian_g_minus_r = Petrosian_g - Petrosian_r
Petrosian_g_minus_r_observation_err = (Petrosian_g_AB_err ** 2 + Petrosian_r_AB_err ** 2) ** 0.5

Petrosian_g_minus_i = Petrosian_g - Petrosian_i
Petrosian_g_minus_i_observation_err = (Petrosian_g_AB_err ** 2 + Petrosian_i_AB_err ** 2) ** 0.5

Model_u_minus_r = Model_u - Model_r
Model_u_minus_r_observation_err = (Model_u_AB_err ** 2 + Model_r_AB_err ** 2) ** 0.5

Model_g_minus_r = Model_g - Model_r
Model_g_minus_r_observation_err = (Model_g_AB_err ** 2 + Model_r_AB_err ** 2) ** 0.5

Model_g_minus_i = Model_g - Model_i
Model_g_minus_i_observation_err = (Model_g_AB_err ** 2 + Model_i_AB_err ** 2) ** 0.5

TotalHostMassEstimate_Petro = 1.15 + 0.7 * Petrosian_g_minus_i - 0.4 * Petrosian_i_AB
TotalHostMassEstimate_Model = 1.15 + 0.7 * Model_g_minus_i - 0.4 * Model_i_AB

Mass_Systematic_Error = 0.1  # Source: Taylor2011 - This is added in quadrature to the colour mass error
Mass_Colour_Error_Petro = ((Petrosian_g_minus_i_observation_err ** 2) + (Petrosian_i_AB_err ** 2)) ** 0.5
Mass_Colour_Error_Model = ((Model_g_minus_i_observation_err ** 2) + (Model_i_AB_err ** 2)) ** 0.5

TotalHostMassEstimate_Err_Petro = (Mass_Systematic_Error ** 2 + Mass_Colour_Error_Petro ** 2) ** 0.5
TotalHostMassEstimate_Err_Model = (Mass_Systematic_Error ** 2 + Mass_Colour_Error_Model ** 2) ** 0.5

# Petrosian magnitudes should be preferred for this calculation but there is at least one case in the confirmed
# candidates where the petrosian radii and magnitudes have truly biblical uncertainties so I'm adding this calculation
# as a background and for comparison purposes

Petro_AB_i_Check = []
AGN_Check = []

SDSS_Explore_Links = []
SDSS_Navigate_Links = []

for ii, item in enumerate(DR16_SpectroscopicID):

    if -21.3 < Petrosian_i_AB[ii] < 18.8:
        Petro_AB_i_Check.append(1)
    else:
        Petro_AB_i_Check.append(0)

    SDSS_Explore_Links.append(f'{SDSS_Explore_url}{DR16_PhotometricID[ii]}')
    SDSS_Navigate_Links.append(f'{SDSS_Navigate_url}{RA[ii]}&dec={DEC[ii]}')

Data = Hirogen_Functions.database_connection(user=Database_User, password=Database_Password, database=Database)
cursor = Data.cursor()
cursor.execute(f"SELECT * FROM {Database}.{TableID}")
rows = cursor.fetchall()

print(f'Number of rows prior to row entry is: {len(rows)}')

if debug:
    print(rows)

# Inserting new rows if required
for ii, item in enumerate(DR16_SpectroscopicID):
    cmd = f"INSERT IGNORE INTO `{Database}`.`{TableID}` (DR16_Spectroscopic_ID) " \
          f"VALUES ({DR16_SpectroscopicID[ii]});"
    cursor.execute(cmd)

cursor.execute(f"SELECT DR16_Spectroscopic_ID FROM {Database}.{TableID}")
rows = cursor.fetchall()
print(f'Number of rows is now: {len(rows)}')

if debug:
    print(rows)

########################
# Generation of the files needed to download the spectra in 1000 object blocks + storage of data within the database
########################

print("Adding new object data to the database and generating the files required to download the new spectra")

Spec_to_Download_ID = 0
Download_File_Number = 1

while Spec_to_Download_ID < len(DR16_SpectroscopicID):

    Spectra_Download_File = open(f"Spectra_to_Download_{Input_Data_File}_{Download_File_Number}.txt", 'w+')

    while Spec_to_Download_ID < len(DR16_SpectroscopicID):

        cmd = f"update {Database}.{TableID} " \
              f"set Right_Ascension = '{RA[Spec_to_Download_ID]}'," \
              f"Declination = '{DEC[Spec_to_Download_ID]}'," \
              f"RA_HMS = '{RA_HMS[Spec_to_Download_ID]}'," \
              f"DEC_DMS = '{DEC_DMS[Spec_to_Download_ID]}'," \
              f"SDSS_ShortName = '{SDSS_Shortnames[Spec_to_Download_ID]}'," \
              f"DR16_ParentID = '{DR16_ParentID[Spec_to_Download_ID]}'," \
              f"DR16_Photometric_ID = '{DR16_PhotometricID[Spec_to_Download_ID]}'," \
              f"SDSS_DR16_Explore_Link ='{SDSS_Explore_Links[Spec_to_Download_ID]}'," \
              f"SDSS_DR16_Navigate_Link ='{SDSS_Navigate_Links[Spec_to_Download_ID]}'," \
              f"survey = '{Survey[Spec_to_Download_ID]}'," \
              f"run2d = '{Run2D[Spec_to_Download_ID]}'," \
              f"field = '{Field[Spec_to_Download_ID]}'," \
              f"spec_Plate = '{Plate[Spec_to_Download_ID]}'," \
              f"spec_MJD = '{MJD[Spec_to_Download_ID]}' ," \
              f"spec_FiberID = '{FiberID[Spec_to_Download_ID]}'," \
              f"SDSS_Type = '{SDSS_Type[Spec_to_Download_ID]}'," \
              f"spec_class_SDSS = '{Spec_Classification[Spec_to_Download_ID]}'," \
              f"spec_subclass_SDSS = '{Spec_Sub_Classification[Spec_to_Download_ID]}'," \
              f"spec_Human_Comments = '{Spec_Human_Comments[Spec_to_Download_ID]}'," \
              f"z_SDSS_spec = '{Spec_z[Spec_to_Download_ID]}'," \
              f"z_err_SDSS_spec = '{Spec_z_err[Spec_to_Download_ID]}'," \
              f"z_warning_SDSS_spec = '{Spec_z_warning[Spec_to_Download_ID]}'," \
              f"Distance_MPC = '{Distance[Spec_to_Download_ID]}'," \
              f"Max_Distance_MPC = '{MaxDistance[Spec_to_Download_ID]}'," \
              f"Min_Distance_MPC = '{MinDistance[Spec_to_Download_ID]}'," \
              f"Distance_Modulus='{Distance_Modulus[Spec_to_Download_ID]}'," \
              f"Distance_Modulus_Err='{Distance_Modulus_Err[Spec_to_Download_ID]}'," \
              f"median_SNR_SDSS_spec = '{Median_SNR[Spec_to_Download_ID]}'," \
              f"u_extinction = '{Extinction_u[Spec_to_Download_ID]}'," \
              f"g_extinction = '{Extinction_g[Spec_to_Download_ID]}'," \
              f"r_extinction = '{Extinction_r[Spec_to_Download_ID]}'," \
              f"i_extinction = '{Extinction_i[Spec_to_Download_ID]}'," \
              f"z_extinction = '{Extinction_z[Spec_to_Download_ID]}'," \
              f"generalised_extinction = '{Mean_E_BminusV[Spec_to_Download_ID]}'," \
              f"u_petro_mag = '{Petrosian_u[Spec_to_Download_ID]}'," \
              f"u_petro_mag_err = '{Petrosian_u_err[Spec_to_Download_ID]}'," \
              f"u_petro_AB_mag = '{Petrosian_u_AB[Spec_to_Download_ID]}'," \
              f"u_petro_AB_mag_err = '{Petrosian_u_AB_err[Spec_to_Download_ID]}'," \
              f"u_petro_rad= '{Petrosian_u_rad[Spec_to_Download_ID]}'," \
              f"u_petro_rad_err= '{Petrosian_u_rad_err[Spec_to_Download_ID]}'," \
              f"g_petro_mag = '{Petrosian_g[Spec_to_Download_ID]}'," \
              f"g_petro_mag_err = '{Petrosian_g_err[Spec_to_Download_ID]}'," \
              f"g_petro_AB_mag = '{Petrosian_g_AB[Spec_to_Download_ID]}'," \
              f"g_petro_AB_mag_err = '{Petrosian_g_AB_err[Spec_to_Download_ID]}'," \
              f"g_petro_rad= '{Petrosian_g_rad[Spec_to_Download_ID]}'," \
              f"g_petro_rad_err= '{Petrosian_g_rad_err[Spec_to_Download_ID]}'," \
              f"r_petro_mag = '{Petrosian_r[Spec_to_Download_ID]}'," \
              f"r_petro_mag_err = '{Petrosian_r_err[Spec_to_Download_ID]}'," \
              f"r_petro_AB_mag = '{Petrosian_r_AB[Spec_to_Download_ID]}'," \
              f"r_petro_AB_mag_err = '{Petrosian_r_AB_err[Spec_to_Download_ID]}'," \
              f"r_petro_rad= '{Petrosian_r_rad[Spec_to_Download_ID]}'," \
              f"r_petro_rad_err= '{Petrosian_r_rad_err[Spec_to_Download_ID]}'," \
              f"i_petro_mag = '{Petrosian_i[Spec_to_Download_ID]}'," \
              f"i_petro_mag_err = '{Petrosian_i_err[Spec_to_Download_ID]}'," \
              f"i_petro_AB_mag = '{Petrosian_i_AB[Spec_to_Download_ID]}'," \
              f"i_petro_AB_mag_err = '{Petrosian_i_AB_err[Spec_to_Download_ID]}'," \
              f"i_petro_rad= '{Petrosian_i_rad[Spec_to_Download_ID]}'," \
              f"i_petro_rad_err= '{Petrosian_i_rad_err[Spec_to_Download_ID]}'," \
              f"z_petro_mag = '{Petrosian_z[Spec_to_Download_ID]}'," \
              f"z_petro_mag_err = '{Petrosian_z_err[Spec_to_Download_ID]}'," \
              f"z_petro_AB_mag = '{Petrosian_z_AB[Spec_to_Download_ID]}'," \
              f"z_petro_AB_mag_err = '{Petrosian_z_AB_err[Spec_to_Download_ID]}'," \
              f"z_petro_rad= '{Petrosian_z_rad[Spec_to_Download_ID]}'," \
              f"z_petro_rad_err= '{Petrosian_z_rad_err[Spec_to_Download_ID]}'," \
              f"Petro_Total_Host_Mass_Estimate = '{TotalHostMassEstimate_Petro[Spec_to_Download_ID]}'," \
              f"Petro_Total_Host_Mass_Estimate_err = '{TotalHostMassEstimate_Err_Petro[Spec_to_Download_ID]}'," \
              f"Model_Total_Host_Mass_Estimate = '{TotalHostMassEstimate_Model[Spec_to_Download_ID]}'," \
              f"Model_Total_Host_Mass_Estimate_err = '{TotalHostMassEstimate_Err_Model[Spec_to_Download_ID]}'," \
              f"u_model_mag = '{Model_u[Spec_to_Download_ID]}'," \
              f"u_model_mag_err = '{Model_u_err[Spec_to_Download_ID]}'," \
              f"u_model_AB_mag = '{Model_u_AB[Spec_to_Download_ID]}'," \
              f"u_model_AB_mag_err = '{Model_u_AB_err[Spec_to_Download_ID]}'," \
              f"g_model_mag = '{Model_g[Spec_to_Download_ID]}'," \
              f"g_model_mag_err = '{Model_g_err[Spec_to_Download_ID]}'," \
              f"g_model_AB_mag = '{Model_g_AB[Spec_to_Download_ID]}'," \
              f"g_model_AB_mag_err = '{Model_g_AB_err[Spec_to_Download_ID]}'," \
              f"r_model_mag = '{Model_r[Spec_to_Download_ID]}'," \
              f"r_model_mag_err = '{Model_r_err[Spec_to_Download_ID]}'," \
              f"r_model_AB_mag = '{Model_r_AB[Spec_to_Download_ID]}'," \
              f"r_model_AB_mag_err = '{Model_r_AB_err[Spec_to_Download_ID]}'," \
              f"i_model_mag = '{Model_i[Spec_to_Download_ID]}'," \
              f"i_model_mag_err = '{Model_i_err[Spec_to_Download_ID]}'," \
              f"i_model_AB_mag = '{Model_i_AB[Spec_to_Download_ID]}'," \
              f"i_model_AB_mag_err = '{Model_i_AB_err[Spec_to_Download_ID]}'," \
              f"z_model_mag = '{Model_z[Spec_to_Download_ID]}'," \
              f"z_model_mag_err = '{Model_z_err[Spec_to_Download_ID]}'," \
              f"z_model_AB_mag = '{Model_z_AB[Spec_to_Download_ID]}'," \
              f"z_model_AB_mag_err = '{Model_z_AB_err[Spec_to_Download_ID]}'," \
              f"u_minus_r_petro = '{Petrosian_u_minus_r[Spec_to_Download_ID]}'," \
              f"u_minus_r_petro_observation_err = '{Petrosian_u_minus_r_observation_err[Spec_to_Download_ID]}'," \
              f"g_minus_r_petro = '{Petrosian_g_minus_r[Spec_to_Download_ID]}'," \
              f"g_minus_r_petro_observation_err = '{Petrosian_g_minus_r_observation_err[Spec_to_Download_ID]}'," \
              f"g_minus_i_petro = '{Petrosian_g_minus_i[Spec_to_Download_ID]}'," \
              f"g_minus_i_petro_observation_err = '{Petrosian_g_minus_i_observation_err[Spec_to_Download_ID]}'," \
              f"u_minus_r_model = '{Model_u_minus_r[Spec_to_Download_ID]}'," \
              f"u_minus_r_model_observation_err = '{Model_u_minus_r_observation_err[Spec_to_Download_ID]}'," \
              f"g_minus_r_model = '{Model_g_minus_r[Spec_to_Download_ID]}'," \
              f"g_minus_r_model_observation_err = '{Model_g_minus_r_observation_err[Spec_to_Download_ID]}'," \
              f"g_minus_i_model = '{Model_g_minus_i[Spec_to_Download_ID]}'," \
              f"g_minus_i_model_observation_err = '{Model_g_minus_i_observation_err[Spec_to_Download_ID]}'," \
              f"clean = '{Clean_Flag[Spec_to_Download_ID]}'," \
              f"flags = '{Phot_Flags[Spec_to_Download_ID]}'," \
              f"Phot_MJD = '{Phot_MJD[Spec_to_Download_ID]}'," \
              f"Phot_skyVersion = '{Phot_skyVersion[Spec_to_Download_ID]}'," \
              f"Phot_Run = '{Phot_Run[Spec_to_Download_ID]}',"\
              f"Phot_Rerun = '{Photo_Rerun[Spec_to_Download_ID]}'," \
              f"ABS_Petro_i_Expected_Mag_Range = '{Petro_AB_i_Check[Spec_to_Download_ID]}'" \
              f" where DR16_Spectroscopic_ID = '{DR16_SpectroscopicID[Spec_to_Download_ID]}' "

        if Survey[Spec_to_Download_ID] in ['sdss', 'segue1', 'segue2']:

            Spectra_Download_File.write(
                f"/sdss/spectro/redux/{Run2D[Spec_to_Download_ID]}/spectra/{Plate[Spec_to_Download_ID]}/spec-{Plate[Spec_to_Download_ID]}-"
                f"{MJD[Spec_to_Download_ID]}-{FiberID[Spec_to_Download_ID]}.fits\n"
            )

        elif Survey[Spec_to_Download_ID] in ['eboss', 'boss']:

            Spectra_Download_File.write(
                f"/eboss/spectro/redux/{Run2D[Spec_to_Download_ID]}/spectra/full/{Plate[Spec_to_Download_ID]}/spec-{Plate[Spec_to_Download_ID]}-"
                f"{MJD[Spec_to_Download_ID]}-{FiberID[Spec_to_Download_ID]}.fits\n"
            )

        else:
            print("Not sure how to handle this survey - exiting for now")
            print(f"{DR16_SpectroscopicID[Spec_to_Download_ID]}\t{Survey[Spec_to_Download_ID]}")
            sys.exit()

        cursor.execute(cmd)
        Spec_to_Download_ID = Spec_to_Download_ID + 1

        if Spec_to_Download_ID % 2500 == 0:
            Download_File_Number = Download_File_Number + 1
            Spectra_Download_File.close()
            break

Spectra_Download_File.close()

print("New object data added.\nClosing connection to database and exiting.")

Hirogen_Functions.commit_and_close(Data, cursor)

endT = time.time()
executionT = endT - startT

sys.exit()
