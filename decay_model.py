import enum
import sys
import Hirogen_Functions

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from statistics import mean

NERSC_Flag = False # True - Running on NERSC, FALSE - Running locally
Local_DESI = False
# Both flags should be False to run on local SDSS spectra

if not NERSC_Flag:

    User = 'Joe'
    User_Config = Hirogen_Functions.user_config(user=User)

    #TableID = "SDSS_Test_QSO_Sample"
    #TableID = "SDSS_Test_Sample"
    #TableID = "SDSS_Galaxy_Spectra"
    TableID = "SDSS_ECLE_Spectra"
    #TableID = 'DESI_Local'

    # Search for 'QUERY' to find the actual database access location where parameters can be adjusted

    Database = User_Config[0]
    Database_User = User_Config[1]
    Database_Password = User_Config[2]
    Main_Spectra_Path = User_Config[3]

########################
# Object Info - From mySQL Database
########################

Data = Hirogen_Functions.database_connection(user=Database_User, password=Database_Password, database=Database)

print("\n")
cursor = Data.cursor()

very_shortnames = ["'SDSS J0748'", "'SDSS J0938'", "'SDSS J0952'", "'SDSS J1055'", "'SDSS J1241'", "'SDSS J1342'", "'SDSS J1350'"]

very_shortname = very_shortnames[0]

invalid_lines = []

# QUERY

cursor.execute(
    f"SELECT spec_MJD, lin_con_pEQW_FeVII6008, lin_con_pEQW_FeX6376, "
    f"lin_con_pEQW_FeXI7894, lin_con_pEQW_FeXIV5304, SDSS_Very_ShortName, lin_con_pEQW_OIII5007 "
    f"FROM `{Database}`.`{TableID}`"
    #f"WHERE Manually_Inspected_Flag != -10 AND lin_con_LineFlux_Halpha is NULL"
    #f"WHERE lin_con_pEQW_Halpha is NULL"
    #f"WHERE z_SDSS_spec >= 0.2"
    f"WHERE SDSS_Very_ShortName = {very_shortname} AND Follow_Up_ID = 0"
    #f"WHERE Follow_Up_ID = 0"
    #f"WHERE Follow_Up_ID in (1,2,3)"
    #f"WHERE Nickname = 'Leafeon' AND Follow_Up_ID in (2,3)"
    #f"WHERE Standard_Inclusion = 1"
    #f"WHERE DR16_Spectroscopic_ID IN {IDs} and Standard_Inclusion = 1"
    #f"WHERE DR16_Spectroscopic_ID = {IDs}"
)

Candidate_Data = cursor.fetchall()
# Candidate_Data = cursor.fetchmany(200)

#print(Candidate_Data)

if len(Candidate_Data) >= 1:

    og_MJD_List = np.array([item[0] for item in Candidate_Data])
    og_Fe_pEQW_List = [[item[1] for item in Candidate_Data], [item[2] for item in Candidate_Data], [item[3] for item in Candidate_Data], [item[4] for item in Candidate_Data]]
    og_Short_Name = np.array([item[5] for item in Candidate_Data])
    og_OIII_pEQW = np.array([item[6] for item in Candidate_Data])

else:
    print("Candidate_Data Length error: Check and try again.")

    sys.exit()

clean_Og_Fe_List = []

for i, line in enumerate(og_Fe_pEQW_List):
    clean_Og_Fe_List.append(line[0])
    if line[0] <= -999 and i not in invalid_lines:
        invalid_lines.append(i)

clean_Og_Fe_Array = np.array(clean_Og_Fe_List)

"""for i, eqw in enumerate(clean_Og_Fe_Array):
    if eqw > 0:
        clean_Og_Fe_Array[i] = 0"""

clean_Og_Fe_List.append(np.mean(clean_Og_Fe_Array[clean_Og_Fe_Array != -999]))
#clean_Og_Fe_List.append(np.mean(clean_Og_Fe_Array[np.logical_and(clean_Og_Fe_Array != -999, clean_Og_Fe_Array < 0)]))

fin_Og_MJD = og_MJD_List - (og_MJD_List-1)


cursor.execute(
    f"SELECT spec_MJD, lin_con_pEQW_FeVII6008, lin_con_pEQW_FeX6376, "
    f"lin_con_pEQW_FeXI7894, lin_con_pEQW_FeXIV5304, SDSS_Very_ShortName, lin_con_pEQW_OIII5007 "
    f"FROM `{Database}`.`{TableID}`"
    #f"WHERE Manually_Inspected_Flag != -10 AND lin_con_LineFlux_Halpha is NULL"
    #f"WHERE lin_con_pEQW_Halpha is NULL"
    #f"WHERE z_SDSS_spec >= 0.2"
    f"WHERE SDSS_Very_ShortName = {very_shortname} and Follow_Up_ID = 1"
    #f"WHERE Follow_Up_ID = 0"
    #f"WHERE Follow_Up_ID in (1,2,3)"
    #f"WHERE Nickname = 'Leafeon' AND Follow_Up_ID in (2,3)"
    #f"WHERE Standard_Inclusion = 1"
    #f"WHERE DR16_Spectroscopic_ID IN {IDs} and Standard_Inclusion = 1"
    #f"WHERE DR16_Spectroscopic_ID = {IDs}"
)

Candidate_Data = cursor.fetchall()
# Candidate_Data = cursor.fetchmany(200)

#print(Candidate_Data)

if len(Candidate_Data) >= 1:

    middle_MJD_List = np.array([item[0] for item in Candidate_Data])
    middle_Fe_pEQW_List = [[item[1] for item in Candidate_Data], [item[2] for item in Candidate_Data], [item[3] for item in Candidate_Data], [item[4] for item in Candidate_Data]]
    mid_Short_Name = np.array([item[5] for item in Candidate_Data])
    mid_OIII_pEQW = np.array([item[6] for item in Candidate_Data])

else:
    print("Candidate_Data Length error: Check and try again.")

    sys.exit()

clean_Mid_Fe_List = []

for i, line in enumerate(middle_Fe_pEQW_List):
    clean_Mid_Fe_List.append(line[0])
    if line[0] <= -999 and i not in invalid_lines:
        invalid_lines.append(i)

clean_Mid_Fe_Array = np.array(clean_Mid_Fe_List)

"""for i, eqw in enumerate(clean_Mid_Fe_Array):
    if eqw > 0:
        clean_Mid_Fe_Array[i] = 0"""

clean_Mid_Fe_List.append(np.mean(clean_Mid_Fe_Array[clean_Mid_Fe_Array != -999]))
#clean_Mid_Fe_List.append(np.mean(clean_Mid_Fe_Array[np.logical_and(clean_Mid_Fe_Array != -999, clean_Mid_Fe_Array < 0)]))

fin_Mid_MJD = middle_MJD_List - (og_MJD_List-1)


cursor.execute(
    f"SELECT spec_MJD, lin_con_pEQW_FeVII6008, lin_con_pEQW_FeX6376, "
    f"lin_con_pEQW_FeXI7894, lin_con_pEQW_FeXIV5304, SDSS_Very_ShortName, lin_con_pEQW_OIII5007 "
    f"FROM `{Database}`.`{TableID}`"
    #f"WHERE Manually_Inspected_Flag != -10 AND lin_con_LineFlux_Halpha is NULL"
    #f"WHERE lin_con_pEQW_Halpha is NULL"
    #f"WHERE z_SDSS_spec >= 0.2"
    f"WHERE SDSS_Very_ShortName = {very_shortname} and Follow_Up_ID = 2"
    #f"WHERE Follow_Up_ID = 0"
    #f"WHERE Follow_Up_ID in (1,2,3)"
    #f"WHERE Nickname = 'Leafeon' AND Follow_Up_ID in (2,3)"
    #f"WHERE Standard_Inclusion = 1"
    #f"WHERE DR16_Spectroscopic_ID IN {IDs} and Standard_Inclusion = 1"
    #f"WHERE DR16_Spectroscopic_ID = {IDs}"
)

Candidate_Data = cursor.fetchall()
# Candidate_Data = cursor.fetchmany(200)

#print(Candidate_Data)

if len(Candidate_Data) >= 1:

    new_MJD_List = np.array([item[0] for item in Candidate_Data])
    new_Fe_pEQW_List = [[item[1] for item in Candidate_Data], [item[2] for item in Candidate_Data], [item[3] for item in Candidate_Data], [item[4] for item in Candidate_Data]]
    new_Short_Name = np.array([item[5] for item in Candidate_Data])
    new_OIII_pEQW = np.array([item[6] for item in Candidate_Data])

else:
    print("Candidate_Data Length error: Check and try again.")

    sys.exit()

clean_New_Fe_List = []

for i, line in enumerate(new_Fe_pEQW_List):
    clean_New_Fe_List.append(line[0])
    if line[0] <= -999 and i not in invalid_lines:
        invalid_lines.append(i)

clean_New_Fe_Array = np.array(clean_New_Fe_List)

"""for i, eqw in enumerate(clean_New_Fe_Array):
    if eqw > 0:
        clean_New_Fe_Array[i] = 0"""

clean_New_Fe_List.append(np.mean(clean_New_Fe_Array[clean_New_Fe_Array != -999]))
#clean_New_Fe_List.append(np.mean(clean_New_Fe_Array[np.logical_and(clean_New_Fe_Array != -999, clean_New_Fe_Array < 0)]))

fin_New_MJD = new_MJD_List - (og_MJD_List-1)


def sigmoid(x, A, K, B, v, Q):
    return A+((K-A)/((1+Q*np.exp(-B*x))**(1/v)))

def power_law(t, M, N):
    return -M*t**(-N)


"""fig, axes = plt.subplots(2,3, figsize=(15,10), sharey = True)

fig.suptitle(f"{og_Short_Name[0]}")

axes[1][2].set_visible(False)

axes[1][0].set_position([0.24,0.125,0.228,0.343])
axes[1][1].set_position([0.55,0.125,0.228,0.343])

l, k = 0, 0

titles = ['FeVII', 'FeX', 'FeXI', 'FeXIV', 'Average of Coronal']

MJD_range = np.arange(fin_Og_MJD, fin_New_MJD, 1)


for j, line in enumerate(clean_Og_Fe_List):
    
    times = np.array([fin_Og_MJD[0], fin_Mid_MJD[0], fin_New_MJD[0]])
    EQWs = np.array([clean_Og_Fe_List[j], clean_Mid_Fe_List[j], clean_New_Fe_List[j]])

    axes[k][l].invert_yaxis()
    axes[k][l].set_ylabel(r"EQW, $\AA$")
    axes[k][l].set_xlabel("MJD")
    axes[k][l].set_title(titles[j])

    if j in invalid_lines:
        times = np.array([fin_Og_MJD[0], fin_New_MJD[0]])
        EQWs = np.array([clean_Og_Fe_List[j], clean_New_Fe_List[j]])

    for i, eqw in enumerate(EQWs):
        if eqw > 0:
            EQWs[i] = 0

    par, cov = curve_fit(power_law, times, EQWs, maxfev = 10000, p0 = (1,1))

    np.savetxt("power_params.csv", par)

    EQW_change = power_law(MJD_range, par[0], par[1])

    axes[k][l].text(3000, -4, par[1])

    axes[k][l].plot(times, EQWs, 'o', linestyle = 'None')

    axes[k][l].plot(MJD_range, EQW_change)

    l += 1

    if l == 3:
        k = 1
        l = 0

#axes[1][2].invert_yaxis()
axes[1][2].set_ylabel(r"EQW, $\AA$")
axes[1][2].set_xlabel("MJD")
axes[1][2].set_title("OIII")

times = np.array([fin_Og_MJD[0], fin_Mid_MJD[0], fin_New_MJD[0]])
EQWs = np.array([og_OIII_pEQW[0], mid_OIII_pEQW[0], new_OIII_pEQW[0]])

#axes[1][2].plot(times, EQWs, 'o', linestyle = 'None')

#plt.savefig(f"{og_Short_Name[0]} breakdown.png")"""

MJD_range = np.arange(fin_Og_MJD, fin_New_MJD, 1)

fig, ax = plt.subplots()

ax.invert_yaxis()
ax.set_ylabel(r"EQW, $\AA$")
ax.set_xlabel(f"MJD - {og_MJD_List[0]}")

ax.set_ylim(0, -4)

times = np.array([fin_Og_MJD[0], fin_Mid_MJD[0], fin_New_MJD[0]])
EQWs = np.array([clean_Og_Fe_List[4], clean_Mid_Fe_List[4], clean_New_Fe_List[4]])

par, cov = curve_fit(power_law, times, EQWs, maxfev = 10000, p0 = (1,1))

np.savetxt("power_params.csv", par)

EQW_change = power_law(MJD_range, par[0], par[1])

if very_shortname in "'SDSS J0952', 'SDSS J1241'":
    a, k, b, v, q = np.genfromtxt('FeVII_parameters.csv')
elif very_shortname in "'SDSS J0748', 'SDSS J1342', 'SDSS J1350'":
    a, k, b, v, q = np.genfromtxt('Non_FeVII_parameters.csv')

efficiency = sigmoid(EQW_change, a, k, b, v, q)

evo = efficiency * EQW_change

ax.plot(times, EQWs, 'o', linestyle = 'None')

#ax.plot(MJD_range, EQW_change, 'k')

ax.plot(MJD_range, evo, 'k')

plt.savefig("eff_evo_example.png")

plt.show()