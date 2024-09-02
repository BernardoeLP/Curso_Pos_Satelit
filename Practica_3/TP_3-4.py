# pylint: disable= C0209, C0301, W0611

import os
import platform
from datetime import datetime
import pandas as pd
#import plotly.graph_objects as go
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as mdates

pd.options.mode.chained_assignment = None  # default='warn'

c = 299792458          # m/s  de ITRF
f0 = 10.23E6           # c/s
f1 = 1575.42           # Mc/s
f2 = 1227.60           # Mc/s



if platform.system() == "Linux":
    #import readchar # type: ignore
    os.system('clear')
elif platform.system() == "Windows":
    #import msvcrt   # type: ignore
    os.system('cls')

#
# Importamos el archivo .csv:
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Columnas: "FechaHora", "C1", "P2", "L1", "L2", "D1", "D2"
# Pero dejamos sin importar los registros que no tienen valores de P2 (que tampoco tienen L2, etc.)
#
sat7 = (pd.read_csv("Practica_3\\sat7.csv", parse_dates=True, dayfirst=True)).query("P2 > 1")

# Elimino algunos puntos anómalos . . . .
sat7 = sat7.drop(sat7.loc[sat7["FechaHora"] == "24/05/2002 01:37:00.002000"].index)
sat7 = sat7.drop(sat7.loc[sat7["FechaHora"] == "24/05/2002 12:00:00.016000"].index)
sat7 = sat7.drop(sat7.loc[sat7["FechaHora"] == "24/05/2002 12:01:00.016000"].index)

# ------------------------------------------------------------------------------------
# Busco los bloques contiguos
#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
sat7["timestamp"]= pd.to_datetime(sat7["FechaHora"], dayfirst=True)
sat7["gap"] = sat7["timestamp"].diff().dt.seconds > 70

sat7["Grupo"] = 9
grupo = 1
for index, row in sat7.iterrows():
    if row["gap"]:
        grupo +=1
    sat7.at[index, "Grupo"] = grupo
#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -



# Hago la diferencia
sat7["P2-C1"]=sat7["P2"]-sat7["C1"]

# Puedo ver cual es el valor medio o constante de la diferencia . . .
valor_med = sat7["P2-C1"].mean()
# y armar otra columan sin el valor constante . . .
sat7["P2-C1_debiased"] = sat7["P2-C1"]-valor_med


# Hago la diferencia Phase
sat7["L1L2"]= c / f0 * (sat7["L1"]/154-sat7["L2"]/120)   # metros

print()
print(sat7.head())
print('\n')

# ------------------------------------------------------------------------------------
# Armo nuevos dataframes con los bloques contiguos
#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -

G1 = sat7[sat7["Grupo"] == 1]
med = G1["L1L2"].mean()
G1["L1L2-deb"] = G1["L1L2"] - med

valor_med = G1["P2-C1"].mean()
G1["P2-C1_debiased"] = G1["P2-C1"]-valor_med

G1["PS3"] = G1["P2-C1_debiased"] - G1["L1L2-deb"]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Olvidarse del grupo 2 que es muy corto...

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
G3 = sat7[sat7["Grupo"] == 3]
med = G3["L1L2"].mean()
G3["L1L2-deb"] = G3["L1L2"] - med

valor_med = G3["P2-C1"].mean()
G3["P2-C1_debiased"] = G3["P2-C1"]-valor_med

G3["PS3"] = G3["P2-C1_debiased"] - G3["L1L2-deb"]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
G4 = sat7[sat7["Grupo"] == 4]
med = G4["L1L2"].mean()
G4["L1L2-deb"] = G4["L1L2"] - med

valor_med = G4["P2-C1"].mean()
G4["P2-C1_debiased"] = G4["P2-C1"]-valor_med

G4["PS3"] = G4["P2-C1_debiased"] - G4["L1L2-deb"]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
G5 = sat7[sat7["Grupo"] == 5]
med = G5["L1L2"].mean()
G5["L1L2-deb"] = G5["L1L2"] - med

valor_med = G5["P2-C1"].mean()
G5["P2-C1_debiased"] = G5["P2-C1"]-valor_med

G5["PS3"] = G5["P2-C1_debiased"] - G5["L1L2-deb"]



fig, ax = plt.subplots(4, gridspec_kw={'hspace': 0.4})#, sharex=False, sharey=False, gridspec_kw={'hspace': 0})  # nuevo
fig.set_size_inches(9, 6)   # w , h

#plt.rcParams['axes.grid'] = True   # nuevo

formatter = ticker.FormatStrFormatter('%1.2f')
tformatter = mdates.DateFormatter('%H:%M')

dotsize = 0.7

ax[0].scatter(G1["FechaHora"].astype('datetime64[us]'), G1["PS3"], c='red', s=dotsize)
"""
ax1 = ax[0].twinx()
ax1.scatter(G1["FechaHora"].astype('datetime64[us]'), G1["P2-C1"], c='green', s=dotsize)
ax1.tick_params(axis="y", labelsize=7)
"""
ax[1].scatter(G3["FechaHora"].astype('datetime64[us]'), G3["PS3"], c='green', s=dotsize)
"""
ax1 = ax[1].twinx()
ax1.scatter(G3["FechaHora"].astype('datetime64[us]'), G3["P2-C1"], c='orange', s=dotsize)
ax1.tick_params(axis="y", labelsize=7)
"""
ax[2].scatter(G4["FechaHora"].astype('datetime64[us]'), G4["PS3"], c='blue', s=dotsize)
"""
ax1 = ax[2].twinx()
ax1.scatter(G4["FechaHora"].astype('datetime64[us]'), G4["P2-C1"], c='red', s=dotsize)
ax1.tick_params(axis="y", labelsize=7)
"""
ax[3].scatter(G5["FechaHora"].astype('datetime64[us]'), G5["PS3"], c='orange', s=dotsize)
"""
ax1 = ax[3].twinx()
ax1.scatter(G5["FechaHora"].astype('datetime64[us]'), G5["P2-C1"], c='blue', s=dotsize)
ax1.tick_params(axis="y", labelsize=7)
"""

for i in range(4):
    #ax[i].set(ylabel='L1L2 [m]')
    ax[i].tick_params(axis="x", labelsize=7)
    ax[i].tick_params(axis="y", labelsize=6)
    ax[i].yaxis.set_major_formatter(formatter)
    ax[i].xaxis.set_major_formatter(tformatter)

fig.suptitle("TP 3 - 4:   P2-C1 - L1L2 [m]",fontsize=13)

plt.show()
print()
