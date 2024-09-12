# pylint: disable= C0209, C0301, C0413, W0603, W0611

import os
import platform

if platform.system() == "Linux":
    #import readchar # type: ignore
    os.system('clear')
elif platform.system() == "Windows":
    #import msvcrt   # type: ignore
    os.system('cls')
print("\n >> Iniciando . . .")

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

#
# Importamos el archivo .csv:
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Columnas: "FechaHora", "C1", "P2", "L1", "L2", "D1", "D2"
# Pero dejamos sin importar los registros que no tienen valores de P2 (que tampoco tienen L2, etc.)
#
print("Cargando archivo de Observaciones . . .")
sat7 = (pd.read_csv("Practica_3\\sat7.csv", parse_dates=True, dayfirst=True)).query("P2 > 1")

# Elimino algunos puntos anómalos . . . .
sat7 = sat7.drop(sat7.loc[sat7["FechaHora"] == "24/05/2002 01:37:00.002000"].index)
sat7 = sat7.drop(sat7.loc[sat7["FechaHora"] == "24/05/2002 12:00:00.016000"].index)
sat7 = sat7.drop(sat7.loc[sat7["FechaHora"] == "24/05/2002 12:01:00.016000"].index)
print("Finalizado !\n")

# ------------------------------------------------------------------------------------
# Busco los bloques contiguos
#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
print("Buscando los bloques contiguos en el dataframe . . .")
sat7["timestamp"]= pd.to_datetime(sat7["FechaHora"], dayfirst=True)
sat7["gap"] = sat7["timestamp"].diff().dt.seconds > 70

sat7["Grupo"] = 9
grupo = 1
for index, row in sat7.iterrows():
    if row["gap"]:
        grupo +=1
    sat7.at[index, "Grupo"] = grupo
print("Finalizado !\n")
#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -



# Hago la diferencia del punto 1
sat7["P2-C1"]=sat7["P2"]-sat7["C1"]


# Hago la diferencia del punto 2 (Phase)
sat7["L1L2"]= c / f0 * (sat7["L1"]/154-sat7["L2"]/120)   # metros

# ------------------------------------------------------------------------------------
# Armo nuevos dataframes con los bloques contiguos
#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
print("Armando nuevos dataframes con los bloques encontrados . . .")

G1 = sat7[sat7["Grupo"] == 1]
med = G1["L1L2"].mean()
G1["L1L2-deb"] = G1["L1L2"] - med
# max1izq = G1["L1L2-deb"].max()

valor_med = G1["P2-C1"].mean()
G1["P2-C1_debiased"] = G1["P2-C1"]-valor_med
max1der = G1["P2-C1_debiased"].max()
max1izq = max1der
min1der = G1["P2-C1_debiased"].min()
#min1izq = max1izq * min1der/max1der
min1izq = min1der

G1["PS3"] = G1["P2-C1_debiased"] - G1["L1L2-deb"]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Olvidarse del grupo 2 que es muy corto...

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
G3 = sat7[sat7["Grupo"] == 3]
med = G3["L1L2"].mean()
G3["L1L2-deb"] = G3["L1L2"] - med
# max3izq = G3["L1L2-deb"].max()

valor_med = G3["P2-C1"].mean()
G3["P2-C1_debiased"] = G3["P2-C1"]-valor_med
max3der = G3["P2-C1_debiased"].max()
max3izq = max3der
min3der = G3["P2-C1_debiased"].min()
#min3izq = max3izq * min3der/max3der
min3izq = min3der

G3["PS3"] = G3["P2-C1_debiased"] - G3["L1L2-deb"]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
G4 = sat7[sat7["Grupo"] == 4]
med = G4["L1L2"].mean()
G4["L1L2-deb"] = G4["L1L2"] - med

valor_med = G4["P2-C1"].mean()
G4["P2-C1_debiased"] = G4["P2-C1"]-valor_med
max4der = G4["P2-C1_debiased"].max()
max4izq = max4der
min4der = G4["P2-C1_debiased"].min()
min4izq = min4der

G4["PS3"] = G4["P2-C1_debiased"] - G4["L1L2-deb"]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
G5 = sat7[sat7["Grupo"] == 5]
med = G5["L1L2"].mean()
G5["L1L2-deb"] = G5["L1L2"] - med

valor_med = G5["P2-C1"].mean()
G5["P2-C1_debiased"] = G5["P2-C1"]-valor_med
max5der = G5["P2-C1_debiased"].max()
max5izq = max5der
min5der = G5["P2-C1_debiased"].min()
min5izq = min5der

G5["PS3"] = G5["P2-C1_debiased"] - G5["L1L2-deb"]

print("Finalizado !\n")


# ------------------------------------------------------------------------------------
# Ploteando . . . . .
#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
print("Ploteando . . .")

fig, ax = plt.subplots(4, gridspec_kw={'hspace': 0.3})#, sharex=False, sharey=False, gridspec_kw={'hspace': 0})  # nuevo
fig.set_size_inches(10, 7)   # w , h

#plt.rcParams['axes.grid'] = True   # nuevo

formatter = ticker.FormatStrFormatter('%1.2f')
tformatter = mdates.DateFormatter('%H:%M')

dotsize = 0.9
margen = 1.1

ax[0].plot(G1["FechaHora"].astype('datetime64[us]'), G1["L1L2-deb"], 'tab:red')
ax[0].set(ylim=(min1izq * margen, max1izq * margen))
ax1 = ax[0].twinx()
ax1.scatter(G1["FechaHora"].astype('datetime64[us]'), G1["P2-C1_debiased"], c='green', s=dotsize)
ax1.set(ylim=(min1der * margen, max1der * margen))
ax1.tick_params(axis="y", labelsize=7)

ax[1].plot(G3["FechaHora"].astype('datetime64[us]'), G3["L1L2-deb"], 'tab:green')
ax[1].set(ylim=(min3izq * margen, max3izq * margen))
ax1 = ax[1].twinx()
ax1.scatter(G3["FechaHora"].astype('datetime64[us]'), G3["P2-C1_debiased"], c='orange', s=dotsize)
ax1.set(ylim=(min3der * margen, max3der * margen))
ax1.tick_params(axis="y", labelsize=7)

ax[2].plot(G4["FechaHora"].astype('datetime64[us]'), G4["L1L2-deb"], 'tab:blue')
ax[2].set(ylim=(min4izq * margen, max4izq * margen))
ax[2].set(ylabel="L1L2 [m]")
ax1 = ax[2].twinx()
ax1.scatter(G4["FechaHora"].astype('datetime64[us]'), G4["P2-C1_debiased"], c='red', s=dotsize)
ax1.set(ylim=(min4der * margen, max4der * margen))
ax1.tick_params(axis="y", labelsize=7)
ax1.set(ylabel='P2 - C1 [m]')

ax[3].plot(G5["FechaHora"].astype('datetime64[us]'), G5["L1L2-deb"], 'tab:orange')
ax[3].set(ylim=(min5izq * margen, max5izq * margen))
ax1 = ax[3].twinx()
ax1.scatter(G5["FechaHora"].astype('datetime64[us]'), G5["P2-C1_debiased"], c='blue', s=dotsize)
ax1.set(ylim=(min5der * margen, max5der * margen))
ax1.tick_params(axis="y", labelsize=7)

for i in range(4):
    #ax[i].set(ylabel='L1L2 [m]')
    ax[i].tick_params(axis="x", labelsize=7)
    ax[i].tick_params(axis="y", labelsize=6)
    ax[i].yaxis.set_major_formatter(formatter)
    ax[i].xaxis.set_major_formatter(tformatter)

fig.suptitle("TP 3 - 2y3 (V2 constantes removidas):   L1L2 vs P2-C1 [m]",fontsize=13)

plt.show()
print()
