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



# Hago la diferencia 
sat7["P2-C1"]=sat7["P2"]-sat7["C1"]

# Puedo ver cual es el valor medio o constante de la diferencia . . .
valor_med = sat7["P2-C1"].mean()
# y armar otra columan sin el valor constante . . .
sat7["P2-C1_debiased"] = sat7["P2-C1"]-valor_med


# Hago la diferencia Phase
sat7["L1L2"]= c / f0 * (sat7["L1"]/154-sat7["L2"]/120)   # metros

"""
time_axis = list(sat7["FechaHora"].astype('datetime64[us]'))
dif_1 = sat7["P2-C1"].tolist()
dif_2 = sat7["L1L2"].tolist()

fig = go.Figure([go.Scatter(x=time_axis, y=dif_2,mode="markers")],layout=go.Layout(
        title=go.layout.Title(text="TP 3-1")))
fig.update_xaxes(title_text="Horas")
fig.update_yaxes(title_text="P2-C1")
fig.show()

"""

print("Armando nuevos dataframes con los bloques encontrados . . .")


G1 = sat7[sat7["Grupo"] == 1]
# Olvidarse del grupo 2 que es muy corto...
G3 = sat7[sat7["Grupo"] == 3]
G4 = sat7[sat7["Grupo"] == 4]
G5 = sat7[sat7["Grupo"] == 5]


print("Finalizado !\n")


# ------------------------------------------------------------------------------------
# Ploteando . . . . .
#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
print("Ploteando . . .")

fig, ax = plt.subplots(4, gridspec_kw={'hspace': 0.4})#, sharex=False, sharey=False, gridspec_kw={'hspace': 0})  # nuevo
fig.set_size_inches(9, 6)   # w , h

#plt.rcParams['axes.grid'] = True   # nuevo

formatter = ticker.FormatStrFormatter('%1.2f')
tformatter = mdates.DateFormatter('%H:%M')

dotsize = 0.7

ax[0].plot(G1["FechaHora"].astype('datetime64[us]'), G1["L1L2"], 'tab:red')
ax1 = ax[0].twinx()
ax1.scatter(G1["FechaHora"].astype('datetime64[us]'), G1["P2-C1"], c='green', s=dotsize)
ax1.tick_params(axis="y", labelsize=7)

ax[1].plot(G3["FechaHora"].astype('datetime64[us]'), G3["L1L2"], 'tab:green')
ax1 = ax[1].twinx()
ax1.scatter(G3["FechaHora"].astype('datetime64[us]'), G3["P2-C1"], c='orange', s=dotsize)
ax1.tick_params(axis="y", labelsize=7)

ax[2].plot(G4["FechaHora"].astype('datetime64[us]'), G4["L1L2"], 'tab:blue')
ax[2].set(ylabel="L1L2 [m]")
ax1 = ax[2].twinx()
ax1.scatter(G4["FechaHora"].astype('datetime64[us]'), G4["P2-C1"], c='red', s=dotsize)
ax1.tick_params(axis="y", labelsize=7)
ax1.set(ylabel='P2 - C1 [m]')

ax[3].plot(G5["FechaHora"].astype('datetime64[us]'), G5["L1L2"], 'tab:orange')
ax1 = ax[3].twinx()
ax1.scatter(G5["FechaHora"].astype('datetime64[us]'), G5["P2-C1"], c='blue', s=dotsize)
ax1.tick_params(axis="y", labelsize=7)

for i in range(4): 
    #ax[i].set(ylabel='L1L2 [m]')
    ax[i].tick_params(axis="x", labelsize=7)
    ax[i].tick_params(axis="y", labelsize=6)
    ax[i].yaxis.set_major_formatter(formatter)
    ax[i].xaxis.set_major_formatter(tformatter)

fig.suptitle("TP 3 - 2y3:   L1L2 vs P2-C1 [m]",fontsize=13)

plt.show()
