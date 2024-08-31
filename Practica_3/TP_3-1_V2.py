# pylint: disable= C0209, C0301, W0611

import os
import platform
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as mdates


if platform.system() == "Linux":
    os.system('clear')
elif platform.system() == "Windows":
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

print()
print(sat7.head())
print('\n')

time_axis = list(sat7["FechaHora"].astype('datetime64[us]'))
dif_1 = sat7["P2-C1"].tolist()


G1 = sat7[sat7["Grupo"] == 1]
# Olvidarse del grupo 2 que es muy corto...
G3 = sat7[sat7["Grupo"] == 3]
G4 = sat7[sat7["Grupo"] == 4]
G5 = sat7[sat7["Grupo"] == 5]



fig, ax = plt.subplots(5, sharex=False, sharey=False, gridspec_kw={'hspace': 0.4})  # nuevo
#fig.set_size_inches(8, 6)   # w , h

#plt.style.use("fivethirtyeight")
#plt.rcParams['axes.grid'] = True   # nuevo

formatter = ticker.FormatStrFormatter('%1.2f')
tformatter = mdates.DateFormatter('%H:%M')
dotsize = 0.7
ax[0].scatter(sat7["FechaHora"].astype('datetime64[us]'), sat7["P2-C1"], c='grey', s=0.5)
ax[1].scatter(G1["FechaHora"].astype('datetime64[us]'), G1["P2-C1"], c='red', s=dotsize)
ax[2].scatter(G3["FechaHora"].astype('datetime64[us]'), G3["P2-C1"], c='green', s=dotsize)
ax[3].scatter(G4["FechaHora"].astype('datetime64[us]'), G4["P2-C1"], c='orange', s=dotsize)
ax[4].scatter(G5["FechaHora"].astype('datetime64[us]'), G5["P2-C1"], c='blue', s=dotsize)

"""
ax[1].plot(G1["FechaHora"].astype('datetime64[us]'), G1["P2-C1"], 'tab:red')
ax[2].plot(G3["FechaHora"].astype('datetime64[us]'), G3["P2-C1"], 'tab:green')
ax[3].plot(G4["FechaHora"].astype('datetime64[us]'), G4["P2-C1"], 'tab:orange')
ax[4].plot(G5["FechaHora"].astype('datetime64[us]'), G5["P2-C1"], 'tab:blue')
"""
for i in range(5):
    #ax[i].set(ylabel='P2-C1 [m]')
    ax[i].tick_params(axis="x", labelsize=7)#,labelrotation=45)
    ax[i].tick_params(axis="y", labelsize=6)#,labelrotation=45)
    ax[i].yaxis.set_major_formatter(formatter)
    ax[i].xaxis.set_major_formatter(tformatter)

fig.suptitle("TP 3 - 1   P2-C1 [m]",fontsize=13)

plt.show()
