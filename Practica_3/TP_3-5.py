"""
TP 3 - ej 5
calcule las coordenadas del satélite 7 para todas las épocas disponibles
y con las coordenadas de la estación, asigne la elevación del satélite
a cada observación para graficar los pseudo observables del último ítem
de cada ejercicio respecto de la elevación del satélite


  4331297.3380   567555.6400  4633133.7180                  APPROX POSITION XYZ

"""
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
from math import sqrt ,sin ,cos ,asin ,acos ,atan2, pi
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

a = 6378137.0
f = 1 / 298.257223563
e2 = f*(2-f)
e= sqrt(e2)

Estacion = [          # Coord. de la Estación [m]
    4331297.3380,     # X
    567555.6400,      # Y
    4633133.7180]     # Z

est_orbitas = []


def ecef_to_geodetic(x, y, z):
    '''
    Olson, D. K. (1996). Converting Earth-Centered, Earth-Fixed Coordinates to Geodetic Coordinates. 
    IEEE Transactions on Aerospace and Electronic Systems. https://doi.org/10.1109/7.481290

    Converted to Python and modified by Steven Ward.  No rights reserved.
    '''
    w2 = x * x + y * y
    w = sqrt(w2)
    z2 = z * z
    lon_rad = atan2(y, x)

    a1 = a * e2
    a2 = a1 * a1
    a3 = a1 * e2 / 2
    a4 = 2.5 * a2
    a5 = a1 + a3
    #a6 = (1 - self.e2)

    r2 = w2 + z2
    r = sqrt(r2)

    s2 = z2 / r2
    c2 = w2 / r2
    u = a2 / r
    v = a3 - a4 / r

    s = 0
    cc = 0
    ss = 0

    # cos(45°)² == ½
    if c2 > 0.5: # Equatorial
        s = (z / r) * (1 + c2 * (a1 + u + s2 * v) / r)
        lat_rad = asin(s)
        ss = s * s
        cc = sqrt(1 - ss)
    else: # Polar
        cc = (w / r) * (1 - s2 * (a5 - u - c2 * v) / r)
        lat_rad = acos(cc)
        ss = 1 - cc * cc
        s = sqrt(ss)

        if z < 0:
            lat_rad = -lat_rad
            s = -s

    d2 = 1 - e2 * ss
    Rn = a / sqrt(d2)
    Rm = (1 - e2) * Rn / d2
    rf = (1 - e2) * Rn
    u = w - Rn * c
    v = z - rf * s
    ff = c * u + s * v
    m = c * v - s * u
    p = m / (Rm + ff)

    lat_rad += p

    ht = f + m * p / 2

    return (lat_rad, lon_rad, ht)

def procesa():
    global est_orbitas
    est_orbitas = []           # reseteo la lista
    elev = []
    Xest = Estacion[0]
    Yest = Estacion[1]
    Zest = Estacion[2]
    est_orbitas=[]             # abro una lista vacía
    for horalist in orbitas:
        deltaX = float(horalist[1]) - Xest
        deltaY = float(horalist[2]) - Yest
        deltaZ = float(horalist[3]) - Zest
        dist = sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ)
        xprima = R[0][0] * deltaX + R[0][1] * deltaY + R[0][2] * deltaZ
        yprima = R[1][0] * deltaX + R[1][1] * deltaY + R[1][2] * deltaZ
        zprima = R[2][0] * deltaX + R[2][1] * deltaY + R[2][2] * deltaZ
        A = atan2(yprima,xprima)
        Z = acos(zprima/dist)
        E = pi / 2 -Z
        # la lista de la entrada al dict "satname" está compuesta por estas listas,
        #    una para cada entrada de FechaHora --> es decir cada 15'
        est_orbitas.append([horalist[0],xprima,yprima,zprima,dist,A,Z])
        elev.append([horalist[0],E])
    return elev


############################### -  -  -  -  -  -  -  -  -  -  -  -  -
#   I N I C I O    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
############################### -  -  -  -  -  -  -  -  -  -  -  -  -

print("\n Coordenadas de la Estación:")
print(Estacion)
print()


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


# Leo coordenadas del satélite 7
# en los momentos de las observaciones
# y las guardo en la lista "orbitas"
print("Cargando efemérides del satélite 7 . . .")
orbitas = []
with open("Practica_3\\sat7_coord.csv") as forbs:
    l = 0
    for l_obs in forbs:
        if l>0:
            orbitas.append(l_obs[:-1].split(','))
        l +=1
print("Finalizado !\n")

#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
# Calcula las coordenadas gerodéticas para la estación y la matriz de transformación
x = Estacion[0]
y = Estacion[1]
z = Estacion[2]
p = sqrt(x*x+y*y)
r = sqrt(x*x+y*y+z*z)
λ = atan2(y,x)  # declinacion = longitud
φ = atan2(z,p)  # inclinacion
pol_est = [φ,λ,r]
coord = ecef_to_geodetic(x,y,z)
Ҩ = coord[0]  # latitud
lon = coord[1]
h = coord[2]
geo_est = [Ҩ,lon,h]
R = [   [-sin(Ҩ)*cos(λ),-sin(Ҩ)*sin(λ),cos(Ҩ)],
        [-sin(λ)       , cos(λ)       ,0     ],
        [ cos(Ҩ)*cos(λ), cos(Ҩ)*sin(λ),sin(Ҩ)]]


# ------------------------------------------------------------------------------------
#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
# Hago la diferencia del punto 1
sat7["P2-C1"]=sat7["P2"]-sat7["C1"]

# Hago la diferencia del punto 2 (Phase)
sat7["L1L2"]= c / f0 * (sat7["L1"]/154-sat7["L2"]/120)   # metros


#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
# Obtengo las elevaciones a partir de las coordenadas del satélite
print("Calculando elevaciones para los tiempos de las observaciones . . .")
elevaciones = procesa()
print("Finalizado !\n")


# ------------------------------------------------------------------------------------
# Agrego las elevaciones al dataframe
#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
print("Agregando las elevaciones al dataframe . . .")
sat7["Elev"]=0.0
for index, row in sat7.iterrows():
    for HoraElev in elevaciones:
        if row["FechaHora"]==HoraElev[0]:
            sat7.at[index, "Elev"] = HoraElev[1] * 180/ pi
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

# ------------------------------------------------------------------------------------
# Armo nuevos dataframes con los bloques contiguos
#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
print("Armando nuevos dataframes con los bloques encontrados . . .")

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

dotsize1 = 1
dotsize2 = 0.7

ax[0].scatter(G1["FechaHora"].astype('datetime64[us]'), G1["PS3"], c='red', s=dotsize1)
ax1 = ax[0].twinx()
ax1.scatter(G1["FechaHora"].astype('datetime64[us]'), G1["Elev"], c='green', s=dotsize2)
ax1.tick_params(axis="y", labelsize=7)

ax[1].scatter(G3["FechaHora"].astype('datetime64[us]'), G3["PS3"], c='green', s=dotsize1)
ax1 = ax[1].twinx()
ax1.scatter(G3["FechaHora"].astype('datetime64[us]'), G3["Elev"], c='orange', s=dotsize2)
ax1.tick_params(axis="y", labelsize=7)

ax[2].scatter(G4["FechaHora"].astype('datetime64[us]'), G4["PS3"], c='blue', s=dotsize1)
ax[2].set(ylabel="P2-C1 - L1L2 [m]")
ax1 = ax[2].twinx()
ax1.scatter(G4["FechaHora"].astype('datetime64[us]'), G4["Elev"], c='red', s=dotsize2)
ax1.tick_params(axis="y", labelsize=7)
ax1.set(ylabel='Elevacion [º]')

ax[3].scatter(G5["FechaHora"].astype('datetime64[us]'), G5["PS3"], c='orange', s=dotsize1)
ax1 = ax[3].twinx()
ax1.scatter(G5["FechaHora"].astype('datetime64[us]'), G5["Elev"], c='blue', s=dotsize2)
ax1.tick_params(axis="y", labelsize=7)


for i in range(4):
    ax[i].tick_params(axis="x", labelsize=7)
    ax[i].tick_params(axis="y", labelsize=6)
    ax[i].yaxis.set_major_formatter(formatter)
    ax[i].xaxis.set_major_formatter(tformatter)
    ax[i].axhline(0,color='black',linestyle='-')


fig.suptitle("TP 3 - 5:   Combinación de Observables [m] vs Elevaciones [º]",fontsize=13)

plt.show()
print()
