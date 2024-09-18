"""
Páctica 4 Ej-1:

Construya las dobles dsaterencias de todos los observables de código (C1, P1 y P2)
combinando lpg2-lpgs y sat02-sat27, sat08-sat02, sat10-sat08 y sat28-sat10.

Construya gráficos para cada observable 
y todas las series de dobles dsaterencias en función del tiempo.
"""
# pylint: disable= C0209, C0301, W0611

import os
import platform
#from datetime import datetime
from numpy import nan
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go

if platform.system() == "Linux":
    #import readchar # type: ignore
    os.system('clear')
elif platform.system() == "Windows":
    #import msvcrt   # type: ignore
    os.system('cls')

#
# Importamos el archivo .csv:
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Columnas:  Seconds,lpg2-27-L1,lpg2-27-L2,lpg2-27-C1,lpg2-27-P1,lpg2-27-P2,lpg2-02-L1,lpg2-02-L2,lpg2-02-C1,lpg2-02-P1,lpg2-02-P2,
#                    lpg2-08-L1,lpg2-08-L2,lpg2-08-C1,lpg2-08-P1,lpg2-08-P2,lpg2-10-L1,lpg2-10-L2,lpg2-10-C1,lpg2-10-P1,lpg2-10-P2,
#                    lpg2-28-L1,lpg2-28-L2,lpg2-28-C1,lpg2-28-P1,lpg2-28-P2,
#                    lpgs-27-L1,lpgs-27-L2,lpgs-27-C1,lpgs-27-P1,lpgs-27-P2,lpgs-02-L1,lpgs-02-L2,lpgs-02-C1,lpgs-02-P1,lpgs-02-P2,
#                    lpgs-08-L1,lpgs-08-L2,lpgs-08-C1,lpgs-08-P1,lpgs-08-P2,lpgs-10-L1,lpgs-10-L2,lpgs-10-C1,lpgs-10-P1,lpgs-10-P2,
#                    lpgs-28-L1,lpgs-28-L2,lpgs-28-C1,lpgs-28-P1,lpgs-28-P2
#
tabla = pd.read_csv("Practica_4\\tabla.csv")
#tabla = tabla.set_index("Seconds")
#print(tabla)

ests = ["lpg2","lpgs"]
obss = ["C1","P1","P2"]
# sat02-sat27, sat08-sat02, sat10-sat08 y sat28-sat10.
dsats = { "DS1" : ["02","27"],
          "DS2" : ["08","02"],
          "DS3" : ["10","08"],
          "DS4" : ["28","10"] }

# Hago las 1º diferencias
# 1º cada código para cada combinación de satélites

ddif = pd.DataFrame()
ddif["Seconds"] = tabla["Seconds"]
for est in ests:
    for dsat in dsats:
        for obs in obss:
            print(est+"-"+dsat+"-"+obs+" = "+dsats[dsat][0]+"-"+dsats[dsat][1]+" en "+obs)
            ddif[est+"-"+dsat+"-"+obs] = tabla[est+"-"+ dsats[dsat][0] +"-"+obs]-tabla[est+"-"+ dsats[dsat][1] +"-"+obs]
        print()
    print()

cont = 1
graficos = []
for dif in dsats:
    for obs in obss:
        columna = "DDIF{:02d}".format(cont)
        columna_desc = ests[0]+"-("+dsats[dif][0]+"-"+dsats[dif][1]+")-"+obs+" menos "+ests[1]+"-("+dsats[dif][0]+"-"+dsats[dif][1]+")-"+obs
        entrada=[columna,columna_desc]
        graficos.append(entrada)
        print(columna + "= " + columna_desc)
        ddif[columna] = ddif[ests[0]+"-"+dif+"-"+obs]-ddif[ests[1]+"-"+dif+"-"+obs]
        cont += 1
    print()

#print(ddif)

# Quitando los "Outliers"
#~~~~~~~~~~~~~~~~~~~~~~~~~

ddif.loc[ddif["DDIF02"].abs() > 10e6, "DDIF02"] = nan
ddif.loc[ddif["DDIF05"].abs() > 10e6, "DDIF05"] = nan
ddif.loc[ddif["DDIF08"].abs() > 10e6, "DDIF08"] = nan
ddif.loc[ddif["DDIF11"].abs() > 10e6, "DDIF11"] = nan

"""


# Puedo ver cual es el valor medio o constante de la diferencia . . .
valor_med = tabla["P2-C1"].mean()
# y armar otra columan sin el valor constante . . .
#tabla["P2-C1_debiased"] = tabla["P2-C1"]-valor_med

print()
print(tabla.head())
print('\n')

"""

time_axis = list(tabla["Seconds"].astype('datetime64[s]'))
#time_axis = list(ddif["Seconds"])


fig = make_subplots(rows=12, cols=1,shared_xaxes=True,vertical_spacing=0.005)

fila = 1
for columns in graficos:
    fig.add_trace(
        go.Scatter(x=time_axis, y=ddif[columns[0]],mode="markers", name=columns[1]),
        row=fila, col=1
    )
    fila += 1
    """
    fig = go.Figure([go.Scatter(x=time_axis, y=dsat_1,mode="markers")],layout=go.Layout(
            title=go.layout.Title(text="TP 3-1")))
    """
    #fig.update_xaxes(title_text="Horas")
    #fig.update_yaxes(title_text=col)

fig.show()
