# pylint: disable= C0209, C0301, W0611

import os
import platform
#from datetime import datetime
import pandas as pd
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
print(tabla)


# Hago la diferencia
tabla["P2-C1"]=tabla["lpg2-27-P2"]-tabla["lpg2-27-C1"]

# Puedo ver cual es el valor medio o constante de la diferencia . . .
valor_med = tabla["P2-C1"].mean()
# y armar otra columan sin el valor constante . . .
#tabla["P2-C1_debiased"] = tabla["P2-C1"]-valor_med

print()
print(tabla.head())
print('\n')

time_axis = list(tabla["Seconds"].astype('datetime64[us]'))
dif_1 = tabla["P2-C1"].tolist()

fig = go.Figure([go.Scatter(x=time_axis, y=dif_1,mode="markers")],layout=go.Layout(
        title=go.layout.Title(text="TP 3-1")))
fig.update_xaxes(title_text="Horas")
fig.update_yaxes(title_text="P2-C1")
fig.show()

"""
print("x - Exit")
if platform.system() == "Linux":
    respuesta = readchar.readchar().decode('utf-8')
elif platform.system() == "Windows":
    respuesta = msvcrt.getch().decode('utf-8')
"""
