# pylint: disable= C0209, C0301, W0611

import os
import platform
from datetime import datetime
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
# Columnas: "FechaHora", "C1", "P2", "L1", "L2", "D1", "D2"
# Pero dejamos sin importar los registros que no tienen valores de P2 (que tampoco tienen L2, etc.)
#
sat7 = (pd.read_csv("Practica_3\\sat7.csv", parse_dates=True, dayfirst=True)).query("P2 > 1")

# Elimino algunos puntos anómalos . . . .
sat7 = sat7.drop(sat7.loc[sat7["FechaHora"] == "24/05/2002 01:37:00.002000"].index)
sat7 = sat7.drop(sat7.loc[sat7["FechaHora"] == "24/05/2002 12:00:00.016000"].index)
sat7 = sat7.drop(sat7.loc[sat7["FechaHora"] == "24/05/2002 12:01:00.016000"].index)

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
