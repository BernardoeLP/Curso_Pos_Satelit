# pylint: disable= C0209, C0301, W0611

import os
import platform
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as dates



# Columnas: "FechaHora", "C1", "P2", "L1", "L2", "D1", "D2"
sat7 = (pd.read_csv("Practica_3\\sat7.csv", parse_dates=True, dayfirst=True)).query("P2 > 10")
sat7 = sat7.drop(sat7.loc[sat7["FechaHora"] == "24/05/2002 01:37:00.002000"].index)
time_axis = sat7["FechaHora"].tolist()
#sat7 = sat7.set_index("FechaHora")

sat7["C1-P2"]=sat7["C1"]-sat7["P2"]
valor_med = sat7["C1-P2"].mean()
sat7["C1-P2_debiased"] = sat7["C1-P2"]-valor_med
#print('\n')
#print(sat7.head(15))

dif_1 = sat7["C1-P2"].tolist()

# plotting the points
plt.plot(time_axis, dif_1)

# naming the x axis  
plt.xlabel("Fecha y hora")
# naming the y axis
plt.ylabel("C1-P2")

# giving a title to my graph
plt.title("TP_3-1")

# function to show the plot
plt.show()
