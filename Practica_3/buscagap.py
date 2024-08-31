import pandas as pd


sat7 = (pd.read_csv("Practica_3\\sat7.csv", parse_dates=True, dayfirst=True)).query("P2 > 1")

# Elimino algunos puntos anómalos . . . .
sat7 = sat7.drop(sat7.loc[sat7["FechaHora"] == "24/05/2002 01:37:00.002000"].index)
sat7 = sat7.drop(sat7.loc[sat7["FechaHora"] == "24/05/2002 12:00:00.016000"].index)
sat7 = sat7.drop(sat7.loc[sat7["FechaHora"] == "24/05/2002 12:01:00.016000"].index)
#sat7["times"] = sat7["FechaHora"].astype('datetime64[us]')

sat7["timestamp"]= pd.to_datetime(sat7["FechaHora"], dayfirst=True)
sat7["gap"] = sat7["timestamp"].diff().dt.seconds > 70

sat7["Grupo"] = 9
grupo = 1


for index, row in sat7.iterrows():
    if row["gap"]:
        grupo +=1
        print(row["FechaHora"], grupo)
    sat7.at[index, "Grupo"] = grupo

print()
print(sat7.head())
print()
print(sat7.tail())
print('\n')
