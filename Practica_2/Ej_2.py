
from math import sqrt
from random import random

c = 299792458         # m/s  de ITRF

PD = {                # PseudoDist [m]
"28": 23334619.807,
"13": 22586189.129,
"01": 25167667.280,
"27": 20873169.266,
"24": 23371141.291,
"10": 21505922.486,
"08": 20958408.428
}

cant_sat = len(PD)

Precisas  ={          # Efemérides Precisas [Km] , [us]
    #        X             Y              Z           clk
"01" : [  581.886423, 25616.666528 ,  7088.545471, 169.092800],
"08" : [22018.953984,  2878.718252 , 14451.124018,   9.709180],
"10" : [10103.948910, -10925.429662, 22009.912003,   1.148951],
"13" : [ 7525.432597,  20488.591201, 15216.097471,  -0.655216],
"24" : [22368.646126, -12657.086060,  6934.928617,  36.698468],
"27" : [15057.427636,   9402.947329, 20171.667340,  14.763242],
"28" : [-5895.039751,  14576.928529, 21538.074040,  14.267922]
}

Estacion = [          # Coord Estación [m]
    3370658.6942,     # X
    711877.0150,      # Y
    5349786.8637]     # Z

# Estacion inicial, agregar una diferencia 

Iniciales = [(i+random()*1000) for i in Estacion]
# print(Iniciales)



difs = {}
A = []
def calculo():
    resu ={}
    for st in Precisas:
        s= Precisas[st]
        dX = s[0] * 1000 - Estacion[0]
        dY = s[1] * 1000 - Estacion[1]
        dZ = s[2] * 1000 - Estacion[2]
        ρ = sqrt(dX*dX+dY*dY+dZ*dZ)
        fila = [dX/ρ,dY/ρ,dZ/ρ,c]
        A.append(fila)
        resu[st] = PD[st] - ρ

    return resu

difs = calculo()
"""
print()
for sat in difs:
    print(sat, difs[sat])

print()
for i in range(cant_sat):
    print(A[i])
"""