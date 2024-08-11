""" Practica 2-e
    para el punto e) se agrega la corrección del reloj del satélite
    4ta. columna de las coord. precisas de cada satélite
"""
# pylint: disable= C0103, C0206, C0301, C0209, W0105, W0602, W0603, W0621

import os
import platform

from math import sqrt
from random import random
from numpy import transpose, linalg

c = 299792458         # m/s  de ITRF

L = []
A = []
C = []


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

Estacion = [          # Coord. precisas de la Estación [m]
    3370658.6942,     # X
    711877.0150,      # Y
    5349786.8637]     # Z


# Para tomar como coordenadas a-priori de la estacion.
#    se le agrega una diferencia a las precisas
Coord = [(i+(random()-0.5)*5000) for i in Estacion]
#Coord = Estacion

# en la segunda iteracion el error de reloj va a estar estimado,
# para lo cual se va a necesitar agregar un elemento
# al vector posición
Coord.append(0)


def arma_matriz():
    """     al armar la matriz de diseño se puede poner c * Delta_t en la
            4ta. columna en vez de C, así los resultados 
            dan en Distancia, en lugar de Delta_t.
            Para eso, pongo todos unos en vez de c,
            así la incógnita incluye a c y quedan números más manejables
    """
    global L
    global A
    global C
    L = []
    A = []
    C = []
    j = 0
    for st in Precisas:
        s= Precisas[st]
        dX = s[0] * 1000 - Coord[0]
        dY = s[1] * 1000 - Coord[1]
        dZ = s[2] * 1000 - Coord[2]
        ρ = sqrt(dX*dX+dY*dY+dZ*dZ)
        fila = [dX/ρ,dY/ρ,dZ/ρ,1]  # opcion con incognita c * Delta_t

        A.append(fila)
        # diferencia Observado - Calculado
        err = PD[st] - ρ + float(Coord[3]) + c * s[3] / 1E6  # Si s[3] > 0 el satélite atrasa con respecto a GPS time, entonces "sumo" error en distancia?
        L.append(err)
        linea_C =[0 for i in range(cant_sat)]
        #linea_C[j]= ρ - float(Coord[3]) + c * s[3] / 1E6 - sqrt(Coord[0]*Coord[0]+Coord[1]*Coord[1]+Coord[2]*Coord[2])
        linea_C[j]= err * err
        j +=1
        C.append(linea_C)

def imprime_resu():
    """Imprime resultados parciales"""

    print("\nObservado- calculado\n")
    for dif in L:
        print("{:15.6f}  ".format(dif))
    print()

    """
    for linea in P:
        print(linea)
    print()
    """

    print("\nMatriz de diseño\n")

    for i in range(cant_sat):
        linea =""
        for i in A[i]:
            if i == int(i):
                linea += "{:5d}  ".format(i)
            else:
                linea += "{:20.16f}  ".format(i)
        print(linea)

    print("\n")
    """
    print()
    print("Coord Calculada")
    linea =""
    for i in Coord:
        linea += "{:20.16f}  ".format(i)
    print(linea)
    print()
    """
    print(" - - - - - - - - - -")
    print("Diferencias Calculadas para corregir las coordenadas")
    linea =""
    for i in X1:
        linea += "{:20.16f}  ".format(i)
    print(linea)
    dift = list(X1[:3])
    print("\nModulo de la diferencia calculada: {:20.16f} m".format(linalg.norm(dift)))
    dift.append(X1[3]*c)
    print("Modulo de la diferencia calculada              incluyendo el reloj: {:10.6f} m\n".format(linalg.norm(dift)))
    print(" - - - - - - - - - -\n")

def imprime_Correg():
    """ Imprime una vez recalculado"""
    print("Coord Corregida,    Precisa,                Diferencia (Calculada - Precisa)")
    j = 0
    acu = 0
    for i in Coord:
        if j<3:
            linea = "{:15.5f}  {:15.5f}  -->{:15.5f} m".format(i,Estacion[j],i-Estacion[j])
            acu += (i-Estacion[j])*(i-Estacion[j])
        else:
            linea = "\n Dif. entre sitios: {:15.5f} m\n".format(sqrt(acu))
            linea += " Delta_t:           {:15.5f} useg\n\n".format(i/c*1E6)
        j +=1
        print(linea)

############################### -  -  -  -  -  -  -  -  -  -  -  -  -
#   I N I C I O    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
############################### -  -  -  -  -  -  -  -  -  -  -  -  -

if platform.system() == "Linux":
    # import readchar # type: ignore
    os.system('clear')
elif platform.system() == "Windows":
    # import msvcrt   # type: ignore
    os.system('cls')


imprime_Correg()   # Primero muestra la condición inicial desde donde partimos
for paso in range(3):
    print("--------------------------------------------------------")
    print("----> Paso: {:4d}".format(paso))
    print()
    arma_matriz()
    P = linalg.inv(C)
    #X1 = linalg.inv(transpose(A) @ P @ A) @ transpose(A) @ P @ L
    X1 = linalg.inv(transpose(A) @ A) @ transpose(A) @ L
    imprime_resu()
    Coord=[(Coord[i] - X1[i]) for i in range( len(Coord))]  # a more 'pythonic' way
    imprime_Correg()

    #Cxyz = linalg.inv(transpose(A) @ P @ A)
    #print(Cxyz)
