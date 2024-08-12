""" Practica 2-e
    para el punto e) se agrega la corrección del reloj del satélite
    4ta. columna de las coord. precisas de cada satélite
"""
# pylint: disable= C0103, C0200, C0206, C0301, C0209, W0105, W0602, W0603, W0621, W0640

import os
import platform

from math import sin, cos, tan, atan, sqrt
from datetime import datetime, timedelta
from random import random
from numpy import  matmul, transpose, linalg

c = 299792458         # m/s  de ITRF
μ = 3.986005E14 # m3/s2  Earth gravitational constant
ωe =   7.2921151467E-5 # radians/s Angular Velocity of the Earth
tGPS0 =  datetime(1980,1,6,0,0,0)

L = []
A = []
C = []
satord=[]        # es un 'parche' sólo para mostrar a que satélite cooresponde cada correccion

PD = {                # PseudoDist [m]
"28": 23334619.807,
"13": 22586189.129,
"01": 25167667.280,
"27": 20873169.266,
"24": 23371141.291,
"10": 21505922.486,
"08": 20958408.428
}
# Obtengo tiempo de tránsito [s]
TT = {s:PD[s]/c for s in PD}

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

Delta_ant = {          # Corrección de centro de fase de antena
"28":  1.04280 ,       #     para cada sat [m]
"13":  1.38950 ,
"01":  2.38080 ,
"27":  2.63340 ,
"24":  2.60380 ,
"10":  2.54650 ,
"08":  2.57810 
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

#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
#  Levanto las efemérides transmitidas . . .
# ------------------------------------------------------------------------------------
filename = "Practica_2\\ifag0780.01n"
start = False
mensajes = {}
satlist = []

# para hacerla fácil, anoto el tiempo correspondiente a las
#  efemérides transmitidas que voy a usar para cada satélite ...
tefem = {} # datetime(2001,3,19,0,0,0)


linea = 0
sat=""
sata=sat
with open(filename,encoding="utf-8") as f:
    for line in f:
        if start:
            if line[5]!='.':
                if sat != "":
                    if sat not in satlist:
                        satlist.append(sat)
                        #print("Nuevo satelite: ",sat)
                        mensajes[sata] = []
                        tefem[sata] = fechaHora
                    mensajes[sata].append({"FechaHora": fechaHora,"a0": a0 ,"a1": a1,"a2": a2,
                                    "T0e": Toe,"GPSweek": GPS_Week,"sqrtA" : sqrtA,"e": e, "M0": M0,
                                    "omega":omega,"i0": i0,"OMEGA": OMEGA,"Delta_n": Delta_n,
                                    "idot": idot,"OMEGA_DOT": OMEGA_DOT,"Cuc": Cuc,"Cus": Cus,
                                    "Crc": Crc,"Crs": Crs,"Cic": Cic,"Cis": Cis})
                    GPSweek = GPS_Week
                sat = line[:2]
                if sat[0]==' ':
                    sata='0'+sat[1]
                else:
                    sata=sat
                #print(sat,sata)
                seg = line[18:22].split('.')
                if seg[0]=='  ':
                    seg[0]='0'
                AA = line[3:5]
                if int(AA)<10:
                    AA = "0" + line[4]
                    #print(AA)
                    #print(seg)
                #                                 AA                   MM              DD
                fechaHora =  datetime(int("20"+AA),int(line[6:8]),int(line[9:11]),int(line[12:14])
                                      #      HH           mm             ss
                                      ,int(line[15:17]),int(seg[0]),int(seg[1]))
                a0 = float(line[22:41].replace('D','E'))
                a1 = float(line[41:60].replace('D','E'))
                a2 = float(line[60:79].replace('D','E'))
                linea = 0
            elif line.startswith("   "):
                linea += 1
                if linea==1:
                    IODE   = float(line[ 3:22].replace('D','E'))
                    Crs    = float(line[22:41].replace('D','E'))
                    Delta_n= float(line[41:60].replace('D','E'))
                    M0     = float(line[60:79].replace('D','E'))
                elif linea==2:
                    Cuc   = float(line[ 3:22].replace('D','E'))
                    e     = float(line[22:41].replace('D','E'))
                    Cus   = float(line[41:60].replace('D','E'))
                    sqrtA = float(line[60:79].replace('D','E'))
                elif linea==3:
                    Toe   = float(line[ 3:22].replace('D','E'))
                    Cic   = float(line[22:41].replace('D','E'))
                    OMEGA = float(line[41:60].replace('D','E'))
                    Cis   = float(line[60:79].replace('D','E'))
                elif linea==4:
                    i0       = float(line[ 3:22].replace('D','E'))
                    Crc      = float(line[22:41].replace('D','E'))
                    omega    = float(line[41:60].replace('D','E'))
                    OMEGA_DOT= float(line[60:79].replace('D','E'))
                elif linea==5:
                    idot    = float(line[ 3:22].replace('D','E'))
                    c2      = float(line[22:41].replace('D','E'))
                    GPS_Week= float(line[41:60].replace('D','E'))
                    L2_P    = float(line[60:79].replace('D','E'))
                elif linea==6:
                    SVa  = float(line[ 3:22].replace('D','E'))
                    SVh  = float(line[22:41].replace('D','E'))
                    TGD  = float(line[41:60].replace('D','E'))
                    IODC = float(line[60:79].replace('D','E'))
                elif linea==7:
                    TxToM  = float(line[ 3:22].replace('D','E'))
                    FitInt = float(line[22:41].replace('D','E'))

        if "LEAP SECONDS" in line:
            L_sec = int(line[:8])
            #print(L_sec)
        if "END OF HEADER" in line:
            start = True

tGPS0 += timedelta(days=GPSweek*7)
# el momento en que se realiza la observación
tobs = datetime(2001,3,19,0,15,0)
tobs_seg = int((tobs-tGPS0).total_seconds())


# ------------------------------------------------------------------------------------
#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -


def newton(fn,Df,x0,epsilon,max_iter):
    '''Approximate solution of f(x)=0 by Newton's method.

    Parameters
    ----------
    fn : function
        Function for which we are searching for a solution f(x)=0.
    Df : function
        Derivative of f(x).
    x0 : number
        Initial guess for a solution f(x)=0.
    epsilon : number
        Stopping criteria is abs(f(x)) < epsilon.
    max_iter : integer
        Maximum number of iterations of Newton's method.

    Returns
    -------
    xn : number
        Implement Newton's method: compute the linear approximation
        of f(x) at xn and find x intercept by the formula
            x = xn - f(xn)/Df(xn)
        Continue until abs(f(xn)) < epsilon and return xn.
        If Df(xn) == 0, return None. If the number of iterations
        exceeds max_iter, then return None.

    Examples
    --------
    >>> f = lambda x: x**2 - x - 1
    >>> Df = lambda x: 2*x - 1
    >>> newton(f,Df,1,1e-8,10)
    Found solution after 5 iterations.
    1.618033988749989
    '''
    xn = x0
    for _ in range(0,max_iter):
        fxn = fn(xn)
        if abs(fxn) < epsilon:
            #print('Found solution after',n,'iterations.')
            return xn
        Dfxn = Df(xn)
        if Dfxn == 0:
            print('Zero derivative. No solution found.')
            return None
        xn = xn - fxn/Dfxn
    print('Exceeded maximum iterations. No solution found.')
    return None

def calPos(h,inst,dt):
    a = h["sqrtA"] * h["sqrtA"]
    a_cubo = a*a*a

    M = h["M0"] + (sqrt(μ/a_cubo)+h["Delta_n"])*dt

    fE  = lambda x: M + h["e"]*sin(x) - x
    dfE = lambda x: h["e"]*cos(x)-1
    E = newton(fE,dfE,0,1E-8,4)
    r0 = a * (1 - h["e"] * cos(E))
    f = 2*atan(sqrt(1+h["e"])/sqrt(1-h["e"])*tan(E/2))
    u0 = h["omega"] + f

    Ω = h["OMEGA"] + h["OMEGA_DOT"] * dt
    ω = h["omega"] + h["Cuc"] * cos(2*u0) + h["Cus"] * sin(2*u0)
    r = r0 + h["Crc"] * cos(2*u0) + h["Crs"] * sin(2*u0)
    i = h["i0"] + h["Cic"] * cos(2*u0) + h["Cis"] * sin(2*u0) + h["idot"] * dt
    ϴ = ωe * inst
    u = ω + f
    """
    print()
    print(HoraCalc.strftime("Hora: %d/%m/%Y %H:%M:%S.%f"))
    print("dt:",dt)
    print("M:",M)
    print("E:",E)
    print("f:",f)
    if platform.system() == "Linux":
        respuesta = readchar.readchar().decode('utf-8')
    elif platform.system() == "Windows":
        respuesta = msvcrt.getch().decode('utf-8')

    if respuesta =="x":
        exit()    

    """
    r += Delta_ant
    xyz = [ [r*cos(u), 0, 0],
            [r*sin(u), 0, 0],
            [    0   , 0, 0]]

    R31 = [ [  cos(-Ω),sin(-Ω), 0],
            [ -sin(-Ω),cos(-Ω), 0],
            [      0 ,      0 , 1]]

    R32 = [ [  cos(ϴ), sin(ϴ),  0],
            [ -sin(ϴ), cos(ϴ),  0],
            [      0 ,     0 ,  1]]

    R1 =  [ [ 1,      0 ,       0],
            [ 0, cos(-i), sin(-i)],
            [ 0,-sin(-i), cos(-i)]]

    R3R3 = matmul(R32,R31)
    RR   = matmul(R3R3,R1)
    xyz_prima  = matmul(RR,xyz)

    return [xyz_prima[0][0],xyz_prima[1][0],xyz_prima[2][0],r,ϴ-Ω,-i]

def coord_obs():
    """ Arma la matriz con coordenadas de cada satélite
         en el tiempo de observación
    """
    Calc = {}
    for sati in TT:
        s= mensajes[sati]
        for h in s:
            if h["FechaHora"]==tefem[sati]:
                tefem_seg = int((tefem[sati]-tGPS0).total_seconds())
                #HoraCalc = tefem + timedelta(seconds=tobs_seg-TT[sati])
                a = h["sqrtA"] * h["sqrtA"]
                a_cubo = a*a*a
                # tiempo entre la efemeride transmitida y
                # el momento en que quiero saber las coordenadas del satélite
                Delta_t = (tobs_seg - TT[sati]) - tefem_seg

                M = h["M0"] + (sqrt(μ/a_cubo)+h["Delta_n"])*Delta_t

                fE  = lambda x: M + h["e"]*sin(x) - x
                dfE = lambda x: h["e"]*cos(x)-1
                E = newton(fE,dfE,0,1E-6,4)
                r0 = a * (1 - h["e"] * cos(E))
                f = 2*atan(sqrt(1+h["e"])/sqrt(1-h["e"])*tan(E/2))
                u0 = h["omega"] + f

                Ω = h["OMEGA"] + h["OMEGA_DOT"] * Delta_t
                ω = h["omega"] + h["Cuc"] * cos(2*u0) + h["Cus"] * sin(2*u0)
                r = r0 + h["Crc"] * cos(2*u0) + h["Crs"] * sin(2*u0)
                i = h["i0"] + h["Cic"] * cos(2*u0) + h["Cis"] * sin(2*u0) + h["idot"] * Delta_t
                ϴ = ωe * (tobs_seg - TT[sati])
                u = ω + f
                r += Delta_ant[sati]
                xyz = [ [r*cos(u), 0, 0],
                        [r*sin(u), 0, 0],
                        [    0   , 0, 0]]

                R31 = [ [  cos(-Ω),sin(-Ω), 0],
                        [ -sin(-Ω),cos(-Ω), 0],
                        [      0 ,      0 , 1]]

                R32 = [ [  cos(ϴ), sin(ϴ),  0],
                        [ -sin(ϴ), cos(ϴ),  0],
                        [      0 ,     0 ,  1]]

                R1 =  [ [ 1,      0 ,       0],
                        [ 0, cos(-i), sin(-i)],
                        [ 0,-sin(-i), cos(-i)]]

                R3R3 = matmul(R32,R31)
                RR   = matmul(R3R3,R1)

                xyz_prima  = matmul(RR,xyz)
                x = float(xyz_prima[0][0])
                y = float(xyz_prima[1][0])
                z = float(xyz_prima[2][0])
                Calc[sati] = [x, y,z,Precisas[sati][3]]
    return Calc

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
    global satord
    L = []
    A = []
    C = []
    satord=[]
    j = 0
    for st in Calculadas:
        s= Calculadas[st]
        """
        dX = s[0] * 1000 - Coord[0]
        dY = s[1] * 1000 - Coord[1]
        dZ = s[2] * 1000 - Coord[2]
        """
        dX = s[0] - Coord[0]
        dY = s[1] - Coord[1]
        dZ = s[2] - Coord[2]
        ρ = sqrt(dX*dX+dY*dY+dZ*dZ)
        fila = [dX/ρ,dY/ρ,dZ/ρ,1]  # opcion con incognita c * Delta_t

        A.append(fila)
        # diferencia Observado - Calculado
        err = PD[st] - ρ + float(Coord[3]) + c * s[3] / 1E6  # Si s[3] > 0 el satélite atrasa con respecto a GPS time, entonces "sumo" error en distancia?
        L.append(err)
        satord.append(st)
        linea_C =[0 for i in range(cant_sat)]
        #linea_C[j]= ρ - float(Coord[3]) + c * s[3] / 1E6 - sqrt(Coord[0]*Coord[0]+Coord[1]*Coord[1]+Coord[2]*Coord[2])
        linea_C[j]= err * err
        j +=1
        C.append(linea_C)

def imprime_resu():
    """Imprime resultados parciales"""

    print("\nObservado-calculado")
    print("SAT         Dif\n")
    for j in range(len(L)):
            print(satord[j],"  {:15.6f}  ".format(L[j]))
    print()

    """
    for dif in L:
        print("{:15.6f}  ".format(dif))
    print()

    for linea in P:
        print(linea)
    print()
    """
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
    """
    #print("\n")
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
    print("           X                   Y                   Z                   t")
    linea =""
    for i in X1:
        linea += "{:20.16f}  ".format(i)
    print(linea)
    dift = list(X1[:3])
    print("\nMódulo de la corrección calculada: {:20.16f} m".format(linalg.norm(dift)))
    dift.append(X1[3]*c)
    print("Módulo de la corrección calculada\n              incluyendo el reloj: {:10.6f} m\n".format(linalg.norm(dift)))
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
"""

os.system('cls')
print()
for s in sorted(mensajes):
    for item in mensajes[s]:
        print(s, item["FechaHora"])
    print('------------------------------------------------------\n')
"""
if platform.system() == "Linux":
    # import readchar # type: ignore
    os.system('clear')
elif platform.system() == "Windows":
    # import msvcrt   # type: ignore
    os.system('cls')
"""
print(satlist)
print()
"""
# consigo las coordenadas de cada satélite
#         en el tiempo de observación
Calculadas = coord_obs()
"""
for s in Calculadas:
    print("Calc: ",Calculadas[s])
    print("Prec: ",Precisas[s])
print()

for s in TT:
    print(s,tefem[s])

"""
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
