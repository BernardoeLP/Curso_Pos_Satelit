""" Practica 2-i
    para el punto i) ..... ??

# pylint: disable= C0103, C0200, C0206, C0301, C0209, E0601, W0105, W0602, W0603, W0621, W0640
"""
# pylint: disable= C0103, C0200, C0206, C0301, C0209, E0601, W0105, W0602, W0603, W0621, W0640

import os
import platform
from datetime import datetime, timedelta
from math import sin, cos, tan, atan, sqrt #, acos, pi
from numpy import  matmul #, transpose, linalg, dot, std, subtract, cross


c = 299792458          # m/s  de ITRF
μ = 3.986005E14        # m3/s2  Earth gravitational constant
ωe =   7.2921151467E-5 # radians/s Angular Velocity of the Earth
ωE = [ 0, 0, ωe]       # Velocity of the Earth (vector)
tGPS0 =  datetime(1980,1,6,0,0,0)


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

    #r += Delta_ant
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


############################### -  -  -  -  -  -  -  -  -  -  -  -  -
#   I N I C I O    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
############################### -  -  -  -  -  -  -  -  -  -  -  -  -
if platform.system() == "Linux":
    # import readchar # type: ignore
    os.system('clear')
elif platform.system() == "Windows":
    # import msvcrt   # type: ignore
    os.system('cls')


#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
#  Levanto las efemérides transmitidas . . .
# ------------------------------------------------------------------------------------
filename = "Practica_3\\vill1440.02n"
start = False
mensajes = {}
satlist = []

# para hacerla fácil, anoto el tiempo correspondiente a las
#  efemérides transmitidas que voy a usar para cada satélite ...
tefem = {}
linea = 0
sat=""
sata=sat
with open(filename,encoding="utf-8") as fnav:
    for line in fnav:
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
# ------------------------------------------------------------------------------------
#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -



for s in satlist:
    print("Satélite: ",s)


fs = open("Practica_3\\sat7_nav.csv",'w')
fs.write("FechaHora,X,Y,Z\n")
"""
for entrada in mensajes["G07"]:
    posicion = calPos(entrada,,)
    fila = entrada["FechaHora"].strftime("%d/%m/%Y %H:%M:%S.%f")+",{:15.5f},{:15.5f},{:15.5f},{:15.5f},{:15.5f},{:15.5f}\n".format(
        entrada["C1"], entrada["P2"], entrada["L1"], entrada["L2"], entrada["D1"], entrada["D2"])
    fs.write(fila)
"""
fs.close()
