# pylint: disable=W0640,W0611,E0601,C0209,C0301

import os
import platform
from datetime import datetime, timedelta
from math import sin, cos, tan, atan, sqrt, pi
from numpy import matmul,dot
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
myFmt = mdates.DateFormatter('%H:%M')

μ = 3.986005E14 # m3/s2  Earth gravitational constant
ωe =   7.2921151467E-5 # radians/s Angular Velocity of the Earth
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


if platform.system() == "Linux":
    import readchar
    os.system('clear')
elif platform.system() == "Windows":
    import msvcrt
    os.system('cls')

"""
Aqui se lee el archivo de efemérides precisas para todos los satélites
(podría hacerse sólo para el que nos interesa también)
Arma un diccionario de orbitas precisas (precorbitas) que tiene una entrada para cada satélite
Para cada satélite (x ej " 1" o "22") tiene cada 15' como viene en el archivo los valores de
 la órbita en x, y , z

"""
filename = "igs11061.sp3"
precsatlist = []
precorbitas = {}
with open(filename,encoding="utf-8") as f:
    for line in f:
        if line.startswith("*  2001"):
            seg = line[19:].split('.')
            fechaHora =  datetime(int(line[1:7]),int(line[7:10]),int(line[10:13]),int(line[13:16]),int(line[16:19]),int(seg[0]),int(seg[1]))
        elif line.startswith("P"):
            sat = line[2:4]
            x = float(line[5:19])*1000
            y = float(line[19:33])*1000
            z = float(line[33:47])*1000
            fila = [fechaHora,x,y,z]
            if not sat in precsatlist:
                precsatlist.append(sat)
                precorbitas[sat]=[]
            precorbitas[sat].append(fila)

#print(precsatlist)

"""
Aqui se lee el archivo de efemérides transmitidas para todos los satélites
(podría hacerse sólo para el que nos interesa también)
Arma un diccionario de mensajes de navegación (mensajes) que tiene una entrada para cada satélite
Para cada satélite (x ej " 1" o "22") tiene para cada tiempo reportado un diccionario con entradas 
correspondientes a cada uno de los parámetros orbitales transmitidos

"""
filename = "ifag0780.01n"
start = False
mensajes = {}
satlist = []
linea = 0
sat=""
with open(filename,encoding="utf-8") as f:
    for line in f:
        if start:
            if line[5]!='.':
                if sat != "":
                    if sat not in satlist:
                        satlist.append(sat)
                        mensajes[sat] = []
                    mensajes[sat].append({"FechaHora": fechaHora,"a0": a0 ,"a1": a1,"a2": a2,
                                    "T0e": Toe,"GPSweek": GPS_Week,"sqrtA" : sqrtA,"e": e, "M0": M0,
                                    "omega":omega,"i0": i0,"OMEGA": OMEGA,"Delta_n": Delta_n,
                                    "idot": idot,"OMEGA_DOT": OMEGA_DOT,"Cuc": Cuc,"Cus": Cus,
                                    "Crc": Crc,"Crs": Crs,"Cic": Cic,"Cis": Cis})
                    GPSweek = GPS_Week
                sat = line[:2]
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

tref = datetime(2001,3,19,2,0,0)
tref_seg = int((tref-tGPS0).total_seconds())

t_inic = datetime(2001,3,19,0,0,0)
ti_seg = int((t_inic-tGPS0).total_seconds())

t_fin = datetime(2001,3,19,8,15,0)
tf_seg = int((t_fin-tGPS0).total_seconds())

#print(mensajes)

sati = " 1"   # el satélite que me interesa
s= mensajes[sati]   # me quedo sólo con los mensajes del satélite que me interesa


# Estas son listas que armo para después plotearlo si quiero . . . 
# / - - - \/ - - - - \
posicion_x = []
posicion_y = []
dx = []
dy = []
dz = []
times = []
dr = []
# \ - - - /\ - - - - /

# Busco en el archivo de mensajes del satélite de interés
#  los mensajes del tref que me piden --> 19/03/2001 a las 2:00
for h in s:
    if h["FechaHora"]==tref:
        #
        #  Ahora h es un diccionario que tiene solamente los mensajes de navegación 
        #    del satelite 1 a tref cada uno con su nombre . . . 
        #

        #print("Satélite: "+ sati)
        #print(h["FechaHora"])
        #print(h)
        for t in range(ti_seg,tf_seg,15*60):
            # Voy haciendo cada 15' lo las transformaciones
            #   necesarias para calcular las coordendas del satélite
            #   en ese momento (HoraCalc)
            HoraCalc = t_inic + timedelta(seconds=t-ti_seg)
            a = h["sqrtA"] * h["sqrtA"]
            a_cubo = a*a*a
            Delta_t = t - tref_seg

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
            ϴ = ωe * t
            u = ω + f
            """
            print()
            print(HoraCalc)
            print("a : ",a)
            print("M : ",M)
            print("Ω : ",Ω)
            print("E : ",E)
            print("r0: ",r0)
            print("f : ",f)
            print("u0: ",u0)
            print("ω : ",ω)
            print("r : ",r)
            print("i : ",i)
            print("ϴ : ",ϴ)
            print("u : ",u)
            
            print()
            print(HoraCalc)
            print("r: {:14.3f}      u[º]: {:8.3f} ".format(r, u*180/pi))
            """
            xyz_prima = [[r*cos(u), 0, 0],
                         [r*sin(u), 0, 0],
                         [    0   , 0, 0]]

            posicion_x.append(r*cos(u))
            posicion_y.append(r*sin(u))

            R31 = [  [  cos(-Ω),sin(-Ω), 0],
                     [ -sin(-Ω),cos(-Ω), 0],
                     [      0 ,      0 , 1]]

            R32 = [  [  cos(ϴ), sin(ϴ),  0],
                     [ -sin(ϴ), cos(ϴ),  0],
                     [      0 ,     0 ,  1]]

            R1 =  [  [ 1,      0 ,       0],
                     [ 0, cos(-i), sin(-i)],
                     [ 0,-sin(-i), cos(-i)]]
            R3R3 = matmul(R32,R31)
            RR   = matmul(R3R3,R1)
            """
            print()
            print(HoraCalc)
            print(R3R3)
            print()
            print(RR)
            """

            xyz  = matmul(RR,xyz_prima)
            x = xyz[0][0]
            y = xyz[1][0]
            z = xyz[2][0]

            # Ya tengo las coordendas x, y, z transmitidas para el tiempo "t"
            # y las tengo que comparar con las precisas en ese mismo tiempo
            # para lo cual busco en la entrada del diccionario correspondiente
            # al satélite en cuestión, los datos del instante calculado arriba (HoraCalc)

            for j in precorbitas[sati]:
                if j[0]==HoraCalc:
                    dxx = j[1] - x
                    dyy = j[2] - y
                    dzz = j[3] - z
                    dx.append(dxx)
                    dy.append(dyy)
                    dz.append(dzz)
                    rp = sqrt(j[1]*j[1]+j[2]*j[2]+j[3]*j[3])
                    drr = r-rp
                    times.append(HoraCalc)
                    dr.append(drr)

                    print(j[0].strftime("Hora: %d/%m/%Y %H:%M"))#+"   Satélite: " + sati)
                    print("r calculada: {:14.3f}      r precisa: {:14.3f}       r diff: {:14.3f}".format(r,rp,drr))
                    print("X calculada: {:14.3f}      X precisa: {:14.3f}       X diff: {:14.3f}".format(x,j[1],dxx))
                    print("Y calculada: {:14.3f}      Y precisa: {:14.3f}       Y diff: {:14.3f}".format(y,j[2],dyy))
                    print("Z calculada: {:14.3f}      Z precisa: {:14.3f}       Z diff: {:14.3f}".format(z,j[3],dzz))


plt.plot(times,dr)
plt.gca().xaxis.set_major_formatter(myFmt)
plt.xticks(rotation=90)
plt.show()
