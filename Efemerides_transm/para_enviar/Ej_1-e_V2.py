# pylint: disable=W0640,W0611,E0601,C0209,C0301

import os
import platform
from datetime import datetime, timedelta
from math import sin, cos, tan, atan, sqrt, pi
from numpy import matmul,dot
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
myFmt = mdates.DateFormatter('%H:%M')

μ = 3.986004418E14 # m3/s2  Earth gravitational constant
ωe =   7.2921151467E-5 # radians/s Angular Velocity of the Earth
tGPS0 =  datetime(1980,1,6,0,0,0)
sati = " 1"   # el satélite que me interesa lo defino de entrada para filtrar en la lectura de los archivos

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
Aqui se lee el archivo de efemérides precisas sólo para el satélite que nos interesa
Arma una lista de listas orbitas precisas (precorbitas) cada 15' como vienen en
el archivo los valores de la órbita en x, y , z

"""
filename = "igs11061.sp3"
precorbitas = []
with open(filename) as f:
    for line in f:
        if line.startswith("*  2001"):
            seg = line[19:].split('.')
            fechaHora =  datetime(int(line[1:7]),int(line[7:10]),int(line[10:13]),int(line[13:16]),int(line[16:19]),int(seg[0]),int(seg[1]))
        elif line.startswith("P"):
            sat = line[2:4]
            if sat==sati:
                x = float(line[5:19])*1000
                y = float(line[19:33])*1000
                z = float(line[33:47])*1000
                fila = [fechaHora,x,y,z]
                precorbitas.append(fila)

#print(precorbitas)

"""
Aqui se lee el archivo de efemérides transmitidas sólo para el satélite que nos interesa
Arma una lista de mensajes de navegación (mensajes) que tiene para cada tiempo reportado
un diccionario con entradas correspondientes a cada uno de los parámetros orbitales transmitidos

"""
filename = "Efemerides_transm\\ifag0780.01n"
start = False
mensajes = []
linea = 0
sat=""
with open(filename,encoding="utf-8") as f:
    for line in f:
        if start:
            if line[5]!='.':
                if sat != "":
                    if sat==sati:
                        mensajes.append({"FechaHora": fechaHora,"a0": a0 ,"a1": a1,"a2": a2, "T0e": Toe,
                                         "GPSweek": GPS_Week,"sqrtA" : sqrtA,"e": e, "M0": M0,
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

#print(mensajes)

tGPS0 += timedelta(days=GPSweek*7)

tref = datetime(2001,3,19,2,0,0)
tref_seg = int((tref-tGPS0).total_seconds())

t_inic = datetime(2001,3,19,0,0,0)
ti_seg = int((t_inic-tGPS0).total_seconds())
t_fin = datetime(2001,3,19,8,15,0)
tf_seg = int((t_fin-tGPS0).total_seconds())

# Estas son listas que armo para después plotearlo si quiero . . .
# / - - - \/ - - - - \
times = []
dr = []
dx = []
dy = []
dz = []
# \ - - - /\ - - - - /

# Busco en el archivo de mensajes
#  los mensajes correspondientes al tref que me piden
#    --> 19/03/2001 a las 2:00
for h in mensajes:
    if h["FechaHora"]==tref:
        #
        #  Ahora h es un diccionario que tiene solamente los mensajes de navegación
        #    del satelite 1 a tref cada uno con su nombre . . .
        #
        for t in range(ti_seg,tf_seg,15*60):
            # Voy haciendo cada 15' las transformaciones
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

            xyz_prima = [[r*cos(u), 0, 0],
                         [r*sin(u), 0, 0],
                         [    0   , 0, 0]]

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
            xyz  = matmul(RR,xyz_prima)
            x = xyz[0][0]
            y = xyz[1][0]
            z = xyz[2][0]

            # Ya tengo las coordendas x, y, z transmitidas para el tiempo "t"
            #   y las tengo que comparar con las precisas en ese mismo tiempo
            #   para lo cual busco en la lista de órbitas precisas
            #   los datos del instante calculado arriba (HoraCalc)

            for j in precorbitas:
                if j[0]==HoraCalc:
                    dxx = j[1] - x
                    dyy = j[2] - y
                    dzz = j[3] - z
                    rp = sqrt(j[1]*j[1]+j[2]*j[2]+j[3]*j[3])
                    drr = rp - r

                    dx.append(dxx)
                    dy.append(dyy)
                    dz.append(dzz)
                    dr.append(drr)
                    times.append(HoraCalc)

                    print(j[0].strftime("Hora: %d/%m/%Y %H:%M"))
                    print("r calculada: {:14.3f}      r precisa: {:14.3f}       r diff: {:14.3f}".format(r,rp,drr))
                    print("X calculada: {:14.3f}      X precisa: {:14.3f}       X diff: {:14.3f}".format(x,j[1],dxx))
                    print("Y calculada: {:14.3f}      Y precisa: {:14.3f}       Y diff: {:14.3f}".format(y,j[2],dyy))
                    print("Z calculada: {:14.3f}      Z precisa: {:14.3f}       Z diff: {:14.3f}".format(z,j[3],dzz))


print()
print("pulse cualquier tecla para plotear, 'x' o 'c' para cancelar!")    
resp_set = ['x','c','X','C']
if platform.system() == "Linux":
    respuesta = readchar.readchar().decode('utf-8')
elif platform.system() == "Windows":
    respuesta = msvcrt.getch().decode('utf-8')

print()
if not respuesta in resp_set:
    # ---------------------------------------------------------------------------
    #     P L O T E A N D O  ! !
    # ---------------------------------------------------------------------------
    fig, ax = plt.subplots(4, sharex=True, sharey=False, gridspec_kw={'hspace': 0})
    fig.set_size_inches(12, 7)   # w , h
    ax[0].plot(times, dr,'o', c='darkslategray')
    ax[0].set(ylabel= "r [m]")
    ax[1].plot(times, dx,'o', c='darkslategray')
    ax[1].set(ylabel= "x [m]")
    ax[2].plot(times, dy,'o', c='darkslategray')
    ax[2].set(ylabel= "y [m]")
    ax[3].plot(times, dz,'o', c='darkslategray')
    ax[3].set(ylabel= "z [m]")

    for axs in ax.flat:
        axs.label_outer()

    fig.suptitle("Diferencias en órbitas transmitidas vs calculadas el 19/03/2001 de 0:00 a 8:00",fontsize=13)
    plt.xlabel("Hora del día (UTC)",labelpad=4 ,fontsize=13)

    plt.gca().xaxis.set_major_formatter(myFmt)
    plt.xticks(rotation=45, ha='right',fontsize=10)
    plt.yticks(fontsize=6)
    plt.show()
