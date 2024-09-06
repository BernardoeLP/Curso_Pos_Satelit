""" Practica 3
Lee archivo RINEX de Navegación:    
     2.10           N: GPS NAV DATA                         RINEX VERSION / TYPE

Y graba las coordenadas ECEF del satélite "sati" en un archivo ".csv"

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
sati = "07"            # Satelite del que quiero obtener las coordenadas

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
    # Calculate the semi-major axis squared
    a = h["sqrtA"] * h["sqrtA"]
    a_cubo = a * a * a  # Cube of the semi-major axis

    # Mean anomaly calculation
    M = h["M0"] + (sqrt(μ / a_cubo) + h["Delta_n"]) * dt

    # Define the function and its derivative for solving Kepler's equation
    fE = lambda x: M + h["e"] * sin(x) - x
    dfE = lambda x: h["e"] * cos(x) - 1
    E = newton(fE, dfE, 0, 1E-8, 4)  # Solve for eccentric anomaly

    # Calculate the radius and true anomaly
    r0 = a * (1 - h["e"] * cos(E))
    f = 2 * atan(sqrt(1 + h["e"]) / sqrt(1 - h["e"]) * tan(E / 2))
    u0 = h["omega"] + f  # Argument of latitude

    # Calculate the corrected longitude of ascending node
    Ω = h["OMEGA"] + h["OMEGA_DOT"] * dt
    # Calculate the corrected argument of perigee
    ω = h["omega"] + h["Cuc"] * cos(2 * u0) + h["Cus"] * sin(2 * u0)
    # Calculate the corrected radius
    r = r0 + h["Crc"] * cos(2 * u0) + h["Crs"] * sin(2 * u0)
    # Calculate the corrected inclination
    i = h["i0"] + h["Cic"] * cos(2 * u0) + h["Cis"] * sin(2 * u0) + h["idot"] * dt
    ϴ = ωe * inst  # Earth rotation angle
    u = ω + f  # Argument of latitude

    #r += Delta_ant

    # Position in orbital plane
    xyz = [ [r*cos(u), 0, 0],
            [r*sin(u), 0, 0],
            [    0   , 0, 0]]

    # Rotation matrices
    R31 = [ [  cos(-Ω),sin(-Ω), 0],
            [ -sin(-Ω),cos(-Ω), 0],
            [      0 ,      0 , 1]]

    R32 = [ [  cos(ϴ), sin(ϴ),  0],
            [ -sin(ϴ), cos(ϴ),  0],
            [      0 ,     0 ,  1]]

    R1 =  [ [ 1,      0 ,       0],
            [ 0, cos(-i), sin(-i)],
            [ 0,-sin(-i), cos(-i)]]

    # Combine rotation matrices
    R3R3 = matmul(R32,R31)
    RR   = matmul(R3R3,R1)
    xyz_prima  = matmul(RR,xyz)

    # Return the final position in ECEF coordinates
    return [xyz_prima[0][0],xyz_prima[1][0],xyz_prima[2][0]]


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
mensajes = []  # List to store messages for each satellite
fechas = []    # List to store dates
linea = 0      # Line counter for satellite data
cont = 0       # Counter for the number of messages
sat = ""       # Current satellite identifier
with open(filename,encoding="utf-8") as f:
    for line in f:
        if start:
            if line[5]!='.':
                if sat != "":
                    if sat==sati:
                        # Append the date and message to their respective lists
                        fechas.append(fechaHora)
                        mensajes.append({"FechaHora": fechaHora,"a0": a0 ,"a1": a1,"a2": a2, "T0e": Toe,
                                         "GPSweek": GPS_Week,"sqrtA" : sqrtA,"e": e, "M0": M0,
                                         "omega":omega,"i0": i0,"OMEGA": OMEGA,"Delta_n": Delta_n,
                                         "idot": idot,"OMEGA_DOT": OMEGA_DOT,"Cuc": Cuc,"Cus": Cus,
                                         "Crc": Crc,"Crs": Crs,"Cic": Cic,"Cis": Cis})
                        cont += 1
                        GPSweek = GPS_Week
                # Parse the satellite identifier and timestamp
                sat = line[:2].replace(' ','0')
                useg = str((float(line[18:22])-int(line[18:20]))* 1E6)
                seg = []
                seg.append(line[18:20].replace(' ','0'))
                seg.append(useg.split('.',maxsplit=1)[0])

                AA = line[3:5]
                AA.replace(' ','0')
                # Create a datetime object from the parsed date and time components
                #                            AA             MM             DD             HH
                fechaHora =  datetime(int("20"+AA),int(line[6:8]),int(line[9:11]),int(line[12:14])
                                      #        mm             ss           useg
                                      ,int(line[15:17]),int(seg[0]),int(seg[1]))
                # Parse the polynomial coefficients
                a0 = float(line[22:41].replace('D','E'))
                a1 = float(line[41:60].replace('D','E'))
                a2 = float(line[60:79].replace('D','E'))
                linea = 0
            elif line.startswith("   "):
                linea += 1
                # Parse the remaining satellite data based on the line number
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
                    #FitInt = float(line[22:41].replace('D','E'))

        if "LEAP SECONDS" in line:
            L_sec = int(line[:8])  # Parse the leap seconds
            #print(L_sec)
        if "END OF HEADER" in line:
            start = True           # Start processing data after the header

t_ref = tGPS0+timedelta(days=GPSweek*7) # Tiempo en el que inicia la semana GPS actual
                                        # y es la referencia para los tiempos de las efemérides


# ------------------------------------------------------------------------------------
#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -

# Open the output file for writing satellite coordinates
fs = open("Practica_3\\sat7_coord.csv",'w')
fs.write("FechaHora,X,Y,Z\n")
 
with open("Practica_3\\sat7.csv") as fobs:
    l = 0  # Line counter
    for l_obs in fobs:
        if l>0:
            # Parse the observation epoch
            pos_epoch = l_obs.split(',')[0]
            # tiempo al que quiero calcular las coordendas
            epoch = datetime.strptime(l_obs.split(',')[0],"%d/%m/%Y %H:%M:%S.%f")
            epoch_seg = int((epoch-t_ref).total_seconds())
            # tiempo de la efeméride más cercana
            # Find the closest ephemeris time
            t_efem = min(fechas, key=lambda x: abs(x - epoch))
            t_efem_seg = int((t_efem-t_ref).total_seconds())
            Delta_t = epoch_seg - t_efem_seg
            print(epoch,"    -->  ", t_efem)
            # Retrieve the ephemeris data for the closest time
            efemerides = next(item for item in mensajes if item["FechaHora"] == t_efem)
            # Calculate the satellite position
            posicion = calPos(efemerides,epoch_seg,Delta_t)
            print(posicion)
            print()
            # Write the satellite position to the output file
            #  Grabo en el archivo de posiciones del satélite:  "FechaHora, X ,Y ,Z"
            fs.write(epoch.strftime("%d/%m/%Y %H:%M:%S.%f")+",{},{},{}\n".format(float(posicion[0]),float(posicion[1]),float(posicion[2])))
        l += 1

fs.close()  # Close the output file
