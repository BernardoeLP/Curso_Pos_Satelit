import os
from datetime import datetime
from math import sin, cos, asin, acos, atan2, sqrt
import platform
import plotly.graph_objs as go

if platform.system() == "Linux":
    import readchar
    os.system('clear')  
elif platform.system() == "Windows":
    import msvcrt
    os.system('cls')


a = 6378137.0
f = 1 / 298.257223563
e2 = f*(2-f)
e= sqrt(e2)

estaciones = {
'LPGS' : [float("2780089.996"),float("-4437398.174"),float("-3629387.406")],
'TREL' : [float("1938169.196"),float("-4229010.224"),float("-4348880.32")],
'RIO2' : [float("1429907.806"),float("-3495354.833"),float("-5122698.637")]
}


def ecef_to_geodetic(x, y, z):
    '''
    Olson, D. K. (1996). Converting Earth-Centered, Earth-Fixed Coordinates to Geodetic Coordinates. 
    IEEE Transactions on Aerospace and Electronic Systems. https://doi.org/10.1109/7.481290

    Converted to Python and modified by Steven Ward.  No rights reserved.
    '''
    w2 = x * x + y * y
    w = sqrt(w2)
    z2 = z * z
    lon_rad = atan2(y, x)

    a1 = a * e2
    a2 = a1 * a1
    a3 = a1 * e2 / 2
    a4 = 2.5 * a2
    a5 = a1 + a3
    #a6 = (1 - self.e2)

    r2 = w2 + z2
    r = sqrt(r2)

    s2 = z2 / r2
    c2 = w2 / r2
    u = a2 / r
    v = a3 - a4 / r

    s = 0
    c = 0
    ss = 0

    # cos(45°)² == ½
    if c2 > 0.5: # Equatorial
        s = (z / r) * (1 + c2 * (a1 + u + s2 * v) / r)
        lat_rad = asin(s)
        ss = s * s
        c = sqrt(1 - ss)
    else: # Polar
        c = (w / r) * (1 - s2 * (a5 - u - c2 * v) / r)
        lat_rad = acos(c)
        ss = 1 - c * c
        s = sqrt(ss)

        if z < 0:
            lat_rad = -lat_rad
            s = -s

    d2 = 1 - e2 * ss
    Rn = a / sqrt(d2)
    Rm = (1 - e2) * Rn / d2
    rf = (1 - e2) * Rn
    u = w - Rn * c
    v = z - rf * s
    f = c * u + s * v
    m = c * v - s * u
    p = m / (Rm + f)

    lat_rad += p

    ht = f + m * p / 2

    return (lat_rad, lon_rad, ht)


def procesa(est):
    global est_orbitas
    est_orbitas = {}
    Xest = estaciones[est][0]
    Yest = estaciones[est][1]
    Zest = estaciones[est][2]
    for satname in orbitas:
        est_orbitas[satname]=[]
        for horalist in orbitas[satname]:
            deltaX = horalist[1] * 1000 - Xest
            deltaY = horalist[2] * 1000 - Yest
            deltaZ = horalist[3] * 1000 - Zest
            dist = sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ)
            xprima = R[est][0][0] * deltaX + R[est][0][1] * deltaY + R[est][0][2] * deltaZ
            yprima = R[est][1][0] * deltaX + R[est][1][1] * deltaY + R[est][1][2] * deltaZ
            zprima = R[est][2][0] * deltaX + R[est][2][1] * deltaY + R[est][2][2] * deltaZ
            A = atan2(yprima,xprima)
            Z = acos(zprima/dist)
            est_orbitas[satname].append([horalist[0],xprima,yprima,zprima,dist,A,Z])
    return


def plotea(estacion):
    fig = go.Figure()
    for satname in satlist:
        traza = est_orbitas.get(satname)
        fig.add_trace(go.Scatter3d(
            x=[row[1] for row in traza],
            y=[row[2] for row in traza],
            z=[row[3] for row in traza],
            name=satname))
    fig.update_traces(marker_size = 2)
    fig.add_trace(go.Scatter3d(x=[0],y=[0],z=[0], name=estacion, marker=dict(color='black', size=5, sizemode='diameter')))
    fig.show()
    return


pol_est = {}
geo_est = {}
R = {}
orbitas = {}
est_orbitas = {}

filename = "igs11060.sp3"
satlist = []
with open(filename,encoding="utf-8") as f:
    for line in f:
        if line.startswith("*  2001"):
            seg = line[19:].split('.')
            fechaHora =  datetime(int(line[1:7]),int(line[7:10]),int(line[10:13]),int(line[13:16]),int(line[16:19]),int(seg[0]),int(seg[1]))
        elif line.startswith("P"):
            sat = line[1:4]
            x = float(line[5:19])
            y = float(line[19:33])
            z = float(line[33:47])
            fila = [fechaHora,x,y,z]
            if not(sat in satlist):
                satlist.append(sat)
                orbitas[sat]=[]
            orbitas[sat].append(fila)

for est in estaciones:
    x = estaciones[est][0]
    y = estaciones[est][1]
    z = estaciones[est][2]
    p = sqrt(x*x+y*y)
    r = sqrt(x*x+y*y+z*z)
    λ = atan2(y,x)  # declinacion = longitud
    φ = atan2(z,p)  # inclinacion
    pol_est[est] = [φ,λ,r]
    coord = ecef_to_geodetic(x,y,z)
    Ҩ = coord[0]  # latitud
    lon = coord[1]
    h = coord[2]
    geo_est[est] = [Ҩ,lon,h]
    R[est] = [  [-sin(Ҩ)*cos(λ),-sin(Ҩ)*sin(λ),cos(Ҩ)],
                [-sin(λ)       , cos(λ)       ,0     ],
                [ cos(Ҩ)*cos(λ), cos(Ҩ)*sin(λ),sin(Ҩ)]]

respuesta=''
while 1:
    print("Selecciones estación")
    print("1 - La Plata")
    print("2 - Trelew")
    print("3 - Rio Grande")
    #print("4 - Lat, Lon")
    print("x - Exit")
    if platform.system() == "Linux":
        respuesta = readchar.readchar().decode('utf-8')
    elif platform.system() == "Windows":
        respuesta = msvcrt.getch().decode('utf-8')
    if respuesta == "1":
        print()
        print("La PLata")
        procesa("LPGS")
        plotea("LPGS")
    elif respuesta == "2":
        print()
        print("Trelew")
        procesa("TREL")
        plotea("TREL")
    elif respuesta == "3":
        print()
        print("Rio Grande")
        procesa("RIO2")
        plotea("RIO2")
    elif respuesta == "4":
        print()
        print("Ingrese Latitud, Longitud en Grados y Decimales")
        lat = float(input("Latitud ? :"))
        lon = float(input("Longitud? :"))
        print()
        print("Latitud :",lat)
        print("Longitud:",lon)
        print()
        #procesa("OTRO")
    elif respuesta == "x":
        print()
        print("Exit")
        break
    else:
        print()
        continue
