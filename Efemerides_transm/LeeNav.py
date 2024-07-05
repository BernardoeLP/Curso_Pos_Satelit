from datetime import datetime
from math import sin, cos, asin, acos, atan2, sqrt, pi
#   3.986004418E14   m3/s2
μ = 3986004.418E08 # m3/s2  Earth gravitational constant
ωe =   7292115E-11 # radians/s Angular Velocity of the Earth
f = 1 / 298.257223563
e2 = f*(2-f)
e= sqrt(e2)


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
    for n in range(0,max_iter):
        fxn = fn(xn)
        if abs(fxn) < epsilon:
            print('Found solution after',n,'iterations.')
            return xn
        Dfxn = Df(xn)
        if Dfxn == 0:
            print('Zero derivative. No solution found.')
            return None
        xn = xn - fxn/Dfxn
    print('Exceeded maximum iterations. No solution found.')
    return None

filename = "Efemerides_transm\\TREL1610.24n"
filename = "Efemerides_transm\\ifag0770.01n"
filename = "Efemerides_transm\\ifag0780.01n"
filename = "Efemerides_transm\\TREL1630.24n"
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
                                    "T0e": Toe,"sqrtA" : sqrtA,"e": e, "M0": M0,
                                    "omega":omega,"i0": i0,"OMEGA": OMEGA,"Delta_n": Delta_n,
                                    "idot": idot,"OMEGA_DOT": OMEGA_DOT,"Cuc": Cuc,"Cus": Cus,
                                    "Crc": Crc,"Crs": Crs,"Cic": Cic,"Cis": Cis})
                sat = line[:2]
                seg = line[18:22].split('.')
                if seg[0]=='  ':
                    seg[0]='0'
                AA = line[3:5]
                if int(AA)<10:
                    AA = "0" + line[4]
                    print(AA)
                    print(seg)
                #                                 AA             MM               DD                HH           mm             ss
                fechaHora =  datetime(int("20"+AA),int(line[6:8]),int(line[9:11]),int(line[12:14]),int(line[15:17]),int(seg[0]),int(seg[1]))
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
            print(L_sec)
        if "END OF HEADER" in line:
            start = True

position = {}
t=datetime(2001,3,19,0,0,0)
#print(mensajes)

"""
for sati in mensajes:
    s= mensajes[sati]
    for h in s:
        print("Satélite: "+ sati)
        #print(h)
        #print()
        a = h["sqrtA"] * h["sqrtA"]
        a_cubo = a*a*a
        M = h["M0"] + (sqrt(μ/a_cubo)+h["Delta_n"])*(t-h["T0e"])
        fE  = lambda x: M +e*sin(x) - x
        dfE = lambda x: e*cos(x)-1
        E = newton(fE,dfE,0,1E-5,3)
        print(E)
"""
