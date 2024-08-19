#= Practica 2-i
    para el punto i) ..... ??

=#

#=
import os
import platform
=#

using Printf
using Format
using Dates
using LinearAlgebra
using Statistics

c = 299792458          # m/s  de ITRF
μ = 3.986005E14        # m3/s2  Earth gravitational constant
ωe =   7.2921151467E-5 # radians/s Angular Velocity of the Earth
ωE = [ 0, 0, ωe]       # Velocity of the Earth (vector)
tGPS0 =  DateTime(1980,1,6,0,0,0)

L = []
A = []
C = []
satord=Dict()        # este es sólo para mostrar los satélites con su coorespondiente correccion


# Observables incluyendo código P
#          C/A          Fase L1           Fase L2           P1           P2
P2D = Dict(
"28" => [23334619.807 ,-10956241.60549 , -8525575.42946 ,23334619.277 ,23334623.308],
"13" => [22586189.129 ,-12029006.00949 , -9358589.69746 ,22586189.572 ,22586192.629],
"01" => [25167667.280 ,   -19838.71849 ,   275138.51344 ,25167667.160 ,25167670.651],
"27" => [20873169.266 ,-23420787.65349 ,-18223955.59247 ,20873168.733 ,20873173.057],
"24" => [23371141.291 , -8542099.54349 , -6138557.63846 ,23371140.735 ,23371147.385],
"10" => [21505922.486 ,-16684416.38149 ,-12470663.29047 ,21505921.442 ,21505925.535],
"08" => [20958408.428 ,-22772002.73249 ,-17730879.20747 ,20958407.438 ,20958412.414]
)

f1 = 1575.42
f12 = f1 ^2
f2 = 1227.60
f22 = f2 ^2
fdif = f12 - f22
PD =  Dict(key => (f12 * P2D[key][4] - f22 * P2D[key][5])/fdif for key in keys(P2D))

# Obtengo tiempo de tránsito [s]
TT = Dict(key => PD[key]/c for key in keys(PD))

cant_sat = length(PD)


Precisas = Dict(      # Efemérides Precisas [Km] , [us]
    #        X             Y              Z           clk
"01" => [  581.886423, 25616.666528 ,  7088.545471, 169.092800],
"08" => [22018.953984,  2878.718252 , 14451.124018,   9.709180],
"10" => [10103.948910, -10925.429662, 22009.912003,   1.148951],
"13" => [ 7525.432597,  20488.591201, 15216.097471,  -0.655216],
"24" => [22368.646126, -12657.086060,  6934.928617,  36.698468],
"27" => [15057.427636,   9402.947329, 20171.667340,  14.763242],
"28" => [-5895.039751,  14576.928529, 21538.074040,  14.267922]
)

Delta_ant = Dict(     # Corrección de centro de fase de antena
"28" => 1.04280 ,     #     para cada sat [m]
"13" => 1.38950 ,
"01" => 2.38080 ,
"27" => 2.63340 ,
"24" => 2.60380 ,
"10" => 2.54650 ,
"08" => 2.57810 
)

Estacion = [          # Coord. precisas de la Estación [m]
    3370658.6942,     # X
    711877.0150,      # Y
    5349786.8637]     # Z


# Para tomar como coordenadas a-priori de la estacion.
#    se le agrega una diferencia a las precisas
Coord = [(i+(rand()-0.5)*5000) for i in Estacion]
#Coord = Estacion

# en la segunda iteracion el error de reloj va a estar estimado,
# para lo cual se va a necesitar agregar un elemento
# al vector posición
push!(Coord,0)


#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -
#  Levanto las efemérides transmitidas . . .
# ------------------------------------------------------------------------------------
filename = "Practica_2\\ifag0780.01n"
#mensajes = Dict()
#satlist = []

# para hacerla fácil, anoto el tiempo correspondiente a las
#  efemérides transmitidas que voy a usar para cada satélite ...
#tefem = Dict()
#sat=""
#sata = sat
open(filename) do f
    global mensajes = Dict()
    global satlist = []
    global GPSweek = 0.0

    # para hacerla fácil, anoto el tiempo correspondiente a las
    #  efemérides transmitidas que voy a usar para cada satélite ...
    global tefem = Dict()

    linea = 0
    sat=""
    sata = sat
    fechaHora = DateTime(2001,3,19,0,15,0)
    a0 = a1 = a2 = 0.0
    Toe = sqrtA = Delta_n = GPS_Week = 0.0
    e = M0 = i0 = idot = 0.0
    omega = OMEGA = OMEGA_DOT = 0.0
    Cuc = Cus = Crc = Crs = Cic = Cis = 0.0
    start = false
    for line in readlines(f)
        if start
            if line[6]!='.'
                if sat ≠ "" 
                    if sat ∉ satlist
                        push!(satlist,sat)
                        #print("Nuevo satelite: ",sat)
                        mensajes[sata] = []
                        tefem[sata] = fechaHora
                    end
                    push!(mensajes[sata],Dict("FechaHora"=> fechaHora,"a0"=> a0 ,"a1"=> a1,"a2"=> a2,
                                    "T0e"=> Toe,"GPSweek"=> GPS_Week,"sqrtA" => sqrtA,"e"=> e, "M0"=> M0,
                                    "omega"=>omega,"i0"=> i0,"OMEGA"=> OMEGA,"Delta_n"=> Delta_n,
                                    "idot"=> idot,"OMEGA_DOT"=> OMEGA_DOT,"Cuc"=> Cuc,"Cus"=> Cus,
                                    "Crc"=> Crc,"Crs"=> Crs,"Cic"=> Cic,"Cis"=> Cis))
                    GPSweek = GPS_Week
                end    
                sat = line[1:2]
                if sat[1]==' '
                    sata="0"*sat[2]
                else
                    sata=sat
                end
                #print(sat,sata)
                seg = split(line[19:22],'.')
                
                if seg[1]=="  "
                    seg[1]="0"
                end    
                AA = line[3:5]
                if parse(Int64,AA)<10
                    AA = "0" * line[5]
                    #print(AA)
                    #print(seg)
                end    
                #=
                #                             AA           MM              DD             HH
                fechaHora =  datetime(int("20"+AA),int(line[6:8]),int(line[9:11]),int(line[12:14])
                                    #         mm             ss
                                    ,int(line[15:17]),int(seg[0]),int(seg[1]))
                =#
                #                                  AA                    MM                   DD
                fechaHora = DateTime(parse(Int64,"20"*AA),parse(Int64,line[6:8]),parse(Int64,line[9:11])
                                    #               HH                       mm                  ss
                                    ,parse(Int64,line[12:14]),parse(Int64,line[15:17]),parse(Int64,seg[1]),parse(Int64,seg[2]))
                a0 = parse(Float64,replace(line[23:41],'D'=>'E'))
                a1 = parse(Float64,replace(line[42:60],'D'=>'E'))
                a2 = parse(Float64,replace(line[61:79],'D'=>'E'))
                linea = 0
            elseif startswith(line,"   ")
                linea += 1
                if linea==1
                    IODE   = parse(Float64,replace(line[ 4:22],'D'=>'E'))
                    Crs    = parse(Float64,replace(line[23:41],'D'=>'E'))
                    Delta_n= parse(Float64,replace(line[42:60],'D'=>'E'))
                    M0     = parse(Float64,replace(line[61:79],'D'=>'E'))
                elseif linea==2
                    Cuc   = parse(Float64,replace(line[ 4:22],'D'=>'E'))
                    e     = parse(Float64,replace(line[23:41],'D'=>'E'))
                    Cus   = parse(Float64,replace(line[42:60],'D'=>'E'))
                    sqrtA = parse(Float64,replace(line[61:79],'D'=>'E'))
                elseif linea==3
                    Toe   = parse(Float64,replace(line[ 4:22],'D'=>'E'))
                    Cic   = parse(Float64,replace(line[23:41],'D'=>'E'))
                    OMEGA = parse(Float64,replace(line[42:60],'D'=>'E'))
                    Cis   = parse(Float64,replace(line[61:79],'D'=>'E'))
                elseif linea==4
                    i0       = parse(Float64,replace(line[ 4:22],'D'=>'E'))
                    Crc      = parse(Float64,replace(line[23:41],'D'=>'E'))
                    omega    = parse(Float64,replace(line[42:60],'D'=>'E'))
                    OMEGA_DOT= parse(Float64,replace(line[61:79],'D'=>'E'))
                elseif linea==5
                    idot    = parse(Float64,replace(line[ 4:22],'D'=>'E'))
                    c2      = parse(Float64,replace(line[23:41],'D'=>'E'))
                    GPS_Week= parse(Float64,replace(line[42:60],'D'=>'E'))
                    L2_P    = parse(Float64,replace(line[61:79],'D'=>'E'))
                elseif linea==6
                    SVa  = parse(Float64,replace(line[ 4:22],'D'=>'E'))
                    SVh  = parse(Float64,replace(line[23:41],'D'=>'E'))
                    TGD  = parse(Float64,replace(line[42:60],'D'=>'E'))
                    IODC = parse(Float64,replace(line[61:79],'D'=>'E'))
                elseif linea==7
                    TxToM  = parse(Float64,replace(line[ 4:22],'D'=>'E'))
                    FitInt = parse(Float64,replace(line[23:41],'D'=>'E'))
                end
            end
        end    
        if occursin("LEAP SECONDS",line)
            L_sec = parse(Int64,line[1:8])
            #print(L_sec)
        end    
        if occursin("END OF HEADER",line)
            start = true
        end
    end
end                    
tGPS0 += Day(GPSweek*7)
# ------------------------------------------------------------------------------------
#  / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / -


function newton(fn,Df,x0,epsilon,max_iter)
    #=Approximate solution of f(x)=0 by Newton's method.

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
    =#
    xn = x0
    for _ in 0:max_iter
        fxn = fn(xn)
        if abs(fxn) < epsilon
            #print('Found solution after',n,'iterations.')
            return xn
        end    
        Dfxn = Df(xn)
        if Dfxn == 0
            print("Zero derivative. No solution found.")
            return None
        end    
        xn = xn - fxn/Dfxn
    end    
    print("Exceeded maximum iterations. No solution found.")
    return None
end

function coord_obs()
    #= Arma la matriz con coordenadas de cada satélite
         en el tiempo de observación
    =#
    Calc = Dict()
    for sati in keys(TT)
        s = mensajes[sati]
        for h in s
            if h["FechaHora"]==tefem[sati]
                #tefem_seg = int((tefem[sati]-tGPS0).total_seconds())
                tefem_seg = convert(Int64,(Dates.DateTime(tefem[sati]) - Dates.DateTime(tGPS0)) / Millisecond(1000))
                a = h["sqrtA"] * h["sqrtA"]
                a_cubo = a*a*a
                # tiempo entre la efemeride transmitida y
                # el momento en que quiero saber las coordenadas del satélite
                Delta_t = (tobs_seg - TT[sati]) - tefem_seg

                M = h["M0"] + (sqrt(μ/a_cubo)+h["Delta_n"])*Delta_t
                #=
                fE  = lambda x: M + h["e"]*sin(x) - x
                dfE = lambda x: h["e"]*cos(x)-1
                =#
                fE(x)  = M + h["e"]*sin(x) - x
                dfE(x) = h["e"]*cos(x)-1
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

                xyz = [  r*cos(u)   0   0;
                         r*sin(u)   0   0;
                             0      0   0]

                R31 = [  cos(-Ω) sin(-Ω)  0;
                        -sin(-Ω) cos(-Ω)  0;
                              0       0   1]

                R32 = [   cos(ϴ) sin(ϴ)   0;
                         -sin(ϴ) cos(ϴ)   0;
                              0      0    1]

                R1 =  [ 1       0         0;
                        0  cos(-i)  sin(-i);
                        0 -sin(-i)  cos(-i)]
                R3R3 = R32 * R31
                RR   = R3R3 * R1

                xyz_prima  = RR * xyz
                x = xyz_prima[1][1]
                y = xyz_prima[2][1]
                z = xyz_prima[3][1]
                Calc[sati] = [x, y, z, Precisas[sati][4]]
            end
        end
    end       
    return Calc
end

function arma_matriz()
    #=     al armar la matriz de diseño se puede poner c * Delta_t en la
            4ta. columna en vez de C, así los resultados 
            dan en Distancia, en lugar de Delta_t.
            Para eso, pongo todos unos en vez de c,
            así la incógnita incluye a c y quedan números más manejables
    =#
    global C
    global satord
    #global L = []
    global L = zeros(1,cant_sat)
    global A = zeros(Float64,4,cant_sat)
    C = []
    satord=Dict()
    j = 1
    for st in keys(Calculadas)
        s= Calculadas[st]
        rs = s[1:3]       # Coordenadas (x, Y, Z) calculadas del satélite "st"
        dX = Coord[1] - s[1]
        dY = Coord[2] - s[2]
        dZ = Coord[3] - s[3]
        ρ = sqrt(dX*dX+dY*dY+dZ*dZ)
        ρSagnac = dot((rr-rs),cross(ωE, rr)) / c

        @printf("Sat: %s  Coor. x Sagnac: %10.5f\n",st,ρSagnac)

        fila = Vector{Float64}([dX/ρ,dY/ρ,dZ/ρ,1])  # opcion con incognita c * Delta_t
        #println(fila)
        A[:,j] = fila
        # diferencia Observado - Calculado
        err = PD[st] - ρ - float(Coord[4]) + c * s[4] / 1E6  - ρSagnac
        # Si s[4] > 0 el satélite atrasa con respecto a GPS time, entonces "sumo" error en distancia
        L[1,j] = err
        #push!(L,err)
        satord[st] = err
        linea_C =[0 for i in 1:cant_sat]
        #linea_C[j]= err * err
        j +=1
        push!(C,linea_C)
    end
end    

function imprime_dif_sat()
    for s in keys(Calculadas)
        println("\nSat: ",s)
        linea1 =  " Calculadas : "
        linea2 =  " Precisas   : "
        linea3 =  " Diferencia : "
        mod = 0
        for j in 1:3
            linea1 = linea1 * @sprintf("  %15.5f m",Calculadas[s][j])
            linea2 = linea2 * @sprintf("  %15.5f m",Precisas[s][j]*1000)
            d = Calculadas[s][j]-Precisas[s][j]*1000
            linea3 = linea3 * @sprintf("  %15.5f m",d)
            mod += d * d
        end    
        println(linea1)
        println(linea2)
        println(linea3)
        println(@sprintf("  Módulo de la diferencia: %15.5f m",sqrt(mod)))
    end    
end

function imprime_resu()
    #=Imprime resultados parciales
    =#

    println("\nObservado-calculado")
    println("SAT         Dif\n")
    for j in keys(satord)
        @printf("%s  %15.6f\n",j,satord[j])
    end    
    println()
    @printf("Dispersión:  %8.6f\n\n",(Statistics.std(L)))

    println(" - - - - - - - - - -")
    println("Diferencias Calculadas para corregir las coordenadas")
    println("           X                   Y                   Z                   t")
    linea =""
    for i in X1
        linea *= @sprintf("%10.15f  ",i)
    end
    println(linea)
    dift = X1[1:3]
    @printf("\nMódulo de la corrección calculada: %20.16f m\n",norm(dift))
    push!(dift,X1[4]*c)
    @printf("Módulo de la corrección calculada\n              incluyendo el reloj: %10.6f m\n\n",norm(dift))
    println(" - - - - - - - - - -\n")
end

function imprime_Correg()
    #= Imprime una vez recalculado
    =#
    println("Coord Corregida,    Precisa,                Diferencia (Calculada - Precisa)")
    j = 1
    acu = 0
    for i in Coord
        if j<4
            linea = @sprintf("%15.5f  %15.5f  -->%15.5f m",i,Estacion[j],i-Estacion[j])
            acu += (i-Estacion[j])*(i-Estacion[j])
        else
            linea =  @sprintf("\n Dif. entre sitios: %15.5f m\n",sqrt(acu))
            linea *=   @sprintf(" Delta_t:           %15.5f useg",i/c*1E6)
            linea *= @sprintf("\n c * Delta_t:       %15.5f m\n",i)
        end    
        j +=1
        println(linea)
    end
end

function calcula_angulo_dif()
    vect_dif = []
    Coord_est = Estacion[1:3]
    for j in 1:3
        push!(vect_dif,Coord[j]-Estacion[j])
    end    
    modprod = norm(Coord_est)*norm(vect_dif)
    return acos(dot(Coord_est,vect_dif)/modprod)*180/pi
end



############################### -  -  -  -  -  -  -  -  -  -  -  -  -
#   I N I C I O    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
############################### -  -  -  -  -  -  -  -  -  -  -  -  -

if Sys.islinux()
    # import readchar # type: ignore
    run(`clear`)
elseif Sys.iswindows()
    # import msvcrt   # type: ignore
    run(`cmd /C cls`)
end

# el momento en que se realiza la observación
tobs = DateTime(2001,3,19,0,15,0)
#tobs_seg = int((tobs-tGPS0).total_seconds())
tobs_seg = convert(Int64,(Dates.DateTime(tobs) - Dates.DateTime(tGPS0)) / Millisecond(1000))


# consigo las coordenadas de cada satélite
#    en el tiempo de emisión de la señal
#    para usarlas en lugar de las precisas
Calculadas = coord_obs()

# Muestra las diferencias entrre las coordendas precisas del satélite
# al instante de observación y las calculadas al instante de emisión
imprime_dif_sat()
println("\n\n--------------------------------------------------------")
imprime_Correg()   # Primero muestra la condición inicial desde donde partimos

for paso in 1:500
    if paso == 2
        global inicio = now()
    end
    println("--------------------------------------------------------")
    println(@sprintf("----> Paso: %4d",paso))
    println()
    global rr = Coord[1:3]   # Coordenadas (x, Y, Z) calculadas de la estación
    arma_matriz()

    #P = inv(C)
    #MA = real(A')
    MA = A'
    ML = L'
    global X1 = MA\ML
    #println()
    imprime_resu()
    global Coord=[(Coord[i] + X1[i]) for i in eachindex(Coord)]
    imprime_Correg()
    @printf "Angulo: %8.3fº\n\n" calcula_angulo_dif()
end

final = now()
tproceso=(Dates.DateTime(final) - Dates.DateTime(inicio)) / Millisecond(1000)
println("\n------------------------------------------------------")
@printf("\ntiempo de proceso: %8.5fseg\n\n\n",tproceso)
