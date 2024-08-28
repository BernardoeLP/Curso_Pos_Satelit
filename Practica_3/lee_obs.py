"""
Lee archivo RINEX de Observables:
     2              OBSERVATION DATA    G (GPS)             RINEX VERSION / TYPE
Zimmerwald LT88                                             COMMENT
GPS-BASE            SWISSTOPO                               OBSERVER / AGENCY
2691                TRIMBLE 4000SSI     7.29                REC # / TYPE / VERS
99390               TRM29659.00     NONE                    ANT # / TYPE
  4331297.3380   567555.6400  4633133.7180                  APPROX POSITION XYZ
        0.0000        0.0000        0.0000                  ANTENNA: DELTA H/E/N
     1     1     0                                          WAVELENGTH FACT L1/2
     6    C1    P2    L1    L2    D1    D2                  # / TYPES OF OBSERV
    30                                                      INTERVAL
  2002     5    24     0     0    0.000000                  TIME OF FIRST OBS
  2002     5    24    23    59   30.000000                  TIME OF LAST OBS
    28                                                      # OF SATELLITES

"""

from datetime import datetime
import os

os.system("cls")
filename = "Practica_3\\zimm1440.02o.txt"
start = False
mensajes = {}
satlist = []
linea = 0
sat=""
with open(filename,encoding="utf-8") as f:
    for line in f:
        if start:
            if (len(line)> 15) and (line[18]=='.'):
                #print(line[:-1])
                useg = str((float(line[16:26])-int(line[16:18]))* 1E6)
                seg = []
                seg.append(line[16:18].replace(' ','0'))
                seg.append(useg.split('.',maxsplit=1)[0])
                AA = line[1:3]
                AA.replace(' ','0')
                #print(seg)
                #                            AA             MM             DD
                fechaHora =  datetime(int("20"+AA),int(line[4:6]),int(line[6:9]),
                #                       HH              mm              ss         useg
                                int(line[9:12]),int(line[12:15]),int(seg[0]),int(seg[1]))
                #print(fechaHora)
                nsats_li= int(line[30:32])
                sats_li = line[32:-1]
                #print(sats_li)
                #print(nsats_li,len(sats_li)/3)
                s_list=[]
                for i in range(nsats_li):
                    s_list.append(sats_li[i*3:i*3+3])
                #print(s_list)
                #print('\n')
                for s in s_list:
                    if s not in satlist:
                        satlist.append(s)
                        mensajes[s] = []
                linea = 0
                s_count = 0
            elif line.startswith("  "):
                linea += 1
                sat = s_list[s_count]
                #print("Linea: ",linea)
                #print(line)
                if linea==1:
                    C1 = float(line[:16])
                    P2 = float(line[16:32])
                    L1 = float(line[32:48])
                    L2 = float(line[48:64])
                    D1 = float(line[64:-1])
                if linea==2:
                    D2 = float(line.strip())
                    mensajes[sat].append({"FechaHora": fechaHora,"C1":C1,
                                    "P2":P2,"L1":L1,"L2":L2,"D1":D1 ,"D2":D2})
                    linea = 0
                    s_count += 1


        if "# OF SATELLITES" in line:
            N_sat = int(line[:8])
            #print(L_sec)
        if "END OF HEADER" in line:
            start = True

for s in satlist:
    print("Satélite: ",s)

fs = open("Practica_3\\sat7.csv",'w')
fs.write("Fecha y Hora,C1,P2,L1,L2,D1,D2\n")
for hora in mensajes["G07"]:
    fila = hora["FechaHora"].strftime("%d/%m/%Y %H:%M:%S.%f")+",{:15.5f},{:15.5f},{:15.5f},{:15.5f},{:15.5f},{:15.5f}\n".format(
        hora["C1"], hora["P2"], hora["L1"], hora["L2"], hora["D1"], hora["D2"]) 
    fs.write(fila)
    #print(fila[:-1])
fs.close()
