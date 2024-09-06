""" Practica 3
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

    

Y graba los correspondientes al satélite "07" en un archivo ".csv"

"""
# pylint: disable= C0209, C0301

import os
import platform
from datetime import datetime


############################### -  -  -  -  -  -  -  -  -  -  -  -  -
#   I N I C I O    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
############################### -  -  -  -  -  -  -  -  -  -  -  -  -
if platform.system() == "Linux":
    # import readchar # type: ignore
    os.system('clear')
elif platform.system() == "Windows":
    # import msvcrt   # type: ignore
    os.system('cls')



filename = "Practica_3\\zimm1440.02o.txt"
start = False
mensajes = {}  # Dictionary to store messages for each satellite
satlist = []   # List to store unique satellite identifiers
linea = 0      # Line counter for satellite data
sat = ""       # Current satellite identifier
with open(filename,encoding="utf-8") as f:
    for line in f:
        if start:
            if (len(line)> 15) and (line[18]=='.'):
                # Process the line containing satellite data
                useg = str((float(line[16:26])-int(line[16:18]))* 1E6)
                seg = []
                seg.append(line[16:18].replace(' ','0'))
                seg.append(useg.split('.',maxsplit=1)[0])
                AA = line[1:3]
                AA.replace(' ','0')
                # Create a datetime object from the parsed date and time components
                #                            AA             MM             DD
                fechaHora =  datetime(int("20"+AA),int(line[4:6]),int(line[6:9]),
                #                       HH              mm              ss         useg
                                int(line[9:12]),int(line[12:15]),int(seg[0]),int(seg[1]))
                nsats_li= int(line[30:32])  # Number of satellites in the line
                sats_li = line[32:-1]       # Satellite identifiers in the line
                s_list=[]
                for i in range(nsats_li):
                    s_list.append(sats_li[i*3:i*3+3])
                for s in s_list:
                    if s not in satlist:
                        satlist.append(s)
                        mensajes[s] = []
                linea = 0
                s_count = 0
            elif line.startswith("  "):
                linea += 1
                sat = s_list[s_count]
                if linea==1:
                    # Parse the first line of satellite data
                    C1 = float(line[:16])
                    P2 = float(line[16:32])
                    L1 = float(line[32:48])
                    L2 = float(line[48:64])
                    D1 = float(line[64:-1])
                if linea==2:
                    # Parse the second line of satellite data and store it
                    D2 = float(line.strip())
                    mensajes[sat].append({"FechaHora": fechaHora,"C1":C1,
                                    "P2":P2,"L1":L1,"L2":L2,"D1":D1 ,"D2":D2})
                    linea = 0
                    s_count += 1


        if "# OF SATELLITES" in line:
            N_sat = int(line[:8])  # Number of satellites in the header
        if "END OF HEADER" in line:
            start = True           # Start processing data after the header

for s in satlist:
    print("Satélite: ",s)

# Write the parsed data to a CSV file
fs = open("Practica_3\\sat7.csv",'w')
fs.write("FechaHora,C1,P2,L1,L2,D1,D2\n")
for entrada in mensajes["G07"]:
    fila = entrada["FechaHora"].strftime("%d/%m/%Y %H:%M:%S.%f")+",{:15.5f},{:15.5f},{:15.5f},{:15.5f},{:15.5f},{:15.5f}\n".format(
        entrada["C1"], entrada["P2"], entrada["L1"], entrada["L2"], entrada["D1"], entrada["D2"])
    fs.write(fila)
fs.close()
