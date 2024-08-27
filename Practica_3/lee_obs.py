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

filename = "Practica_3\\zimm1440.02o.txt"
start = False
mensajes = {}
satlist = []
linea = 0
sat=""
with open(filename,encoding="utf-8") as f:
    for line in f:
        if start:
            if line[18]!='.':

                if sat != "":

                    if sat not in satlist:
                        satlist.append(sat)
                        mensajes[sat] = []



                nsats= int(line[30:32])
                sats = line[32:]
                seg = line[18:22].split('.')
                if seg[0]=='  ':
                    seg[0]='0'
                AA = line[:3]
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
                    fechaHora =
                    C1 =
                    P2 =
                    L1 = 
                    L2 =
                    D1 =
                if linea==2:
                    mensajes[sat].append({"FechaHora": fechaHora","C1":C1,"P2":P2,"L1":L1,"L2":L2,"D1":D1 ,"D2":D2"})
                    linea = 0

        if "# OF SATELLITES" in line:
            N_sat = int(line[:8])
            #print(L_sec)
        if "END OF HEADER" in line:
            start = True
