from datetime import datetime

estaciones = {
'LPGS' : [float("2780089.996"),float("-4437398.174"),float("-3629387.406")],
'TREL' : [float("1938169.196"),float("-4229010.224"),float("-4348880.32")],
'RIO2' : [float("1429907.806"),float("-3495354.833"),float("-5122698.637")]
}



filename = "igs11060.sp3"
satlist = []
orbitas = dict()
with open(filename) as f:
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

#print(satlist)
#print(orbitas.get(satlist[0]))
