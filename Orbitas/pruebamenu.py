import msvcrt

est = "LPGS"
respuesta=''

while 1:
    print("Selecciones estación")
    print("1 - La Plata")
    print("2 - Trelew")
    print("3 - Rio Grande")
    print("4 - Lat, Lon")
    print("x - Exit")
    respuesta = msvcrt.getch().decode('utf-8')
    # print(respuesta)
    if respuesta == "1":
        print()
        print("La PLata")
    elif respuesta == "2":
        print()
        print("Trelew")
    elif respuesta == "3":
        print()
        print("Rio Grande")
    elif respuesta == "4":
        print()
        print("Ingrese Latitud, Longitud en Grados y Decimales")
        lat = float(input("Latitud ? :"))
        lon = float(input("Longitud? :"))
        print()
        print("Latitud :",lat)
        print("Longitud:",lon)
        print()
    elif respuesta == "x":
        print()
        print("Exit")
        break
    else:
        print()
        continue


