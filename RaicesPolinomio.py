#Mauricio Torres Diaz.      No. De Control: 19161426
#Programa que lee un polinomio de grado n > 0 y localiza puntos notables con el metodo de bairstown

from math import pow, inf
import cmath as cmath
import math

#declaracion de metodos

#Vamos a usar el algoritmo de horner para hacer la division de polinomios 
#horner devuelve el valor del residuo de dividir el polinomio entre un número dado (o sea, el residuo)
#horner tambien sirve para evaluar un polinomio en cierto valor de x

def horner(grado, listaCoeficientes, x):  
    polinomio= listaCoeficientes[grado]
    k = grado - 1
    while (k >= 0):
        polinomio= listaCoeficientes[k] + (polinomio*x)
        k = k - 1
        # Al término de este ciclo WHILE, la variable polinomio tiene el valor del P(x) 
    return polinomio


#devuelve un array con los coeficientes de la derivada del polinomio
def derivar(listaCoeficientes):
    #len(listaCoeficientes)-1 de una segunda derivada del polinomio grado 1, es igual a 0
    #pero debemos hacer que sea 1, así que ponemos un if
    derivada=[]
    i=1
    if((len(listaCoeficientes)) == 1):
        derivada.append(0)
    else:
        #grado del polinomio = len(listaCoeficientes)-1
        while(len(listaCoeficientes)>i): 
            derivada.append(listaCoeficientes[i]*i)
            i=i+1

    return derivada 

#Metodo de Bairstow)
def bairstow(coeficientes, r, s, grado, raiz, tolerancia):
	if(grado<1): #no hay raices
		return None
	if((grado==1) and (coeficientes[1]!=0)): #saca la única raíz que hay
		raiz.append((float(coeficientes[0])/float(coeficientes[1]))*(-1))
		return None
	if(grado==2):
		D = (coeficientes[1]**2.0)-(4.0)*(coeficientes[2])*(coeficientes[0])
		if(D<0):
			X1 = (-coeficientes[1] - cmath.sqrt(D))/(2.0*coeficientes[2])
			X2 = (-coeficientes[1] + cmath.sqrt(D))/(2.0*coeficientes[2])
		else:
			X1 = (-coeficientes[1] - math.sqrt(D))/(2.0*coeficientes[2])
			X2 = (-coeficientes[1] + math.sqrt(D))/(2.0*coeficientes[2])
		raiz.append(X1)
		raiz.append(X2)
		return None
	n = len(coeficientes)
	b = [0]*len(coeficientes)
	c = [0]*len(coeficientes)
	b[n-1] = coeficientes[n-1]
	b[n-2] = coeficientes[n-2] + r*b[n-1]
	i = n - 3
	while(i>=0):
		b[i] = coeficientes[i] + r*b[i+1] + s*b[i+2]
		i = i - 1
	c[n-1] = b[n-1]
	c[n-2] = b[n-2] + r*c[n-1]
	i = n - 3
	while(i>=0):
		c[i] = b[i] + r*c[i+1] + s*c[i+2]
		i = i - 1
	Din = ((c[2]*c[2])-(c[3]*c[1]))**(-1.0)
	r = r + (Din)*((c[2])*(-b[1])+(-c[3])*(-b[0]))
	s = s + (Din)*((-c[1])*(-b[1])+(c[2])*(-b[0]))
	if(abs(b[0])> tolerancia or abs(b[1])>tolerancia):
		return bairstow(coeficientes,r,s,grado,raiz, tolerancia)
	if (grado>=3):
		Dis = ((r)**(2.0))+((4.0)*(1.0)*(s))
		X1 = (r - (cmath.sqrt(Dis)))/(2.0)
		X2 = (r + (cmath.sqrt(Dis)))/(2.0)
		raiz.append(X1)
		raiz.append(X2)
		return bairstow(b[2:],r,s,grado-2,raiz, tolerancia)


#Una ves tenemos todos los metodos que usaremos a lo largo del programa
#Empezamos a hacer las operaciones pidiendo el polinomio
print("-----------------------------------------------------------------------------------------------------------------------")
print("\n Hola que tal, soy el programa que te ayudar a a calcular los puntos notables de tu polinomio, así que comencemos.:D")
print("-----------------------------------------------------------------------------------------------------------------------")

#se asegura de que el grado sea mayor a cero
grado = input("Porfavor, ingresa el grado del polinomio (debe ser mayor a cero): ")
print("\n------------------------------------------------------------------------------------------")
grado=int(grado)
while grado <= 0:
   grado=input("El grado debe ser mayor a cero, porfavor intentelo de nuevo.")
   grado=int(grado)

#Pide los coeficientes al usuario y se los asigna su respectiva x
listaCoeficientes=[]
print("Ahora me daras los valores de tu polinomio.:D")
for i in range (grado+1):
    coef=float(input("Ingresa x^"+str(grado-i)+": "))
    listaCoeficientes.append(coef)
listaCoeficientes.reverse()

#vas a usar el porcentaje de error relativo y las cifras significativas, para estimar el numero de iteraciones que harás en el algoritno de Newton-Raphson
#Pide las cifras significativas con las que hemos estado trabajando
print("\n------------------------------------------------------------------------------------------")
print("Con cuantas cifras significativas quieres trabajar.:D")
cifrasSignificativas=input("Ingrese las cifras significativas: ")
cifrasSignificativas=int(cifrasSignificativas)
cifras = 1/pow(10, cifrasSignificativas)

# La variable TOLERANCIA se calcula con la cantidad de cifras significativas.
tolerancia = 0.5*(10**(2-cifrasSignificativas))

# La variable ERNP guardará el error relativo. Se inicia en 1, para garantizar que es mayor que la tolerancia desde el inicio.
errorRNP = 1

# Primera y segunda derivada
derivada1 = derivar(listaCoeficientes)
derivada2 = derivar(derivada1)
r = listaCoeficientes[0]/listaCoeficientes[-1]
s = r
raiz1 = derivada1[0]/derivada1[-1]
#r2 no importa si el polinomio es de grado 1
if(grado==1):
    raiz2=1
else:
    raiz2 = derivada2[0]/derivada2[-1]

s1 = raiz1
s2 = raiz2

# Raíces de los polinomios
raices=[]
raices1=[]
raices2=[]

bairstow(listaCoeficientes, r, s, grado, raices, tolerancia)
bairstow(derivada1, raiz1, s1, grado-1, raices1, tolerancia) #raices de la primera derivada
bairstow(derivada2, raiz2, s2, grado-2, raices2, tolerancia) #raices de la segunda derivada

#las raices con +0j son normales y nos quedamos solo con la parte real
for i in range (len(raices)): 
    if(complex(raices[i]).imag==0):
        raices[i]=raices[i].real
for i in range (len(raices1)): 
    if(complex(raices1[i]).imag==0):
        raices1[i]=raices1[i].real
for i in range (len(raices2)): 
    if(complex(raices2[i]).imag==0):
        raices2[i]=raices2[i].real

raices1= [x for x in raices1 if type(x) is not complex] #si se vacía, no tiene max, min, ni crecec ni decrece
raices2 = [x for x in raices2 if type(x) is not complex]#si se vacía, no tiene puntos de inflexion ni concavidad

#Calculamos los puntos máximos y mínimos y los guarda en una lista
puntosMaximos = []
puntosMinimos = []

if(len(raices1)!=0):
    for x in raices1:
        izquierda = horner(grado-1, derivada1, (x-cifras))
        derecha = horner(grado-1, derivada1, (x+cifras))

        if ((izquierda < 0) and (derecha > 0)): #-_-
            puntosMinimos.append((x, horner(grado, listaCoeficientes, x)))

        elif ((izquierda > 0) and derecha < 0): #_-_
            puntosMaximos.append((x, horner(grado, listaCoeficientes, x)))
else: #si raices1 está vecío, no hay maximos ni minimos
    puntosMaximos.append("No hay máximos")
    puntosMinimos.append("No hay minimos")


#Calculamos los Puntos de inflexión (donde cambia de direccion la concavidad)(la x es la que lo dice más que nada) y los guardamos en su respectiva lista
puntosDeInflexion = []
puntosArriba = []
puntosAbajo = []

if(len(raices2)!=0):
    #cuando el grado <= 2, no tiene puntos de inflexion ni concavidad
    if(grado > 2):
        for x in raices2: #la coordenada x la saca de raices2
            #si raices2 está vacío, no habrá puntos de inflexion, y por tanto, tampoco concavidad
            puntosDeInflexion.append((x, horner(grado, listaCoeficientes, x))) #calcula la coordenada y

        #Intervalos de concavidad (tiene que ver con los puntos de inflexion)
        raices2.sort() #ordena

        #poner un i=0
        for i in range(1, len(raices2)):
            if horner(grado-2, derivada2, raices2[i]-cifras) < 0: #la concava va hacia abajo cuando la segunda derivada es negativa
                puntosAbajo.append((raices2[i-1], raices2[i]))
                #abajo.append(complex(raices2[i-1], raices2[i]))
            else:
                puntosArriba.append((raices2[i-1], raices2[i])) #lo contrario si es positiva
                #arriba.append(complex(raices2[i-1], raices2[i]))
    
        #ultimo punto de inflexion que apunta a infinito
        if horner(grado-2, derivada2, raices2[-1]+cifras) < 0:
            puntosAbajo.append((raices2[-1], inf))
            #abajo.append(complex(raices2[-1], inf))
        else:
            puntosArriba.append((raices2[-1], inf))
            #arriba.append(complex(raices2[-1], inf))

        #primer punto de inflexion que apunta desde -infinito
        if horner(grado-2, derivada2, raices2[0]-cifras) < 0:
            puntosAbajo.append((-inf, raices2[0]))
            #abajo.append(complex(-inf, raices2[0]))
        else:
            puntosArriba.append((-inf, raices2[0]))
            #arriba.append(complex(-inf, raices2[0]))
    else:
        puntosDeInflexion.append("No hay existen de inflexión en un polinomio de grado <= 2")
        puntosArriba.append("No hay existen de inflexión en un polinomio de grado <= 2")
        puntosAbajo.append("No hay puntos de inflexión en un polinomio de grado <= 2")
else:#cuando la lista raices2[] está vacía, no hay puntos de inflexion ni concavidades
    puntosDeInflexion.append("No hay puntos de inflexion")
    puntosArriba.append("No hay concavidad hacia arriba")
    puntosAbajo.append("No hay concavidad hacia abajo")



#Calculamos los intervalos donde crece o decrece
creciente = []
decreciente = []
raices1.sort() #ordena

#la función no crece ni decrece cuando el grado=1 (index out of range)
if(grado==1):
    creciente.append("La función no crece ni decrece cuando grado = 1")
    decreciente.append("La función no crece ni decrece cuando grado = 1")
else:
    if(len(raices1)!=0):
        #poner un i=0 o i=1
        for i in range(1, len(raices1)):#sumar +1
            if horner(grado-1, derivada1, raices1[i]-cifras) < 0:
                decreciente.append((raices1[i-1], raices1[i]))
            else:
                creciente.append((raices1[i-1], raices1[i]))

        if horner(grado-1, derivada1, raices1[-1]+cifras) < 0:     
            decreciente.append((raices1[-1], inf))
        else:
            creciente.append((raices1[-1], inf))

        if horner(grado-1, derivada1, raices1[0]-cifras) < 0:     
            decreciente.append((-inf, raices1[0]))
        else:
            creciente.append((-inf, raices1[0]))


    else:#si la lista raices1 está vacía, el polinomio no crece ni dcrece
        creciente.append("No crece")
        decreciente.append("No decrece")

#Imprimimos todos los valores obtenidos

#raices
print("\n------------------------------------------------------------------------------------------")
print("\nLas Raices del polinomio son:")
for i in raices:
    print(i)

#Puntos Maximos
print("\n------------------------------------------------------------------------------------------")
print("\nLos puntos maximos son:")
for x in puntosMaximos:
    print(x, end=", ")
#Putnos Minimos
print("\n------------------------------------------------------------------------------------------")
print("\nLos puntos minimos son:")
for x in puntosMinimos:
    print(x, end=", ")

#Puntos de inflexion
print("\n------------------------------------------------------------------------------------------")
print("\nLospuntos de inflexion son (coordenadas):")
for x in puntosDeInflexion:
    print(x, end=", ")

#crecientes
print("\n------------------------------------------------------------------------------------------")
print("\nEs creciente en:")
for x in creciente:
    print(x, end=", ")
#Decreciente
print("\n------------------------------------------------------------------------------------------")
print("\nEs decreciente en:")
for x in decreciente:
    print(x, end=", ")

#Concava hacia arriba
print("\n------------------------------------------------------------------------------------------")
print("\nEs concava hacia arriba en los intervalos:")
for x in puntosArriba:
    print(x, end=", ")
#Concava hacia abajo
print("\n------------------------------------------------------------------------------------------")
print("\nEs concava hacia abajo en los intervalos:")
for x in puntosAbajo:
    print(x, end=", ")

print("\n\nPues eso ha sido todo, agradezco hayas confiado en mí.:D. \nQue tengas un buen día.:D")
