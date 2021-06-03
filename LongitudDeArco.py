# Mauricio Torres Díaz.
# No. de Control: 19161426
# Metodos Numericos: Longitud de arco.
# Programa que lee un polinomio de grado > 1 y determina su longitud de arco por el metodo de simpson/trapecio.

from math import sqrt 

#METODO PARA SACAR XI
def Xi(a,b,deltaX,n):
    for i in range(n+1):
        x.append(a+i*deltaX)

#FUNCION QUE SACA LA DERIVADA DEL POLINOMIO
def derivar(coeficientes):
    derivada = []
    k = 1
    if((len(coeficientes)) == 1):
        derivada.append(0)
    else:
        while(len(coeficientes)>k): 
            derivada.append(coeficientes[k]*k)
            k=k+1
    return derivada 

#METODO DE HORNER
def horner(grado, derivada, x):
    polinomio = 0
    for i in range(grado):
        polinomio = polinomio + (derivada[i] *pow(x,i))
    polinomio = sqrt(1+pow(polinomio,2))
	# Al término de este ciclo WHILE, la variable polinomio tiene el valor del P(x)
    return polinomio

#METODO DE LA REGLA DEL TRAPECIO
def Trapecio(a,b,n,tolerancia, valor_anterior,derivada):
    suma = integral = 0
    deltaX = (b-a)/n # Variable DetlaX

    #Realizamos la suma de Riemman
    for i in range(1, n): #Iniciamos en 1 hasta n
        suma += horner(grado,derivada,(a+(i*deltaX)))

    #Realizamos la regla del Trapecio
    integral = ((b-a)/(2*n))*(horner(grado, derivada, a)+ 2*suma +horner(grado, derivada, b))
    error = ((integral - valor_anterior)/integral)*100;
    print("Aproximacion de la integral: "+str(integral)+" fue obtenida con n = "+str(n)+". El error es: "+str(error))

    if(abs(error)>=tolerancia):
        n +=10
        valor_anterior = integral
        return Trapecio(a,b,n,tolerancia,valor_anterior,derivada)

    print("Resultado de la integral: "+str(integral))
    print("Calculado con "+str(n)+" subintervalos.")

#METODO DE LA REGLA DE SIMPSON
def Simpson(a,b,n,tolerancia,valor_anterior,derivada):
    suma1 = suma2 = integral = 0
    deltaX = (b-a)/n #delta x

    #realizamos a la suma de riemmann
    for i in range(1, n): #iniciamos en 1 hasta n
        if i%2 == 0:
            suma2 += horner(grado,derivada,(a+(i*deltaX)))
        else:
            suma1 += horner(grado,derivada,(a+(i*deltaX)))
    
    #formula de la regla de simpson
    integral = (deltaX/3)*(horner(grado, derivada, a)+ 4*suma1 + 2*suma2 +horner(grado, derivada, b))
    error = ((integral - valor_anterior)/integral)*100;
    print("La aproximacion de la integral es: "+str(integral)+" tomando en cuenta n = "+str(n)+". El error es: "+str(error))

    if(abs(error)>=tolerancia):
        n +=10
        valor_anterior = integral
        return Simpson(a,b,n,tolerancia,valor_anterior,derivada)

    print("Resultado de la integral: "+str(integral))
    print("Calculado con "+str(n)+" subintervalos.")

#METODO DE LONGITUD DE ARCO
def longitudDeArco():
    derivada = derivar(coeficientes)
    print(derivada)
    #REALIZA EL METODO CORREPONDIENTE DEPENDIENDO QUE HAYA PEDIDO EL USUARIO
    if opcion == 1:
        Simpson(a,b,n,tolerancia,valor_anterior,derivada)
    elif opcion == 2:
        Trapecio(a,b,n,tolerancia,valor_anterior,derivada)

#............................................................................................................
#FUNCIONES PRINCIPALES 

#Pedimos el grado del polinomio y este debe ser mayor a 1
grado = 0
print("Este programa te ayudara a obtener la Longitud de Arco.")
#El grado del polinomio tieen que ser mayor a 1
while grado<=0:
    grado = int(input("¿Cual es el grado de tu polinomio? (Tiene que ser mayor a 1): ")) 
    if (grado < 2):
        print("ERROR.:( \nNo te puedo ayudar con ese Polinomio, porfavor ingrese un Polinomio valido.")

#Array donde se guardaran los coeficientes de nuestro polinomio
coeficientes = []
x =[]

# Ingresamos los coeficientes del polinomio comenzando con el termino de mayor grado y terminando 
print("Porfavor, ingrese los coeficeintes del Polinoio (Empezando por el de mayor grado) ")
for i in range(grado+1):
    coeficiente = float(input("Ingresa el coeficiente: "))
    coeficientes.append(coeficiente)
coeficientes.reverse()

a = b = 0
# a y b no pueden ser iguales
while a == b:
    a = float(input("\nLimite inferior: "))
    b = float(input("Limite superior: "))
    if(a == b):
        print("Los limites no pueden ser iguales, porfvaor ingrese los limites bien.")

n = 6 #Le damos 6 a n, ya que es un numero que puede aproximar los metodos de la mejor manera

#Pedimos las cifras significativas
cifrasSig = int(input("Cifras significativas: "))
tolerancia = 0.5*(10**(2-cifrasSig))
valor_anterior = 0
opcion = 0
print("\nMENU: \n1.- Metodo Simpson\n2.- Metodo Trapecios \n Escriba el número del metodo")
while(opcion != 1 or opcion != 2):
    opcion = int(input("¿Con qué metodo quieres resolverlo?: "))
    if(opcion == 1 or opcion == 2):
        print("\nLa tolerancia es "+str(tolerancia)+" para obtener "+str(cifrasSig)+" cifras significativas.\n")
        longitudDeArco()
        print("\nEspero haberte ayudado.:D")
        break
    else: 
        print("Porfavor, ingrese una opcion valida.")