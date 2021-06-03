
#Algoritmo de Gram Schmidt dada una matriz cuadrada
#Mauricio Torres Diaz
#Escribe un programa que lea una matriz cuadrada nxn encuentre sus eigenvalores utilizando el algoritmo QR.
#
#El programa tendrá las entradas

#Orden de la matriz
#Elementos de la matriz
#y entregará como resultado

#Los eigenvalores reales de la matriz.
#El programa deberá imprimir la factorización de las matrices Ai+1 = Qi+1Ri+1 que se obtienen en cada iteración.

#Se pedirá una captura de pantalla de alguna de las factorizaciones.

from math import sqrt #para raiz cuadrada

#declaracion de metodos-----------------------------------------------------------------

#Este metodo comprueba si la matriz que le damos es triangular superior
def esTriangularSuperior(matriz, n):
    cont=0
    respuesta=True
    i=1
    while(i<n):
        cont=cont+1
        for j in range(0,cont):
            if(matriz[i][j]>0.0000000001 or matriz[i][j]<-0.0000000001):
                respuesta=False
                break
        i=i+1
    return respuesta

#Metodo para Algoritmo QR
def algoritmoQR(matriz,matrizOrtogonal,matrizQ,matrizQt,n):
    #definir variables
    triangular=False
    contador = 0

    matrizR=[]
    matrizAi_mas_1=[]
    matrizAi=[]
    for i in range(0,n):
        matrizAi.append([0]*n)
        matrizAi_mas_1.append([0]*n)
        matrizR.append([0]*n)

    #paso 1: Ai=A
    for j in range(0,n):
        for i in range(0,n):
            matrizAi[i][j] =matriz[i][j]

    #paso 2: Repetir hasta que Ai+1
    #c)
    while(triangular==False):
        #a)
        #sacar la matrizQ de Ai
        gram(matrizAi, n, matrizOrtogonal)
        ortonormalizar(matrizOrtogonal, n, matrizQ)

        #sacar la matrizR de Ai
        #primero sacar la transpuesta de Q -> Qt
        hacerTranspuesta(matrizQt, n, matrizQ)
        #multiplicar R=(Qt)(A)
        multiplicacionMatrices(matrizQt, matrizAi, matrizR, n)

        #mostrar la factorización Ai=QiRi para cada iteracion
        print("\n A"+str(contador)+" = Q"+str(contador)+"*R"+str(contador))
        for i in range(0,n):
            for j in range(0,n):
                print("|",round(matrizAi[i][j], 5), end=" |") 
            if(i==0): #solo impribe una vez el = para que se va más bonito UwU
                print(" = ",end='')
            else:
                print("   ",end='')
            for j in range(0,n):
                print("|",round(matrizQ[i][j], 5), end=" |") 
            if(i==0): #soolo imprime una vez el * para que es vea más bonito uwu
                print("*",end='')
            else:
                print(" ",end='')
            for j in range(0,n):
                print("|",round(matrizR[i][j], 5), end=" |")
            print()
        #sacar Ai_mas_1 = (R)(Q)
        multiplicacionMatrices(matrizR, matrizQ, matrizAi_mas_1, n)

        #comprobamos si la matrizAi_mas_1 ya es diagonal superior
        if(esTriangularSuperior(matrizAi_mas_1,n)==True):
            triangular=True #si es así, salimos de las iteraciones
        else:#si no lo es, repetimos el proceso, ahora Ai_mas_1 es la nueva Ai
            contador = contador + 1
            for j in range(0,n):
                for i in range(0,n):
                    matrizAi[i][j] =matrizAi_mas_1[i][j]
                    #limpiamos la matriz resultado y la Ai_mas_1
                    matrizR[i][j]=0
                    matrizAi_mas_1[i][j]=0

    #los eigenvalores son los elementos de la diagonal
    for i in range(0,n):
        print(str(matrizAi_mas_1[i][i]), ", ")

#Metodo para sacar KI
def sacaKI(vi, qj, n):
    #declaracion de variables
    ki=0
    numerador=denominador=0

    #obtiene el numerador
    for i in range (0, n):
        numerador=numerador + (vi[i]*qj[i])

    #obtiene el denominador
    for i in range (0, n):
        denominador=denominador + (qj[i]*qj[i])
    
    ki=numerador/denominador
    return ki

#Metodo para algoritmo de GramSchmidt
def gram(matriz, n, matrizOrtogonal):
    #declaracion de variables
    ki=0
    vi=[]
    qi=[]
    qj=[]
    parentesis=[]
    for i in range(0,n):
        parentesis.append(0)
        vi.append(0)
        qi.append(0)
        qj.append(0)
    contqi=0
    contqj=1
    
    
    #Paso 1: q1=v1
    for i in range(0,n):
         matrizOrtogonal[i][0] = matriz[i][0]
    contqi=1

    #paso 2:
    while(contqi<n):  #while que se encarga de llenar la matriz ortogonal

        #obtener las variables que necesitamos apara aplicar la formula del PDF 29
        #vi
        for i in range (0, n):
            vi[i] = matriz[i][contqi]
        
        #qj
        while(contqj<=contqi): #while que se encarga de sacar a qj #inicializa contqj de nuevo
            for i in range (0, n):
                qj[i] = matrizOrtogonal[i][contqj-1]
            contqj=contqj+1

            #obtener el valor de la fraccion <vi,qj>/<qj,qj> = ki
            ki=sacaKI(vi, qj, n)      
           
            for i in range (0, n):
                parentesis[i] = parentesis[i] - (ki*qj[i])
        
        #se suma vi+parentesis (no es resta porque parentesis ya es negativo)
        for i in range (0, n):
            qi[i] = vi[i]+parentesis[i]

        #agregar qi a la matriz ortogonal
        for i in range(0,n):
            matrizOrtogonal[i][contqi] = qi[i]
        
        #reestablacer valores
        contqi=contqi+1
        contqj=1
        for i in range(0,n):
            parentesis[i]=0


#Proceso de ortonormalizacion de la matriz ortogonal
def ortonormalizar(matrizOrtogonal, n, matrizOrtonormal): 
    #declaracion de variables
    cont=0
    ui=[]
    wi=[]
    wiNorma=0
    for i in range(0,columnas):
        ui.append(0)
        wi.append(0)
    
    for i in range(0,n):
        #obtencion de variables para hacer ui=wi/wiNorma
        #wi
        for b in range (0,n):
            wi[b]=matrizOrtogonal[b][cont]

        #wiNorma
        for b in range (0,n):
            wiNorma=wiNorma+(wi[b]*wi[b])
        wiNorma=sqrt(wiNorma)

        #ui
        for b in range (0,n):
            ui[b]=(wi[b])/wiNorma

        #agregar ui a la matriz ortonormal
        for b in range(0,n):
            matrizOrtonormal[b][cont] = ui[b]

        cont=cont+1
        wiNorma=0
    
#imprimir la matriz
def imprimirMatriz(matriz, n):
    for i in range(0,n):
        for j in range(0,n):
            #imprime el numero redondead0 por 5 cifras
            print(" | ",round(matriz[i][j], 5), end=" | ") 
        print()

#sacar transpuesta
def hacerTranspuesta(matrizTranspuesta, n, matrizBase):
    for i in range(0,n):
        for j in range(0,n):
            matrizTranspuesta[i][j] = matrizBase[j][i]
 
#multiplica dos matrices
def multiplicacionMatrices(matrizIzquierda, matrizDerecha, matrizResultado, n):
    for i in range(0,n):
        for j in range(0,n):
            for k in range(0,n):
                matrizResultado[i][j] += matrizIzquierda[i][k]*matrizDerecha[k][j]

                
#main-------------------------------------------------------------------------------------

#Llamamos a todoso los metodos para poder resolver la matriz
#declaracion de variables
filas=columnas=0

print ("Se creará una matriz de nxn")
print("El numero de vectores será igual al numero n del conjunto R")
n = int(input("Ingrese el valor de n: "))
filas = columnas = n

matrizRandom=[[4,6],[0,8]]

matriz=[]
matrizOrtonormal =[]
matrizOrtogonal =[]
matrizTranspuesta=[]
matrizResultado=[]
for i in range(filas):
    matriz.append([0]*columnas)
    matrizOrtogonal.append([0]*columnas)
    matrizOrtonormal.append([0]*columnas)
    matrizTranspuesta.append([0]*columnas)
    matrizResultado.append([0]*columnas)
    
#rellenar la matriz con vectores acomodados por columnas, no por filas
for j in range(columnas):
    for i in range(filas):
        matriz[i][j] = float(input("Ingresa la componente del vector "+str(j)+" posicion "+str(i)+": "))

# Imprimir matriz original
print("Matriz Original: -----------------------------------------------------------------------------------")
imprimirMatriz(matriz,n)

print ("\n Nueva base ortogonal con GramSchmidt: ----------------------------------------------------------")
gram(matriz, n, matrizOrtogonal)
imprimirMatriz(matrizOrtogonal,n)

print("\n Nueva base ortonormal: --------------------------------------------------------------------------")
ortonormalizar(matrizOrtogonal, n, matrizOrtonormal)
imprimirMatriz(matrizOrtonormal, n)

print("\n Eigenvalores: -----------------------------------------------------------------------------------")
algoritmoQR(matriz,matrizOrtogonal,matrizOrtonormal,matrizTranspuesta,n)

