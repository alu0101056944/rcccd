#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Robótica Computacional - Curso 2014/2015
# Grado en Ingeniería Informática (Cuarto)
# Práctica: Resolución de la cinemática inversa mediante CCD
#           (Cyclic Coordinate Descent).

import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import colorsys as cs

# ******************************************************************************
# Declaración de funciones

def muestra_origenes(O,final=0):
  # Muestra los orígenes de coordenadas para cada articulación
  print('Origenes de coordenadas:')
  for i in range(len(O)):
    print('(O'+str(i)+')0\t= '+str([round(j,3) for j in O[i]]))
  if final:
    print('E.Final = '+str([round(j,3) for j in final]))

def muestra_robot(O,obj):
  # Muestra el robot graficamente
  plt.figure(1)
  plt.xlim(-L,L)
  plt.ylim(-L,L)
  T = [np.array(o).T.tolist() for o in O]
  for i in range(len(T)):
    plt.plot(T[i][0], T[i][1], '-o', color=cs.hsv_to_rgb(i/float(len(T)),1,1))
  plt.plot(obj[0], obj[1], '*')
  plt.show()
  raw_input()
  plt.clf()

def matriz_T(d,th,a,al):
  # Calcula la matriz T (ángulos de entrada en grados)
  
  return [[cos(th), -sin(th)*cos(al),  sin(th)*sin(al), a*cos(th)]
         ,[sin(th),  cos(th)*cos(al), -sin(al)*cos(th), a*sin(th)]
         ,[      0,          sin(al),          cos(al),         d]
         ,[      0,                0,                0,         1]
         ]

def cin_dir(th,a):
  #Sea 'th' el vector de thetas
  #Sea 'a'  el vector de longitudes
  T = np.identity(4)
  o = [[0,0]]
  for i in range(len(th)):
    T = np.dot(T,matriz_T(0,th[i],a[i],0))
    tmp=np.dot(T,[0,0,0,1])
    o.append([tmp[0],tmp[1]])
  return o

# ******************************************************************************
# Cálculo de la cinemática inversa de forma iterativa por el método CCD

# valores articulares arbitrarios para la cinemática directa inicial
th=[0.,0.,0.]   # ángulos de cada punto
a =[5.,5.,5.]   # distancias de cada parte
L = sum(a) # variable para representación gráfica
EPSILON = .01

plt.ion() # modo interactivo

# introducción del punto para la cinemática inversa
if len(sys.argv) != 3:
  sys.exit("python " + sys.argv[0] + " x y")
objetivo=[float(i) for i in sys.argv[1:]]

# O=range(len(th)+1) # Reservamos estructura en memoria
O=range(len(th)+1)
O[0]=cin_dir(th,a) # Calculamos la posicion inicial
print "- Posicion inicial:"
muestra_origenes(O[0])

dist = float("inf")
prev = 0.
iteracion = 1
print "Objetivo", objetivo
while (dist > EPSILON and abs(prev-dist) > EPSILON/100.):
  prev = dist
  indexOfR = len(th) # length = 3
  # Para cada combinación de articulaciones:
  # Por cada iteración, cambiamos las coordenada de todos 
  # los puntos de la articulación:
  for i in range(len(th)): # En este caso, itera desde 0 a 2
    indexOfR -= 1 # previous of last
    print "Iteracion: ", iteracion, "------------------------------------------------"
    print "Posiciones", O
    # print "IndexOfR:", indexOfR
    # cálculo de la cinemática inversa:
    # Cálculo del punto E:
    posXOfR = O[i][indexOfR][0]
    posYOfR = O[i][indexOfR][1]
    # print "posXOfR: ", posXOfR
    # print "posYOfR: ", posYOfR
    # calcular distancias de R a objetivo
    vectRToObj = [objetivo[0] - posXOfR, objetivo[1] - posYOfR]
    distRToObj = sqrt(vectRToObj[0] ** 2 + vectRToObj[1] ** 2)
    # distancia de R a O_último
    vectRToLast = [O[i][3][0] - posXOfR, O[i][3][1] - posYOfR]
    distRToLast = sqrt(vectRToLast[0] ** 2 + vectRToLast[1] ** 2)
    # y distancia de E a objetivo
    # vectLastToObj = [O[i][3][0] - objetivo[0], O[i][3][1] - objetivo[1]]
    # distLastToObj= sqrt(vectLastToObj[0] ** 2 + vectLastToObj[1] ** 2)
    # Producto vectorial
    dividendoVectorial = (vectRToObj[0] * vectRToLast[0]) + (vectRToObj[1] * vectRToLast[1])
    divisorVectorial = distRToLast * distRToObj
    # Imprimimos datos de distancias
    # print "vectRToObj: ", vectRToObj    # (0, 10)
    # print "DistRToObj: ", distRToObj    # 10
    # print "vectRToLast: ", vectRToLast  # ()
    # print "DistRToLast: ", distRToLast  # 5
    # Calcular el ángulo entre E y Objetivo
    print "dividendoVectorial", dividendoVectorial
    print "divisorVectorial: ", divisorVectorial
    thetaLastToObj = np.arccos(dividendoVectorial / divisorVectorial)
    print "thetaLastToObj: ", thetaLastToObj
    # Calcular vector perpendicular por producto vectorial
    # productoVectorial = distRToObj * distRToLast * sin(thetaLastToObj)
    # print "ProductoVectorial: ", productoVectorial
    productoVectorial = np.cross(O[i][3], objetivo)
    if productoVectorial < 0:
      thetaLastToObj = -thetaLastToObj
    # Añadir el ángulo a R
    indexFromLast = len(th) - i - 1
    # Calcula la cinemática directa de la tabla construida, calcula el conjunto de nuevos puntos y los asigna
    th[indexFromLast] += thetaLastToObj
    O[i+1] = cin_dir(th,a)

  print "O[-1]:", O[-1]
  print "O[-1][-1]:", O[-1][-1], "A"
  dist = np.linalg.norm(np.subtract(objetivo,O[-1][-1]))
  print "\n- Iteracion " + str(iteracion) + ':'
  muestra_origenes(O[-1])
  muestra_robot(O,objetivo)
  print "Distancia al objetivo = " + str(round(dist,5))
  iteracion+=1
  O[0]=O[-1]

if dist <= EPSILON:
  print "\n" + str(iteracion) + " iteraciones para converger."
else:
  print "\nNo hay convergencia tras " + str(iteracion) + " iteraciones."
print "- Umbral de convergencia epsilon: " + str(EPSILON)
print "- Distancia al objetivo:          " + str(round(dist,5))
print "- Valores finales de las articulaciones:"
for i in range(len(th)):
  print "  theta" + str(i+1) + " = " + str(round(th[i],3))
for i in range(len(th)):
  print "  L" + str(i+1) + "     = " + str(round(a[i],3))

""" def getAngle(punto1, center, punto2):
  centerToPunto2 = atan2(punto2[1] - center[1], punto2[0]-center[0])
  centerToPunto1 = atan2(punto1[1] - center[1], punto1[0]-center[0])
  return centerToPunto2 - centerToPunto1 """
  
