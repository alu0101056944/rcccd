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

# Modificación:
# Añadir ángulos máximos para cada una de las articulaciones en ambos sentidos
max_th=[3.1416, 3.1416, 3.1416]
# max_th=[3.1416, 3.1416, 3.1416]  # Valores en radianes por defecto. Es decir, permiten flexibilidad plena.

# valores articulares arbitrarios para la cinemática directa inicial
th=[0.,0.,0.]   # ángulos de cada punto
a =[2., 3., 3.]   # distancias de cada parte
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
  # Por cada iteración, cambiamos las coordenada de todos los puntos de la articulación:
  for i in range(len(th)):
    indexOfR -= 1 # previous of last
    # Cálculo de la cinemática inversa:
    posXOfR = O[i][indexOfR][0]
    posYOfR = O[i][indexOfR][1]
    vectRToObj = [objetivo[0] - posXOfR, objetivo[1] - posYOfR]
    distRToObj = sqrt(vectRToObj[0] ** 2 + vectRToObj[1] ** 2)
    vectRToLast = [O[i][3][0] - posXOfR, O[i][3][1] - posYOfR]
    distRToLast = sqrt(vectRToLast[0] ** 2 + vectRToLast[1] ** 2)
    # Producto vectorial
    dividendoVectorial = (vectRToObj[0] * vectRToLast[0]) + (vectRToObj[1] * vectRToLast[1])
    divisorVectorial = distRToLast * distRToObj
    # Calcular el ángulo entre E y Objetivo
    thetaLastToObj = np.arccos(dividendoVectorial / divisorVectorial)
    # Calcular el producto vectorial
    productoVectorial = np.cross([vectRToLast[0], vectRToLast[1], 0], [vectRToObj[0], vectRToObj[1], 0])
    print "productoVectorial: ", productoVectorial
    if productoVectorial[2] < EPSILON:
      print "negative in arm ", i
      print "productoVectorial[2]: ", productoVectorial[2]
      thetaLastToObj = -thetaLastToObj
    # Calcula la cinemática directa de la tabla construida, calcula el conjunto de nuevos puntos y los asigna
    th[indexOfR] += thetaLastToObj
    # if th[indexOfR] < -max_th[indexOfR]:
    #   th[indexOfR] = -max_th[indexOfR]
    # elif th[indexOfR] > max_th[indexOfR]:
    #   th[indexOfR] = max_th[indexOfR]
    # Antes de calcular la cinemática directa, sabemos la distancia que falta respecto al punto
    if max_th[indexOfR] == 0:
      vectLastToObj = [objetivo[0] - O[i][3][0], objetivo[1] - O[i][3][1]]
      distLastToObj = sqrt(vectLastToObj[0] ** 2 + vectLastToObj[1] ** 2)
      if distLastToObj > EPSILON:
        a[indexOfR] += distLastToObj
        print "a: ", a[indexOfR]
    # Calcular la cinemática directa
    O[i+1] = cin_dir(th,a)

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
  
