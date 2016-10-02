# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 11:06:48 2016
@author: tfy12dol
"""
from  scipy import *
from  pylab import *

A = array([[1.,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15],[16,17,18,19,20]]) - 1
A = array([[0,1,2],[3,4,5],[6,7,8]])

Apadded = pad(A,((1,1),(1,1)),'constant',constant_values=(-1, -1))
#Apadded = array([[-1,-1,-1,-1],list(A.T),[-1,-1,-1,-1]])
#Apadded = array([[-1,-1,-1,-1,-1,-1,-1],Apadded.T,[-1,-1,-1,-1,-1,-1,-1]])

edof = zeros((size(A),5)) -1
for i in range(1,shape(Apadded)[0]-1):
    for j in range(1,shape(Apadded)[1]-1):
        edofij = [Apadded[i-1,j],Apadded[i,j-1],Apadded[i,j],Apadded[i,j+1],Apadded[i+1,j]]
#        edof[(i-1)*shape(A)[1]+j-1] = edofij
        edof[Apadded[i,j]] = edofij
        
Operator = zeros((size(A),size(A)))
i = 0
for i in range(0,shape(edof)[0]):
    for point in edof[i,:]:
        if point == -1:
            print('Nope')
        elif point == i:
            Operator[i,int(point)] = -4
        else:
            Operator[i,int(point)] = 1