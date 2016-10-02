# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 12:22:12 2016
@author: tfy12dol
"""
from  scipy import *
from  pylab import *


class Room():
    def __init__(self, pointDensity = 10, b = 1,h = 1, leftCond = ("dirichlet",40), rightCond = ("dirichlet",15), upperCond = ("dirichlet",15), lowerCond = ("dirichlet",5)):
        self.conditions = {'leftCond': leftCond, 'rightCond': rightCond, 'upperCond': upperCond, 'lowerCond': lowerCond}
        self.b = b
        self.h = h
        self.pointDensity = pointDensity
        self.edof,self.Deltax,self.boundary = self.generateGrid() 
        
    def generateGrid(self):
        nx = self.pointDensity * self.b
        ny = self.pointDensity * self.h
        Dx = self.b / (nx - 1)
        
        pointNames = zeros((ny,nx))
        name = 0
        for i in range(0,ny):
            for j in range(0,nx):
                pointNames[i,j] = name
                name += 1
                
        
        boundary = {'left': pointNames[:,0], 'right': pointNames[:,-1], 'upper': pointNames[0,:], 'lower': pointNames[-1,:]}
        pointNamespadded = pad(pointNames,((1,1),(1,1)),'constant',constant_values=(-1, -1))
        
        edof = zeros((size(pointNames),5)) -1
        for i in range(1,shape(pointNamespadded)[0]-1):
            for j in range(1,shape(pointNamespadded)[1]-1):
                edofij = [pointNamespadded[i-1,j],pointNamespadded[i,j-1],pointNamespadded[i,j],pointNamespadded[i,j+1],pointNamespadded[i+1,j]]
                edof[int(pointNamespadded[i,j])] = edofij
                
        return edof,Dx,boundary