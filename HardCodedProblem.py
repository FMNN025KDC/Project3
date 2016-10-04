# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 13:08:13 2016
@author:  tfy12dol
"""
from  scipy import *
from  pylab import *

n = 5 # Nbr of points
Deltax = 1 / (n -1)
h1 = 1
b1 = 1
h3 = 1
b3 = 1
h2 = 2
h2l1 = h2 - h1 
h2l2 = h1
h2r1 = h3
h2r2 = h2 - h3
b2 = 1


def generateGrid(b,h, DomainTwo = False):
    nx = n * b
    ny = n * h
#    Dx = self.b / (nx - 1)
    
    pointNames = zeros((ny,nx))
    name = 0
    for i in range(0,ny):
        for j in range(0,nx):
            pointNames[i,j] = name
            name += 1
            
    if DomainTwo:
        boundary = {'leftUpper': pointNames[0:(n*h2l1)+1,0], 'leftLower': pointNames[(n*h2l1)+1:,0], 'rightUpper': pointNames[0:(n*h2r1),-1], 'rightLower': pointNames[(n*h2r1):,-1], 'upper': pointNames[0,:], 'lower': pointNames[-1,:]}
    else:
        boundary = {'leftUpper': pointNames[:,0], 'leftLower': array([]), 'rightUpper': pointNames[:,-1], 'rightLower': array([]), 'upper': pointNames[0,:], 'lower': pointNames[-1,:]}
    
    pointNamespadded = pad(pointNames,((1,1),(1,1)),'constant',constant_values=(-1, -1))
    
    edof = zeros((size(pointNames),5)) -1
    for i in range(1,shape(pointNamespadded)[0]-1):
        for j in range(1,shape(pointNamespadded)[1]-1):
            edofij = [pointNamespadded[i-1,j],pointNamespadded[i,j-1],pointNamespadded[i,j],pointNamespadded[i,j+1],pointNamespadded[i+1,j]]
            edof[int(pointNamespadded[i,j])] = edofij
            
    return edof,boundary


def domainOne():
    h = h1
    b = b1
    
    conditions = {'leftUpper': ('dirichlet',40), 'leftLower': ('dirichlet',40), 'rightUpper': ('neumann',None), 'rightLower': ('neumann',None), 'upper': ('dirichlet',15), 'lower': ('dirichlet',15)}
    edof,boundary = generateGrid(b,h)
    
    return edof,boundary,conditions,n**2*b*h
    
def domainThree():
    h = h3
    b = b3
    
    conditions = {'leftUpper': ('neumann',None), 'leftLower': ('neumann',None), 'rightUpper': ('dirichlet',40), 'rightLower': ('dirichlet',40), 'upper': ('dirichlet',15), 'lower': ('dirichlet',15)}
    edof,boundary = generateGrid(b,h)
    
    return edof,boundary,conditions,n**2*b*h
    
def domainTwo():
    h = h2
    b = b2
    
    conditions = {'leftUpper': ('dirichlet',15), 'leftLower': ('dirichlet',15), 'rightUpper': ('dirichlet',15), 'rightLower': ('dirichlet',15), 'upper': ('dirichlet',40), 'lower': ('dirichlet',5)}
    edof,boundary = generateGrid(b,h,True)
    
    return edof,boundary,conditions,n**2*b*h
    
    
    
    
    
        
    