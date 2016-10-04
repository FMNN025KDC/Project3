# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 13:08:13 2016
@author:  tfy12dol
"""
from  scipy import *
from  pylab import *

n = 10 # Nbr of points
Deltax = 1 / (n -1)
heigth1 = 1
base1 = 1
heigth3 = 1
base3 = 1
heigth2 = 2
heigth2l1 = heigth2 - heigth1 
heigth2l2 = heigth1
heigth2r1 = heigth3
heigth2r2 = heigth2 - heigth3
base2 = 1


def generateGrid(nx,ny, DomainTwo = False):
#    nx = n * b
#    ny = n * h
#    Dx = self.b / (nx - 1)
#    print((ny,nx))
    pointNames = zeros((ny,nx))
    name = 0
    for i in range(0,ny):
        for j in range(0,nx):
            pointNames[i,j] = name
            name += 1
            
    if DomainTwo:
        boundary = {'leftUpper': pointNames[1:(n*h2l1)+1,0], 'leftLower': pointNames[(n*h2l1)+1:-1,0], 'rightUpper': pointNames[1:(n*h2r1)-1,-1], 'rightLower': pointNames[(n*h2r1)-1:-1,-1], 'upper': pointNames[0,:], 'lower': pointNames[-1,:]}
    else:
        boundary = {'leftUpper': pointNames[1:-1,0], 'leftLower': array([]), 'rightUpper': pointNames[1:-1,-1], 'rightLower': array([]), 'upper': pointNames[0,:], 'lower': pointNames[-1,:]}
    
    pointNamespadded = pad(pointNames,((1,1),(1,1)),'constant',constant_values=(-1, -1))
    
    edof = zeros((size(pointNames),5)) -1
    for i in range(1,shape(pointNamespadded)[0]-1):
        for j in range(1,shape(pointNamespadded)[1]-1):
            edofij = [pointNamespadded[i-1,j],pointNamespadded[i,j-1],pointNamespadded[i,j],pointNamespadded[i,j+1],pointNamespadded[i+1,j]]
            edof[int(pointNamespadded[i,j])] = edofij
            
    return edof,boundary


def domainOne():
    h = heigth1
    b = base1
    nx = 10
    ny = 10    
    
    conditions = {'leftUpper': ('dirichlet',array([40])), 'leftLower': ('dirichlet',array([40])), 'rightUpper': ('neumann',array([1])), 'rightLower': (None,None), 'upper': ('dirichlet',array([15])), 'lower': ('dirichlet',array([15]))}
    edof,boundary = generateGrid(nx,ny)    
    
    return edof,boundary,conditions,nx*ny,b,h
    
def domainThree():
    h = heigth3
    b = base3
    nx = 10
    ny = 10  

    conditions = {'leftUpper': ('neumann',array([1])), 'leftLower': (None,None), 'rightUpper': ('dirichlet',array([40])), 'rightLower': ('dirichlet',array([40])), 'upper': ('dirichlet',array([15])), 'lower': ('dirichlet',array([15]))}
    edof,boundary = generateGrid(nx,ny)
    
    return edof,boundary,conditions,nx*ny,b,h
    
def domainTwo():
    h = heigth2
    b = base2
    nx = 10
    ny = 19    
    
    conditions = {'leftUpper': ('dirichlet',array([15])), 'leftLower': ('dirichlet',array([15])), 'rightUpper': ('dirichlet',array([15])), 'rightLower': ('dirichlet',array([15])), 'upper': ('dirichlet',array([40])), 'lower': ('dirichlet',array([5]))}
    edof,boundary = generateGrid(nx,ny,True)
    
    return edof,boundary,conditions,nx*ny,b,h
    
    
    
    
    
        
    