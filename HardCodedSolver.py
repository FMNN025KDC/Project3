# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 13:44:01 2016
@author:  tfy12dol
"""
from  scipy import *
from  pylab import *
import  HardCodedProblem


def generateOperator(edof,boundary,conditions,nbrOfNodes):
    Operator = zeros((nbrOfNodes,nbrOfNodes))
    for i in range(0,nbrOfNodes):
#        if any(self.edof[i,2] == self.boundary['left']) or any(self.edof[i,2] == self.boundary['right']) or any(self.edof[i,2] == self.boundary['upper']) or any(self.edof[i,2] == self.boundary['lower']):
#            pass
#        else:
    
        if any(edof[i,:] == -1):
            for key in boundary.keys():
                if any(edof[i,2] == boundary[key]):
                    condition = conditions[key]
                    if (condition[0] == 'dirichlet'):
                        break
        else:
            condition = (None, None)
                    
        for point in edof[i,:]:
#            print(condition[0])
            if (condition[0] == 'neumann'):
                if point == -1: 
                    pass
                elif point == i:
                    Operator[i,int(point)] = -3
                else:
                    Operator[i,int(point)] = 1
            else:   
#                print('här')                     
                if point == -1:
                    pass
                elif point == i:
#                    print('där')
                    Operator[i,int(point)] = -4
                else:
                    Operator[i,int(point)] = 1
#            print(Operator)
#            input()
    
#    points = array(list(self.boundary['left']) + list(self.boundary['right']) + list(self.boundary['upper']) + list(self.boundary['lower']))
#        points = array([self.boundary['left'],self.boundary['right'],self.boundary['upper'],self.boundary['lower']])
#    points = unique(points)
    dirichletBoundaries = []
    for key in conditions:
        if (conditions[key][0] == 'dirichlet'):
            dirichletBoundaries += list(boundary[key])
            
    dirichletBoundaries = unique(array(dirichletBoundaries))
#    print(dirichletBoundaries)
    Operator = delete(Operator, (dirichletBoundaries), axis=0)
    Operator = delete(Operator, (dirichletBoundaries), axis=1)        
    
    return Operator
    
def generateB(edof,boundary,conditions,nbrOfNodes):
    b = zeros((nbrOfNodes,))
    Deltax = 1 / (nbrOfNodes -1)
    
    for key in conditions.keys():
        if (conditions[key][0] == 'neumann'):
            b(boundary[key]) = -conditions[key][1] / Del
        elif (conditions[key][0] == 'dirichlet'):
            
    
#    for i in range(0,nbrOfNodes):
#        for point in edof[i,:]:
#            for key in boundary.keys():
#                if any(point == boundary[key]):


    
#edof1, boundary1, conditions1, nbrOfNodes1 = domainOne()
#Operator1 = generateOperator(edof1, boundary1, conditions1, nbrOfNodes1)
#edof2, boundary2, conditions2, nbrOfNodes2 = domainTwo()
#Operator2 = generateOperator(edof2, boundary2, conditions2, nbrOfNodes2)
edof3, boundary3, conditions3, nbrOfNodes3 = domainThree()
Operator3 = generateOperator(edof3, boundary3, conditions3, nbrOfNodes3)










