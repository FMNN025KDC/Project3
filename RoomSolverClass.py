# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 12:55:28 2016
@author: tfy12dol
"""
from  scipy import *
from  pylab import *
import matplotlib


class RoomSolver():
    def __init__(self, room):
        self.Deltax = room.Deltax
        self.edof = room.edof
        self.conditions = room.conditions
        self.boundary = room.boundary
        
    def generateOperator(self):
        Operator = zeros((len(self.boundary['left'])*len(self.boundary['lower']),len(self.boundary['left'])*len(self.boundary['lower'])))
        for i in range(0,shape(self.edof)[0]):
            if any(self.edof[i,2] == self.boundary['left']) or any(self.edof[i,2] == self.boundary['right']) or any(self.edof[i,2] == self.boundary['upper']) or any(self.edof[i,2] == self.boundary['lower']):
                pass
            else:
                for point in self.edof[i,:]:
                    if point == -1:
                        pass
                    elif point == i:
                        Operator[i,int(point)] = -4
                    else:
                        Operator[i,int(point)] = 1
        
        points = array(list(self.boundary['left']) + list(self.boundary['right']) + list(self.boundary['upper']) + list(self.boundary['lower']))
#        points = array([self.boundary['left'],self.boundary['right'],self.boundary['upper'],self.boundary['lower']])
        points = unique(points)
        
        Operator = delete(Operator, (points), axis=0)
        Operator = delete(Operator, (points), axis=1)        
        
        return Operator
        
    def generateB(self):
        b = zeros(((len(self.boundary['left'])*len(self.boundary['lower']),)))

        for i in range(0,shape(self.edof)[0]):
            for point in self.edof[i,:]:
                if any(self.boundary['left'] == point):
                    b[i] += self.conditions['leftCond'][1]
                    
                if any(self.boundary['right'] == point):
                    b[i] += self.conditions['rightCond'][1]
                    
                if any(self.boundary['upper'] == point):
                    b[i] += self.conditions['upperCond'][1]

                if any(self.boundary['lower'] == point):
                    b[i] += self.conditions['lowerCond'][1]  

        points = array(list(self.boundary['left']) + list(self.boundary['right']) + list(self.boundary['upper']) + list(self.boundary['lower']))
#        points = array([self.boundary['left'],self.boundary['right'],self.boundary['upper'],self.boundary['lower']])
        points = unique(points)       
        
        b = delete(b, (points))        
        
        return -b
#    def solveStationary(self):
        
            
        
        
        
        
        
        
        