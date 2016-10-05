# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 13:44:01 2016
@author:  tfy12dol
"""
from  scipy import *
from  pylab import *
import  HardCodedProblem
import scipy.linalg as lin
#from mpi4py import MPI

n = 20
omega = 1
nbrOfIterations = 1000
tolerance = 1e-5

nx = n # Nbr of points
ny = n
ny2 = 2*n -1
Deltax = 1 / (n -1)

def generateOperator(edof,boundary,conditions,nbrOfNodes):
    Operator = zeros((nbrOfNodes,nbrOfNodes))
#    Deltax = 1 / (nbrOfNodes -1)    
    
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
    
    return Operator / Deltax**2
    
def generateB(edof,boundary,conditions,nbrOfNodes):
    b = zeros((nbrOfNodes,))
#    Deltax = 1 / (nbrOfNodes -1)
    newconditions = conditions.copy()
    
    dirichletBoundaries = []
    for key in conditions.keys():
        if (conditions[key][0] == 'neumann'):
            neumannPoints = array(boundary[key], dtype='int')
#            print(neumannPoints)
#            print(conditions[key][1])
            b[neumannPoints] += -conditions[key][1] / Deltax
        elif (conditions[key][0] == 'dirichlet'):
            dirichletBoundaries += list(boundary[key])
            
            for point in boundary[key]:
                neighbourPoints = edof[int(point),:]
                neighbourPoints = delete(neighbourPoints, where(neighbourPoints < 0))
                neighbourPoints = array(neighbourPoints, dtype='int')
#                indexRemove = []                
#                for i in range(0,len(neighbourPoints)):
#                    if (neighbourPoints[i] < 0):
#                        indexRemove.append(i)
#                neighbourPoints = delete(neighbourPoints, indexRemove)
#                print(neighbourPoints)
#                print(key)
                if (len(conditions[key][1]) > 1):
                    index = where(boundary[key] == point)[0]
#                    print(index)
                    b[neighbourPoints] += -conditions[key][1][index] / Deltax**2 
                else:
                    b[neighbourPoints] += -conditions[key][1] / Deltax**2 
                
    dirichletBoundaries = unique(array(dirichletBoundaries))
    dirichletBoundaries = unique(dirichletBoundaries)       
    
    b = delete(b, (dirichletBoundaries))
    
    return b,newconditions
    
#    for i in range(0,nbrOfNodes):
#        for point in edof[i,:]:
#            for key in boundary.keys():
#                if any(point == boundary[key]):

#def generateSolutionStructure(u,boundary,conditions,nodeStructure):
    

edof2, boundary2, conditions2, nbrOfNodes2,base2,height2 = domainTwo()
Operator2 = generateOperator(edof2, boundary2, conditions2, nbrOfNodes2)
b2,_ = generateB(edof2, boundary2, conditions2, nbrOfNodes2)
u2 = solve(Operator2, b2)
u2 = u2.reshape((ny2)-2,(nx)-2)
u2old = u2

domain1flow = (u2[len(boundary2['leftUpper'])-1:,0] - conditions2['leftLower'][1]) / Deltax # The value of the flow into the other domains. Positive --> other domain is heated
domain3flow = (u2[0:len(boundary2['rightUpper']),-1] - conditions2['rightUpper'][1]) / Deltax

edof1, boundary1, conditions1, nbrOfNodes1,base1,height1 = domainOne()
Operator1 = generateOperator(edof1, boundary1, conditions1, nbrOfNodes1)
conditions1['rightUpper'] = (conditions1['rightUpper'][0], domain1flow)
b1,_ = generateB(edof1, boundary1, conditions1, nbrOfNodes1)
u1 = solve(Operator1,b1)
u1 = u1.reshape((nx)-2,(nx)-1)
u1old = u1

edof3, boundary3, conditions3, nbrOfNodes3,base3,height3 = domainThree()
Operator3 = generateOperator(edof3, boundary3, conditions3, nbrOfNodes3)
conditions3['leftUpper'] = (conditions3['leftUpper'][0], domain3flow)
b3,_ = generateB(edof3, boundary3, conditions3, nbrOfNodes3)
u3 = solve(Operator3, b3)
u3 = u3.reshape((nx)-2,(nx)-1) 
u3old = u3
#u1 = 0
#u3 = 0

for i in range(0,nbrOfIterations):
    print('Iteration {} of {}...'.format(i+1,nbrOfIterations))
    u1old = u1
    conditions1old = conditions1
    edof1, boundary1, conditions1, nbrOfNodes1,base1,height1 = domainOne()
    Operator1 = generateOperator(edof1, boundary1, conditions1, nbrOfNodes1)
    conditions1['rightUpper'] = (conditions1['rightUpper'][0], domain1flow)
    b1,_ = generateB(edof1, boundary1, conditions1, nbrOfNodes1)
    
#    lu1 = lin.lu_factor(Operator1)
#    u1 = lin.lu_solve(lu1, b1)
    u1 = solve(Operator1,b1)
    u1 = u1.reshape((nx)-2,(nx)-1)
    u1 = omega * u1 + (1 - omega) * u1old
    
    
    u3old = u3
    conditions3old = conditions3
    edof3, boundary3, conditions3, nbrOfNodes3,base3,height3 = domainThree()
    Operator3 = generateOperator(edof3, boundary3, conditions3, nbrOfNodes3)
    conditions3['leftUpper'] = (conditions3['leftUpper'][0], domain3flow)
    b3,_ = generateB(edof3, boundary3, conditions3, nbrOfNodes3)
    
#    lu3 = lin.lu_factor(Operator3)
#    u3 = lin.lu_solve(lu3, b3)
    u3 = solve(Operator3, b3)
    u3 = u3.reshape((nx)-2,(nx)-1) 
    u3 = omega * u3 + (1 - omega) * u3old
    
    
    u2old = u2
    conditions2old = conditions2
    edof2, boundary2, conditions2, nbrOfNodes2,base2,height2 = domainTwo()
    Operator2 = generateOperator(edof2, boundary2, conditions2, nbrOfNodes2)
    conditions2['leftLower'] = (conditions2['leftLower'][0], u1[:,-1])
    conditions2['rightUpper'] = (conditions2['rightUpper'][0], u3[:,0])
    b2,_ = generateB(edof2, boundary2, conditions2, nbrOfNodes2)
    
#    lu2 = lin.lu_factor(Operator2)
#    u2 = lin.lu_solve(lu2, b2)    
    u2 = solve(Operator2, b2)
    u2 = u2.reshape((ny2)-2,(nx)-2)
#    matshow(u2)
    u2 = omega * u2 + (1 - omega) * u2old
    
    normDiffOld = normDiff 
    normDiff = norm(u2 - u2old) 
    print('normDiff: ', normDiff)
    print('omega: ', omega)
    omegaOld = omega
##    if (normDiff - normDiffOld) > 0:
##        omega = 0.8 * omega / normDiff
    
       
    omega = omega / (normDiff * (n/30))
    if omega > normDiff  / (n/30):
        omega = normDiff  / (n/30)
#    omega = normDiff
#    if (normDiff > normDiffOld):    
#        omega = omegaOld * normDiffOld / 0.8
#    if (omega > normDiff):
#        omega = normDiff
        
    if (normDiff > 1e2):
#        print(normDiff)
#        print(normDiffOld)
        print('Solution might not be converging.')
        omega = 0.9 * omegaOld
        #        print('Setting omega to: ', omega)
        u3 = u3old
        u2 = u2old
        u1 = u1old

        normDiff = normDiffOld        
        
        conditions1 = conditions1old
        conditions2 = conditions2old
        conditions3 = conditions3old
    elif (normDiff < tolerance):
        print('Solution appears to have converged within tolerance: ',tolerance)
        break
        omega = omega * normDiff
#        omega = 0.01 / normDiff        
#        print('Setting omega to: ', omega)
     
    
    domain1flow = (u2[len(boundary2['leftUpper'])-1:,0] - conditions2['leftLower'][1]) / Deltax # The value of the flow into the other domains. Positive --> other domain is heated
    domain3flow = (u2[0:len(boundary2['rightUpper']),-1] - conditions2['rightUpper'][1]) / Deltax

c1 = zeros((nx-1,nx-1)) + nan
c2 = zeros((nx-1,nx-1)) + nan
c1[-1,:] = 15
c1[:,-1] = 15
c2[:,0] = 15
c2[0,:] = 15  
part1 = vstack((c1,u1))
part3 = vstack((u3,c2))
part2 = hstack((part1,u2))
final = hstack((part2,part3))
final = pad(final,((1,1),(1,1)),'constant',constant_values=(nan, nan))
#final[:,0] = nan
final[nx:-1,0] = 40
#final[0,:] = nan
#final[ny2-1,:] = nan
final[-1,1:nx] = 15
final[-1,nx:2*nx-2] = 5
final[1:nx-1,-1] = 40
final[0,nx:ny2-1] = 40
final[0,ny2-1:-1] = 15
fig=matshow(final)



#domain1boundary = array(boundary2['leftLower'], dtype='int')
#domain1boundary = array(edof2[domain1boundary,3], dtype='int')
#domain1flow = (conditions2['leftLower'][1] - u2[domain1boundary]) / Deltax
#
#domain2boundary = array(boundary2['rightUpper'], dtype='int')
#domain2boundary = array(edof2[domain2boundary,1], dtype='int')
#domain2flow = (u2[domain1boundary] - conditions2['rightUpper'][1]) / Deltax









