# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 13:44:01 2016
@author:  tfy12dol
"""
from  scipy import *
from  pylab import *
import  HeatEquationProblem
import scipy.linalg as lin
from mpi4py import MPI
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve





# Point density (per unit length). Has to be changed in HeatEquationProblem
n = HeatEquationProblem.n
# Initial omega. Is modified by the (somewhat unstable) adaptive solver.
omega = 1
# Maximum number of iterations.
nbrOfIterations = 2000

# Tolerance of solution error, compared with previous step. Breaks loop when 
# fulfilled.
tolerance = 1e-12
#initialize the MPI enviroment
comm=MPI.COMM_WORLD
rank=comm.Get_rank()
# Points in x- and y-direction.
nx = n
ny = n
# Points in y-direction for room 2. Has to match the grid of smaller rooms.
ny2 = 2*n -1
Deltax = 1 / (n -1)

def generateOperator(edof,boundary,conditions,nbrOfNodes):
    """Generates an operator matrix based on edof and boundary conditions.
    """
    Operator = zeros((nbrOfNodes,nbrOfNodes))
#    Operator=lil_matrix(Operator)
    
    # Uses edof matrix to generate operator. Works for arbitrary node naming 
    # convention.
    for i in range(0,nbrOfNodes):
        # If connected to 'outside' node (-1) find condition. If both neumann 
        # and dirichlet choose dirichlet.
        if any(edof[i,:] == -1):
            for key in boundary.keys():
                if any(edof[i,2] == boundary[key]):
                    condition = conditions[key]
                    if (condition[0] == 'dirichlet'):
                        break
        else:
            condition = (None, None)
         
        # Go through edof, generate operator.           
        for point in edof[i,:]:
            # For neumann points.
            if (condition[0] == 'neumann'):
                if point == -1: 
                    pass
                elif point == i:
                    Operator[i,int(point)] = -3
                else:
                    Operator[i,int(point)] = 1
            else:   
                # For non-neumann points.    
                if point == -1:
                    pass
                elif point == i:
                    Operator[i,int(point)] = -4
                else:
                    Operator[i,int(point)] = 1
    
    # Do not solve for dirichlet boundary point --> Remove from operator.
    dirichletBoundaries = []
    for key in conditions:
        if (conditions[key][0] == 'dirichlet'):
            dirichletBoundaries += list(boundary[key])
            
    dirichletBoundaries = unique(array(dirichletBoundaries))

    Operator = delete(Operator, (dirichletBoundaries), axis=0)
    Operator = delete(Operator, (dirichletBoundaries), axis=1) 
    #Make operator matrix sparse, timed to be faster despite solution being dense
    Operator=csc_matrix(Operator)
    
    return Operator / Deltax**2
    
def generateB(edof,boundary,conditions,nbrOfNodes):
    """Generates vector b based on edof and boundary conditions.
    """
    b = zeros((nbrOfNodes,))
    
    dirichletBoundaries = []
    for key in conditions.keys():
        # If neumann condition, modify b accordingly.
        if (conditions[key][0] == 'neumann'):
            # Force float to int
            neumannPoints = array(boundary[key], dtype='int')
            b[neumannPoints] += -conditions[key][1] / Deltax
            
        elif (conditions[key][0] == 'dirichlet'):
            # Save all dirichlet boundary point so that they can be removed 
            # from calculation.
            dirichletBoundaries += list(boundary[key])
            
            # All points connected to dirichlet point must have their b 
            # modified.
            for point in boundary[key]:
                neighbourPoints = edof[int(point),:]
                # Delete outside (-1) points.
                neighbourPoints = delete(neighbourPoints, where(neighbourPoints < 0))
                neighbourPoints = array(neighbourPoints, dtype='int')

                if (len(conditions[key][1]) > 1):
                    # If one condition per point in boundary, set b value for 
                    # each point.
                    index = where(boundary[key] == point)[0]
                    b[neighbourPoints] += -conditions[key][1][index] / Deltax**2 
                else:
                    # If one condition for all points in boundary, set b value
                    # for all points.
                    b[neighbourPoints] += -conditions[key][1] / Deltax**2 
    
    # Remove dirichlet points from calculation.            
    dirichletBoundaries = unique(array(dirichletBoundaries))
    dirichletBoundaries = unique(dirichletBoundaries)       
    
    b = delete(b, (dirichletBoundaries))

    return b
    
# define tags for MPI
tagNormDiff=0
tagDomain1flow=1
tagDomain3flow=tagDomain1flow
tagu3=3
tagu1=2
   
# Calculate initial u.
if rank == 2:    
    edof2, boundary2, conditions2, nbrOfNodes2,base2,height2 = HeatEquationProblem.domainTwo()
    Operator2 = generateOperator(edof2, boundary2, conditions2, nbrOfNodes2)
    
    b2 = generateB(edof2, boundary2, conditions2, nbrOfNodes2)
    u2 = spsolve(Operator2, b2)
    #u2 = solve(Operator2, b2)
    u2 = u2.reshape((ny2)-2,(nx)-2)
    u2old = u2
    
    domain1flow = (u2[len(boundary2['leftUpper'])-1:,0] - conditions2['leftLower'][1]) / Deltax # The value of the flow into the other domains. Positive --> other domain is heated
    domain3flow = (u2[0:len(boundary2['rightUpper']),-1] - conditions2['rightUpper'][1]) / Deltax
    
    comm.send(domain1flow, 0, tag = tagDomain1flow)
    comm.send(domain3flow, 1, tag = tagDomain3flow)
    
if rank == 0:
    domain1flow=comm.recv(source=2, tag = tagDomain1flow)

    edof1, boundary1, conditions1, nbrOfNodes1,base1,height1 = HeatEquationProblem.domainOne()
    Operator1 = generateOperator(edof1, boundary1, conditions1, nbrOfNodes1)
    conditions1['rightUpper'] = (conditions1['rightUpper'][0], domain1flow)
    b1 = generateB(edof1, boundary1, conditions1, nbrOfNodes1)
    u1 = spsolve(Operator1,b1)
    #u1 = solve(Operator1, b1)
    u1 = u1.reshape((nx)-2,(nx)-1)
    u1old = u1
    comm.send(u1,2,tag = tagu1)

if rank == 1:
    domain3flow=comm.recv(source=2, tag = tagDomain3flow)
    edof3, boundary3, conditions3, nbrOfNodes3,base3,height3 = HeatEquationProblem.domainThree()
    Operator3 = generateOperator(edof3, boundary3, conditions3, nbrOfNodes3)
    conditions3['leftUpper'] = (conditions3['leftUpper'][0], domain3flow)
    b3 = generateB(edof3, boundary3, conditions3, nbrOfNodes3)
    u3 = spsolve(Operator3, b3)
    #u3 = solve(Operator3, b3)
    u3 = u3.reshape((nx)-2,(nx)-1) 
    u3old = u3
    comm.send(u3,2,tag = tagu3)






# Solve with one room per processor.
for i in range(0,nbrOfIterations):

    if rank == 0:
        
        domain1flow=comm.recv(source=2, tag = tagDomain1flow)       
        
        print('Iteration {} of {}...'.format(i+1,nbrOfIterations))  
        u1old = u1
        conditions1old = conditions1
        edof1, boundary1, conditions1, nbrOfNodes1,base1,height1 = HeatEquationProblem.domainOne()
        conditions1['rightUpper'] = (conditions1['rightUpper'][0], domain1flow)
        b1 = generateB(edof1, boundary1, conditions1, nbrOfNodes1)


        u1 = spsolve(Operator1,b1)
#        u1 = solve(Operator1,b1)        
        u1 = u1.reshape((nx)-2,(nx)-1)
        u1 = omega * u1 + (1 - omega) * u1old
        
        comm.send(u1,2,tag = tagu1)
        normDiff=comm.recv(source=2, tag = tagNormDiff)

    
    if rank == 1:
        domain3flow=comm.recv(source=2,tag = tagDomain3flow)
        
        
        u3old = u3
        conditions3old = conditions3
        edof3, boundary3, conditions3, nbrOfNodes3,base3,height3 = HeatEquationProblem.domainThree()
        conditions3['leftUpper'] = (conditions3['leftUpper'][0], domain3flow)
        b3 = generateB(edof3, boundary3, conditions3, nbrOfNodes3)
        
        u3 = spsolve(Operator3, b3)
#        u3 = solve(Operator3, b3)
        u3 = u3.reshape((nx)-2,(nx)-1) 
        u3 = omega * u3 + (1 - omega) * u3old
        
        comm.send(u3,2, tag = tagu3)
        normDiff=comm.recv(source=2, tag = tagNormDiff)

        
    if rank == 2:   
        u1=comm.recv(source=0, tag = tagu1)
        u3=comm.recv(source=1, tag = tagu3)



        u2old = u2
        conditions2old = conditions2
        edof2, boundary2, conditions2, nbrOfNodes2,base2,height2 = HeatEquationProblem.domainTwo()
        conditions2['leftLower'] = (conditions2['leftLower'][0], u1[:,-1])
        conditions2['rightUpper'] = (conditions2['rightUpper'][0], u3[:,0])
        b2 = generateB(edof2, boundary2, conditions2, nbrOfNodes2)

        

        u2 = spsolve(Operator2, b2)
#        u2 = solve(Operator2, b2)
        u2 = u2.reshape((ny2)-2,(nx)-2)
        u2 = omega * u2 + (1 - omega) * u2old     


                
        domain1flow = (u2[len(boundary2['leftUpper'])-1:,0] - conditions2['leftLower'][1]) / Deltax # The value of the flow into the other domains. Positive --> other domain is heated
        domain3flow = (u2[0:len(boundary2['rightUpper']),-1] - conditions2['rightUpper'][1]) / Deltax
        
        comm.send(domain1flow,0, tag = tagDomain1flow)
        comm.send(domain3flow,1, tag = tagDomain3flow)
     
        # Calculate differance between previous and current to see if tolerance
        # is fulfilled.
        normDiff = norm(u2 - u2old) 
        
        comm.send(normDiff,0,tagNormDiff)
        comm.send(normDiff,1,tagNormDiff)

    omegaOld = omega
    
    # Adaptive omega which seems to work for n <= 50
#    print('omega',omega,'nprmDiff',normDiff,'n',n)    
    omega = omega / (normDiff * (n/30))
#    print('omega',omega,'nprmDiff',normDiff,'n',n)    
    if (omega > normDiff  / (n/30)):
        omega = normDiff  / (n/30)
        
    if (normDiff > 1e1):
        print('Solution changing to fast. Decreasing omega...')
        omega = 0.9 * omegaOld
       
        # If solution is changeing to quickly, take a step back.
        if rank == 0:
            conditions1 = conditions1old
            u1=u1old
        elif rank == 2:
            conditions2 = conditions2old
            u2=u2old
        elif rank == 1:
            u3=u3old
            conditions3 = conditions3old
            
    elif (normDiff < tolerance):
        print('Solution appears to have converged within tolerance: ',tolerance)
        break
   
# Print solution of all rooms. Only done by one processor to prevent muliple plots
if rank == 2:
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
    final[nx:-1,0] = 40
    final[-1,1:nx] = 15
    final[-1,nx:2*nx-2] = 5
    final[1:nx-1,-1] = 40
    final[0,nx:ny2-1] = 40
    final[0,ny2-1:-1] = 15   
    

    fig=matshow(final)
    show()









