# import python modules
import sys
import bempp.api
import numpy as np
from cfdbem.file_interfaces import FileReader

# input data
c = int(sys.argv[1])
freq = int(sys.argv[2])
dataName = sys.argv[3]

epsilon = 1E-5
k = freq * 2 * np.pi / c;
print ("k = ", k)

# numerical constant for combined formulation
muD = 1.0/k


# mesh
reader = FileReader(file_name = dataName)
#reader = FileReader(file_name = "all.msh")

grid = reader.grid

import time
time1 = time.time()

# spaces
piecewise_lin_space = bempp.api.function_space(grid, "P", 1)
piecewise_const_space = bempp.api.function_space(grid, "DP", 0)

domain_space = piecewise_lin_space 
range_space = piecewise_lin_space 
dual_to_range_space = piecewise_lin_space 


# operators
identity = bempp.api.operators.boundary.sparse.identity(
    domain_space, range_space, dual_to_range_space)
dlp = bempp.api.operators.boundary.helmholtz.double_layer(
    domain_space, range_space, dual_to_range_space,k)
slp = bempp.api.operators.boundary.helmholtz.single_layer(
    domain_space, range_space, dual_to_range_space,k)

#print ("ok")

#blp = bempp.api.operators.boundary.helmholtz.single_layer(
#    piecewise_lin_space, piecewise_lin_space, piecewise_lin_space,0)
        
#print ("ok")

# dirichlet BC
dirichlet_fun = bempp.api.GridFunction(piecewise_lin_space, coefficients = reader.node_data)

# equation in combined formulation!!!!!
lhs = (.5*identity + dlp) - 1j * muD * slp

time2 = time.time()

print ("prepared to ", time2 - time1, " s")


# solve
from bempp.api.linalg import gmres

it_count = 0
def iteration_counter(x):
    global it_count
    it_count += 1
    
w_fun,info,resid, it_count = gmres(lhs, dirichlet_fun, tol=epsilon, return_residuals = True, return_iteration_count = True)
print("The linear system was solved in {0} iterations".format(it_count))

time3 = time.time()

print ("solved to ", time3 - time2, " s")



# compute and output result


def result (points):
    from bempp.api.operators.potential import helmholtz as helmholtz_potential

    slp_pot=helmholtz_potential.single_layer(piecewise_lin_space, points, k)
    dlp_pot=helmholtz_potential.double_layer(piecewise_lin_space, points, k)

    res =  1j * muD * slp_pot.evaluate(w_fun) - dlp_pot.evaluate(w_fun)
    
    return res

# plot along X-axis in given time point
def plotPerRadius (time, R1, R2, h):
    
    coord = np.arange(R1,R2+h,h,dtype=np.float)
    points = np.vstack((coord.ravel(), 0*np.ones(coord.size), 0*np.ones(coord.size)))
    resR= result(points)

    fileName = "pressureR.dat"
    outputFile = open(fileName, 'w')

    i = 0
    for num in resR[0]:
	cur = abs(num)
	outputFile.write('%s' % coord[i] + ' ' + '%s' % cur)
	outputFile.write('\n')
	i = i + 1

    outputFile.close()

# sound pressure per time [t0, T + t0] in given point
def plotPerTime (point, t0, T, tau, outFile):
    
    p = np.vstack(point)
    resTime = result(p)

    fOutput = open(outFile,"w")

    for t in np.arange(t0, T + t0 + tau, tau):
	
	cur1 = np.real(resTime[0][0] * np.exp (- 1j * 2 * np.pi * freq * (t)))
	fOutput.write('%s' % t + ' ' + '%s' % cur1)
	fOutput.write('\n')
    
    fOutput.close()

# sound pressure: polar plot in Oxz
def polarplot(R):
    theta = np.arange(0,2*np.pi,0.01,dtype = np.float)

    x1 = R*np.cos(theta)
    z1 = R*np.sin(theta)
    polarPoints = np.vstack(( x1.ravel(), 0*np.ones(theta.size), z1.ravel() ))

    polarRes = result(polarPoints)

    fileName = "polarR" + '%s'%R + ".dat"
    polarFile = open(fileName,"w")

    i = 0
    for num in polarRes[0]:
	cur = abs(num)
	polarFile.write('%s' % theta[i] + ' ' + '%s' % cur)
	polarFile.write('\n')
	i = i + 1

    polarFile.close()
    print ("polar OK")
    
    
plotPerTime((0,0,1),0,0.1,0.0001,"pressureTime1.dat")
plotPerTime((0,0,2),0,0.1,0.0001,"pressureTime2.dat")
plotPerTime((0,0,3),0,0.1,0.0001,"pressureTime3.dat")
plotPerTime((0,0,4),0,0.1,0.0001,"pressureTime4.dat")
plotPerTime((0,0,5),0,0.1,0.0001,"pressureTime5.dat")




