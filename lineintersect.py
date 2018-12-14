# python numpy implementation of the following matlab function:
# https://www.mathworks.com/matlabcentral/fileexchange/37192-intersection-point-of-lines-in-3d-space

import numpy as np

def lineintersects(PA, PB):

	# Find intersection point of lines in 3D space, in the least squares sense.
	# PA :          Nx3-matrix containing starting point of N lines
	# PB :          Nx3-matrix containing end point of N 
	# Returns:
	#  P_Intersect : Best intersection point of the N lines, in least squares sense.
	#  distances   : Distances from intersection point to the input lines
  # Original Implementation in MATLAB by Anders Eikenes, 2012
  # Implementation in Python by Aaron Renfroe 2018

	# element wise subtraction
	vectors = PB - PA # N lines as vectors
	ni = vectors / ((vectors**2).sum(axis=1, keepdims=True))**0.5

	#Transposing xyz
	nx = ni.T[0]
	ny = ni.T[1]
	nz = ni.T[2]

	SXX = (nx**2 - 1).sum()
	SYY = (ny**2 - 1).sum()
	SZZ = (nz**2 - 1).sum()
	SXY = (nx*ny).sum()
	SXZ = (nx*nz).sum()
	SYZ = (ny*nz).sum()

	S = np.array([[SXX, SXY, SXZ],[SXY, SYY, SYZ], [SXZ, SYZ ,SZZ]])

	CX = (PA.T[0]*(nx**2-1) + PA.T[1]*(nx*ny) + PA.T[2]*(nx*nz)).sum()
	CY = (PA.T[0]*(nx*ny) + PA.T[1]*(ny**2-1) + PA.T[2]*(ny*nz)).sum()
	CZ = (PA.T[0]*(nx*nz) + PA.T[1]*(ny*nz) + PA.T[2]*(nz**2-1)).sum()

	C = np.array([CX, CY, CZ]).T

	p_intersect = np.linalg.lstsq(S, C)[0] # Freakin Wizardry, Yup

	N = PA.shape[0]
	distances = []
	for i in range(0, N):
		ui = ((p_intersect - PA[i]) * vectors[i] / ((vectors[i]*vectors[i]).T).sum()).sum()
		d = np.linalg.norm(p_intersect-PA[i]-ui*vectors[i])
		distances.append(d)

	return p_intersect, distances


PA = [[400, 100, 100], [-100, 100, 100]]
PB = [[-100, -100, -100], [100, -100, -100]]

# converting from lists to numpy arrays for easier matrix math
PA = np.array(PA)
PB = np.array(PB)

intersects = lineintersects(PA, PB)

print("\n\nIntersect: {}\nDistance to each line: {}\n\n".format(*intersects))
