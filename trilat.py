# source three_tri.m https://www.ece.lsu.edu/scalzo/FDR_Final.pdf
# implemented in python by Aaron Renfroe

import numpy as np
from numpy.linalg import det, lstsq

def r_divide(a,b):
	return lstsq(b.T, a.T)[0].T 

def three_tri(node1, node2, node3, node4):
	[x1, y1, z1, d1] = node1
	[x2, y2, z2, d2] = node2
	[x3, y3, z3, d3] = node3
	[x4, y4, z4, d4] = node4

	#3-D trilateration. 4 anchor nodes
	# x1 = x coordinate of anchor 1
	# y1 = y coordinate of anchor 1
	# z1 = z coordinate of anchor 1
	# d1 = distance from point to anchor 1
	d11 = 2*(x2-x1)
	d21 = 2*(x3-x1)
	d31 = 2*(x4-x1)
	d12 = 2*(y2-y1)
	d22 = 2*(y3-y1)
	d32 = 2*(y4-y1)
	d13 = 2*(z2-z1)
	d23 = 2*(z3-z1)
	d33 = 2*(z4-z1)
	#bringing M together into [3, 3] matrix
	d = np.array([[d11, d12, d13], [d21, d22, d23], [d31, d32, d33]])

	#x_numerator elements
	x_n11 = (d1**2-d2**2) - (x1**2-x2**2) - (y1**2-y2**2) - (z1**2-z2**2) #sigma
	x_n21 = (d1**2-d3**2) - (x1**2-x3**2) - (y1**2-y3**2) - (z1**2-z3**2) #beta
	x_n31 = (d1**2-d4**2) - (x1**2-x4**2) - (y1**2-y4**2) - (z1**2-z4**2) #phi
	x_n12 = 2*(y2-y1)
	x_n22 = 2*(y3-y1)
	x_n32 = 2*(y4-y1)
	x_n13 = 2*(z2-z1)
	x_n23 = 2*(z3-z1)
	x_n33 = 2*(z4-z1)

	#bringing numerator together for x
	x_n = np.array([[x_n11, x_n12, x_n13], [x_n21, x_n22, x_n23], [x_n31, x_n32, x_n33]])
	#finding x by dividing matrix operation and then determinant
	x = r_divide(x_n, d)
	x = det(x)
	#individual y elements
	y_n11 = 2*(x2-x1)
	y_n21 = 2*(x3-x1)
	y_n31 = 2*(x4-x1)
	y_n12 = x_n11 #sigma
	y_n22 = x_n21 #beta
	y_n32 = x_n31 #phi
	y_n13 = 2*(z2-z1)
	y_n23 = 2*(z3-z1)
	y_n33 = 2*(z4-z1)
	#bringing numerator together for y
	y_n = np.array([[y_n11, y_n12, y_n13], [y_n21, y_n22, y_n23], [y_n31, y_n32, y_n33]])
	#finding y by dividing matrix operation and then determinant
	y = r_divide(y_n, d)
	y = det(y)
	#individual z elements
	z_n11 = 2*(x2-x1)
	z_n21 = 2*(x3-x1)
	z_n31 = 2*(x4-x1)
	z_n12 = 2*(y2-y1)
	z_n22 = 2*(y3-y1)
	z_n32 = 2*(y4-y1)
	z_n13 = x_n11 #sigma
	z_n23 = x_n21 #beta
	z_n33 = x_n31 #phi
	#bringing z numerator together
	z_n = np.array([[z_n11, z_n12, z_n13], [z_n21, z_n22, z_n23], [z_n31, z_n32, z_n33]])
	#finding z by dividing matrix operation and then determinant
	z = r_divide(z_n, d)
	z = det(z)

	return x,y,z
