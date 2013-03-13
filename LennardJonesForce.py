from __future__ import division
import numpy as np
from Container import *

class LennardJonesForce:
	def __init__(self, sigma, epsilon):
		self.sigma = sigma
		self.epsilon = epsilon
	def __call__(self,c,t):
		dx, dy = c.distance_matrix() # dx & dy have 0's in the diagonal
		dr = np.sqrt(dx**2 + dy**2)
		dr += np.identity((dr.size)**0.5) # adds identity matrix, removing 0's on the diagonal for division
		#print dr
		sigma, epsilon = self.sigma, self.epsilon

		m1 = (sigma/dr)**12 - (sigma/dr)**6
		
		m2 = 2*(sigma/dr)**12 - (sigma/dr)**6
		
		# Force magnitude:
		#magnitude = self.G * c.m / dr ** 2
		magnitude = 24*epsilon/dr * m2
		#magnitude -= np.identity((magnitude.size)**0.5)*magnitude # removes nonzero diagonals so that a particle does not interact with itself
		#magnitude -= n

		pe = 4 * epsilon * np.sum(np.triu(m1,k=1))
		c.potential.append(pe)
		
		ke = sum(.5 * (c.vx**2 + c.vy**2))
		c.kinetic.append(ke)

		# Project onto components, sum all forces on each particle 
		ax =  np.sum(-magnitude * dx/dr, axis=1) # over dr cancels
		ay =  np.sum(-magnitude * dy/dr, axis=1)
		return ax, ay
		
if __name__ == "__main__":
	width = 6
	height = 6
	x = np.array([1,2,3,4,5], dtype = "float64")
	y = np.array([5,4,3,2,1], dtype = "float64")
	vx = np.array([0,0,0,0,0], dtype = "float64")
	vy = np.array([0,0,0,0,0], dtype = "float64")
	m = np.array([1,1,1,1,1], dtype = "float64")
	c = Container(width, height, x, y, vx, vy, m)
	#print c.distance_matrix()
	f = LennardJonesForce(1.0, 1.0)
	fj = JesseLennardJones(1.0, 1.0)
	#print f(c,0)
	print fj(c)
