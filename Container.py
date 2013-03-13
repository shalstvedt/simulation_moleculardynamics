from __future__ import division
import numpy as np

class Container:
    
	def __init__(self, Lx, Ly, x=np.array([]), y=np.array([]), vx=np.array([]), vy=np.array([]), ax=np.array([]), ay=np.array([]), r=np.array([]), m=np.array([])):
		self.Lx = Lx            # length in x-direction of box
		self.Ly = Ly            # length in y-direction of box
		#self.Lz = Lz            # length in z-direction of box
		self.x = x              # needs to be in numpy array form
		self.y = y              # needs to be in numpy array form
		#self.z = z             # likewise
		self.vx = vx            # initial velocity in x-direction
		self.vy = vy            # initial velocity in y-direction
		self.ax = ax
		self.ay = ay
		self.r = r
		self.m = m              # mass of particle
		self.potential = []
		self.kinetic = []
	
	def add_particle(self, x, y, vx, vy, ax, ay, r, m):
		self.x = np.hstack((self.x, x))
		self.y = np.hstack((self.y, y))
		self.vx = np.hstack((self.vx, vx))
		self.vy = np.hstack((self.vy, vy))
		self.ax = np.hstack((self.ax, ax))
		self.ay = np.hstack((self.ay, ay))
		self.r = np.hstack((self.r, r))
		self.m = np.hstack((self.m, m))
		
	def remove_particle(self, index):
		self.x = np.append(self.x[:index], self.x[index+1:])
		self.y = np.append(self.y[:index], self.y[index+1:])
		self.vx = np.append(self.vx[:index], self.vx[index+1:])
		self.vy = np.append(self.vy[:index], self.vy[index+1:])
		self.ax = np.append(self.ax[:index], self.ax[index+1:])
		self.ay = np.append(self.ay[:index], self.ay[index+1:])
		self.r = np.append(self.r[:index], self.r[index+1:])
		self.m = np.append(self.m[:index], self.m[index+1:])
	
	def distance_matrix(self):   # Find position differences
		x_matrix = np.tile(self.x, [self.x.size, 1] ) # where x.size is the number of particles...
		dist_x = x_matrix - np.transpose(x_matrix)
		y_matrix = np.tile(self.y, [self.y.size, 1] ) # where y.size is the number of particles...
		dist_y = y_matrix - np.transpose(y_matrix)
		#z_matrix = np.tile(self.z, [self.z.size, 1] ) # where z.size is the number of particles...
		#dist_z = z_matrix - np.transpose(z_matrix)
		#print dist_x
		#print dist_y
		dist_x[dist_x >  self.Lx/2] -= self.Lx
		dist_x[dist_x < -self.Lx/2] += self.Lx
		dist_y[dist_y >  self.Ly/2] -= self.Ly
		dist_y[dist_y < -self.Ly/2] += self.Ly
		#dist_z[dist_z >  self.Lz/2] -= self.Lz
		#dist_z[dist_z < -self.Lz/2] += self.Lz
		return dist_x, dist_y, # dist_z
        

if __name__ == "__main__":
    width = 6
    height = 6
    x = np.array([1,2,3,4,5], dtype = "float64")
    y = np.array([5,4,3,2,1], dtype = "float64")
    vx = np.array([0,0,0,0,0], dtype = "float64")
    vy = np.array([0,0,0,0,0], dtype = "float64")
    m = np.array([1,1,1,1,1], dtype = "float64")
    c = Container(width, height, x, y, vx, vy, m)
    print c.distance_matrix()