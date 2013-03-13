from __future__ import division
from Container import *
from LennardJonesForce import *
import pylab
# Generic Integrator class
class Integrator:
    pass
 
# Keeps track of the states of a container system over time, and interfaces between the container class and the integration function
class ContainerIntegrator:
	def __init__(self, container, force, t0):
		self.container = container
		self.force = force
		self.time = t0
		self.history = [[container.x, container.y, container.vx, container.vy, container.ax, container.ay]]	
	# runs a number of steps of integration
	def run_until(self, terminal):
		times = np.linspace(self.time, terminal, terminal/self.dt)
		for i in times:
			self.time = i
			#self.time = self.t0 + i*self.dt
			self.verlet_integrator()
	def verlet_integrator(self, dt):
		self.time += dt
		new_x = self.container.x + self.container.vx*dt + 0.5*self.container.ax*dt**2
		new_y = self.container.y + self.container.vy*dt + 0.5*self.container.ay*dt**2
		self.container.x = new_x
		self.container.y = new_y
		new_ax, new_ay = self.force(self.container, self.time) # not sure whether time should be +dt
		new_vx = self.container.vx + (self.container.ax + new_ax)*0.5*dt
		new_vy = self.container.vy + (self.container.ay + new_ay)*0.5*dt
		

		self.container.vx = new_vx
		self.container.vy = new_vy
		self.container.ax = new_ax
		self.container.ay = new_ay
		#self.container.x[self.container.x>self.container.Lx/2] -= self.container.Lx
		self.container.x[self.container.x>self.container.Lx/2] = self.container.x[self.container.x>self.container.Lx/2] % -self.container.Lx
		#self.container.y[self.container.y>self.container.Ly/2] -= self.container.Ly
		self.container.y[self.container.y>self.container.Ly/2] = self.container.y[self.container.y>self.container.Ly/2] % -self.container.Ly
		#self.container.x[self.container.x<-self.container.Lx/2] += self.container.Lx
		self.container.x[self.container.x<-self.container.Lx/2] = self.container.x[self.container.x<-self.container.Lx/2] % self.container.Lx
		#self.container.y[self.container.y<-self.container.Ly/2] += self.container.Ly
		self.container.y[self.container.y<-self.container.Ly/2] = self.container.y[self.container.y<-self.container.Ly/2] % self.container.Ly
		self.history.append([new_x, new_y, new_vx, new_vy, new_ax, new_ay])
 
class FileReader:
	def __init__(self, wd="."):
		if wd==".":
			return
		else:
			pass
			# TODO set working dir to wd
 
	def read_array(self, name):
		return np.genfromtxt(open(os.cwd() + name, "U"), delimiter = " ")
	def read_container(self, name):
		data = self.read_array(name)
		return Container(data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5])
	
def LinearContainerInitializer(Lx, Ly):
	c = Container(Lx, Ly)
	ys = np.linspace(-Ly/2.+1.5, c.Ly/2-1.5, c.Ly/3)
	print "ys: ", ys
	for i in ys:
		if i == 5.0:
			c.add_particle(0, i, 1., 0, 0, 0, 1.5, 2.0)
		else:
			c.add_particle(0, i, 1., 0, 0, 0, 1.5, 2.0)
	return c
	
def TwoParticleInitializer(Lx, Ly):
	c = Container(Lx, Ly)
	dist = c.Lx / 5.
	vel = dist / 5.
	vel = 0.2
	c.add_particle(-dist,0.,vel,0,0.,0.,1.,1.)
	c.add_particle(dist,0.,-vel,0,0.,0.,1.,1.)
	return c
	
def OneParticleInitializer(Lx, Ly):
	c = Container(Lx, Ly)
	dist = c.Lx / 5.
	vel = dist / 5.
	c.add_particle(0, dist, 0.3, 0.5, 0, 0,1.,1.)
	return c
	
def FourParticleInitializer(Lx, Ly):
	c = Container(Lx, Ly)
	dist = c.Lx / 5.
	vel = dist / 5.
	c.add_particle(-dist, 0., vel, 0.1, 0, 0, 1., 1.)
	c.add_particle(dist, 0., -vel, 0.12, 0, 0, 1., 1.)
	c.add_particle(0., dist, 0., -vel, 0, 0, 1., 1.)
	c.add_particle(0., -dist, 0., vel, 0, 0, 1., 1.)
	return c
	
def EightParticleInitializer(Lx, Ly):
	c = Container(Lx, Ly)
	dist = c.Lx / 5.
	vel = dist / 5.
	c.add_particle(-dist,0.,vel,0.,0.,0.,1.,1.)
	c.add_particle(dist,0.,-vel,0.,0.,0.,1.,1.)
	c.add_particle(0.,dist,0,-vel,0.,0.,1.,1.)
	c.add_particle(0.,-dist,0.,vel,0.,0.,1.,1.)
 
	c.add_particle(dist/np.sqrt(2),dist/np.sqrt(2),-vel/np.sqrt(2),-vel/np.sqrt(2),0.,0.,1.,1.)
	c.add_particle(-dist/np.sqrt(2),dist/np.sqrt(2),vel/np.sqrt(2),-vel/np.sqrt(2),0.,0.,1.,1.)
	c.add_particle(-dist/np.sqrt(2),-dist/np.sqrt(2),vel/np.sqrt(2),vel/np.sqrt(2),0.,0.,1.,1.)
	c.add_particle(dist/np.sqrt(2),-dist/np.sqrt(2),-vel/np.sqrt(2),vel/np.sqrt(2),0.,0.,1.,1.)
	
	return c
	
def SquareLatticeInitializer(Lx, Ly):
	c = Container(Lx, Ly)
	
	N = 8             # Particles per row
	d = 2.**(1/6.)    # Particle diameter
	x = np.linspace(-c.Lx/2+d/2.,c.Lx/2-d/2,N)
	y = np.linspace(-c.Lx/2+d/2.,c.Lx/2-d/2,N)
	for i in range(x.size):
		for j in range(y.size):
			c.add_particle(x[i],y[j],0,0,0,0,1,1)
			
	return c
	
def GammaInitializer(Lx, Ly):
	Lx = 10.
	Ly = 10.
	c = Container(Lx, Ly)
	gamma = 1e-6
#	gamma = 0
	for i in range(11):
		if i ==5:
			c.add_particle(c.Lx / 2., (i-.5) * c.Ly / 11., 1.-gamma,gamma,0,0,1.,1.)
		else:
			c.add_particle(c.Lx / 2., (i-.5) * c.Ly / 11.,1.,0.,0,0,1.,1.)
	return c
    
def TriangularLatticeInitializer(Lx, Ly):
    Ly = np.sqrt(3) / 2. * Lx  # Set this based on Lx
    c = Container(Lx, Ly)
    N = 8                       # particles per row
    d = 2.**(1/6.)              # diameter
    x =  np.linspace(-c.Lx/2 + 3.*d/4.,c.Lx/2. - 1.*d/4., N) # Unstaggered
    xs = np.linspace(-c.Lx/2 + d/4.   ,c.Lx/2. - 3.*d/4., N) # Staggered
    y =  np.linspace(-c.Ly/2 + d/2.,c.Ly/2  - d/2, N)
 
    for i in range(N):
        for j in range(N):
            if np.mod(i,2)==0:
                c.add_particle(x[j],y[i],0,0,0,0,1,1)
            else:
                c.add_particle(xs[j],y[i],0,0,0,0,1,1)    
    return c
				
def draw_circle( xy, radius, color="lightsteelblue", facecolor="red", alpha=1.0, ax=None ):
	e = pylab.Circle( xy=xy, radius=radius )
	if ax is None:
		ax = pylab.gca()  # ax = subplot( 1,1,1 )
	ax.add_artist(e)
	e.set_clip_box(ax.bbox)
	e.set_edgecolor( color )
	e.set_linewidth(3)
	e.set_facecolor( facecolor )  # "none" not None
	e.set_alpha( alpha )
	
def plot_integrate(container, integrator):    
	pylab.figure(1, figsize=(10,10))
	pylab.clf()
	pylab.ion()
	pylab.xlim((-container.Lx/2,container.Lx/2))
	pylab.ylim((-container.Ly/2,container.Ly/2))
	pylab.grid()
	pylab.title('Molecular Dynamics')
	pylab.xlabel('Horizontal Distance')
	pylab.ylabel('Vertical Distance')
	ax = pylab.gca()
	pylab.show()
	
	colors = ["red", "blue", "green", "yellow", "orange", "pink", "maroon", "cyan"]
	
	colors = colors*128

	count = 0

	for i in range((int(Tend/dt))):
		integrator.verlet_integrator(dt) # Step forward in time
		if np.mod(count,pack_interval) == 0:
			# c is the container object.
			if count*dt > 4. and np.min([c.Lx,c.Ly]) > 8.*2.**(1/6.):
				# Squeeze the box!
				c.Lx *= squeeze_factor
				c.Ly *= squeeze_factor
				c.x  *=  squeeze_factor
				c.y  *=  squeeze_factor

		print "SQUEEZE", c.Lx, c.Ly
		count+=1
		if np.mod(count,draw_interval) == 0:
			ax.clear()
			for j in range(c.x.size):
				fcol = colors[j]
				r = np.sqrt((container.x - container.x[j])**2 + (container.y - container.y[j])**2)
				r[r==0] = 2.
				r = np.min(r)
				#print "|i|[j]r: |", i, "|[", j, ']', r
				if (r <= 2**(1/6.)):
					fcol = "black"
				draw_circle((container.x[j],container.y[j]),radius = .5*2**(1/6.),ax=ax,facecolor=fcol);
			pylab.title('Molecular Dynamics - Time Forward (' + str(i*dt) + ')')
			pylab.xlabel('Horizontal Distance')
			pylab.ylabel('Vertical Distance')
			pylab.draw()
            			
	# for i in range((int(Tend/dt))-1, -1, -1):
	# 	integrator.verlet_integrator(-dt) # Step forward in time
	# 	count+=1
	# 	if np.mod(count,draw_interval) == 0:
	# 		ax.clear()
	# 		for j in range(c.x.size):
	# 			fcol = colors[j]
	# 			r = np.sqrt((container.x - container.x[j])**2 + (container.y - container.y[j])**2)
	# 			r[r==0] = 2.
	# 			r = np.min(r)
	# 			print "|i|[j]r: |", i, "|[", j, ']', r
	# 			if (r <= 2**(1/6.)):
	# 				fcol = "black"
	# 			draw_circle((container.x[j],container.y[j]),radius = .5*2**(1/6.),ax=ax,facecolor=fcol);
	# 		pylab.title('Molecular Dynamics - Time Backward (' + str(i*dt) + ')')
	# 		pylab.xlabel('Horizontal Distance')
	# 		pylab.ylabel('Vertical Distance')
	# 		pylab.draw()
	   
	pylab.ioff()
	pylab.figure(2,figsize=(20,5))
	pylab.plot(c.potential)
	pylab.title("Potential Energy")
	pylab.show()
	pylab.figure(3,figsize=(20,5))
	pylab.plot(c.kinetic)
	pylab.title("Kinetic Energy")
	pylab.show()
	
### BEGIN MAIN ###
squeeze_factor = 0.98
draw_interval = 5
pack_interval = 20
Tend = 15
t0 = 0.
dt = 0.01
Lx = 15.
Ly = 15.


#fr = FileReader()
#c = fr.read_container("test_file")
f = LennardJonesForce(1.0, 1.0) #sigma, epsilon

c = SquareLatticeInitializer(Lx, Ly)
#c = LinearContainerInitializer(Lx, Ly)
#c = TwoParticleInitializer(Lx, Ly) 
#c = FourParticleInitializer(Lx, Ly)
#c = EightParticleInitializer(Lx, Ly)
#c = GammaInitializer(Lx, Ly) 
#c = TriangularLatticeInitializer(Lx, Ly)
  
integrator = ContainerIntegrator(c, f, t0)

plot_integrate(c, integrator)
 
#integrator.run_until(50.0)

print c.x, c.y

#xplot1 = []
#yplot1 = []
#xplot2 = []
#yplot2 = []
#for i in integrator.history:
#	xplot1.append(i[0][0])
#	yplot1.append(i[0][1])
#	xplot2.append(i[1][0])
#	yplot2.append(i[1][1])
	
	
#pylab.plot(xplot1, yplot1)
#pylab.plot(xplot2, yplot2)
#pylab.figure(1)
#draw_circle((xplot1[-1],yplot1[-1]), f.sigma*2**(1/6))
#draw_circle((xplot2[-1],yplot2[-1]), f.sigma*2**(1/6))
#pylab.xlim((-c.Lx/2,c.Lx/2))
#pylab.ylim((-c.Ly/2,c.Ly/2))
#pylab.show()

#pylab.figure(2)
#draw_circle((xplot1[0],yplot1[0]), f.sigma*2**(1/6))
#draw_circle((xplot2[0],yplot2[0]), f.sigma*2**(1/6))
#pylab.xlim((-c.Lx/2,c.Lx/2))
#pylab.ylim((-c.Ly/2,c.Ly/2))
#pylab.show()
