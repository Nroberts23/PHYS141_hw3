"""
Name: Nathan Roberts
PID: A14384608
"""
#imports
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from celluloid import Camera
from mpl_toolkits.mplot3d import Axes3D

#function for returning samples within a given distribution.  Here, I use it to select different radii to place stars at
def gen_samps(x0, x1, func, nsamps):
    smp = []
    while (len(smp) < nsamps):
        x = np.random.uniform(low=x0,high=x1)
        p = func(x)
        assert p>=0 and p<=1
        
        if np.random.uniform(low=0,high=1) <= p:
            smp += [x]
    return smp

def radius_to_stars(radii):
	stars = []
	for r in radii:
		X2 = np.random.uniform(0,1)
		X3 = np.random.uniform(0,1)

		z = (1 - 2 * X2) * r
		x = ((r ** 2 - z ** 2) ** 0.5) * np.cos(2 * np.pi * X3)
		y = ((r ** 2 - z ** 2) ** 0.5) * np.sin(2 * np.pi * X3)

		stars.append([x,y,z])

	return np.array(stars)

def radius_to_vel(radii):
	veloc = []
	for r in radii:
		V_e = np.sqrt(2) * (1 + r ** 2) ** (-0.25)

		g = lambda q: (q ** 2) * (1 - q ** 2) ** (7/2)

		X4 = 0
		X5 = 0
		#keep finding new values until this equality is satisfied
		while (0.1 * X5 >= g(X4)):
			X4 = np.random.uniform(0,1)
			X5 = np.random.uniform(0,1)

		V = X4 * V_e # V / Ve = q = X4

		#convert V to velocity components with random vector
		X6 = np.random.uniform(0,1)
		X7 = np.random.uniform(0,1)
		w = (1 - 2 * X6) * V
		u = (V ** 2 - w ** 2) ** 0.5 * np.cos(2 * np.pi * X7)
		v = (V ** 2 - w ** 2) ** 0.5 * np.sin(2 * np.pi * X7)

		veloc.append([w, u, v])

	return np.array(veloc)

def plot_3d_dist(star_x):
	fig = plt.figure(figsize=(13,13))
	ax = fig.add_subplot(111, projection='3d')

	for star in star_x:
		x = star[0]
		y = star[1]
		z = star[2]

		ax.scatter3D(xs=[x], ys=[y], zs=[z], color='white', s=1)
		#ax.plot(xs=x_history[i][-1], ys=v_history[i][-1], zs=t_history[-1])
	ax.set_facecolor('black')
	ax.set_title('Plummer Sphere')

	plt.savefig('plummer_pt1.png')

def sphere_density_plot(rs, R):
	#given the radii of all of our stars, plot the density of our model
	fig = plt.figure(figsize=(6,4))
	ax = fig.add_subplot(111)
	ax.hist(rs, bins=200)
	ax.set_xlabel('Radius (pc)')
	ax.set_ylabel('Mass (2x10^7 M_solar)')
	ax.set_xlim(0, max(rs))
	ax.vlines(R, ymin=0, ymax=200)
	ax.set_ylim(0,170)
	ax.set_title('Mass Distribution of Generated Sphere')
	plt.show()


def main(print_imgs=False):
	#initial conditions

	R = 1400 #pc
	M = 2*10**11 #solar masses

	x0 = 0
	x1 = 2 * R
	dens_func = lambda x: (3/(4* np.pi) * M * R ** (-3) * (1 + (x / R) ** 2)** (-5/2)) / (3/(4* np.pi) * M * R ** (-3)) #Density function normalized so the max = 1
	n_stars = 10000
	m = M / n_stars # the mass of a single point in solar masses

	star_rs = gen_samps(x0, x1, dens_func, n_stars)
	star_x = radius_to_stars(star_rs)
	star_v = radius_to_vel(star_rs)

	x_history = [[i] for i in star_x]
	v_history = [[i] for i in star_v]
	t_history = [0]

	plot_3d_dist(star_x)

	theor_dens_func = lambda x: (3/(4* np.pi) * M * R ** (-3) * (1 + (x / R) ** 2)** (-5/2))	
	X = np.linspace(0,2*R,30000)
	fig = plt.figure(figsize=(6,4))
	ax = fig.add_subplot(111)
	ax.plot(X,theor_dens_func(X))
	ax.vlines(R, ymin=0, ymax=18)
	ax.set_xlim(0,max(star_rs))
	ax.set_ylim(0,18)
	ax.set_xlabel('Radius (pc)')
	ax.set_ylabel('Density (Msolar / pc^3)')
	ax.set_title('Theoretical Plummer Sphere Density')
	plt.show()

	sphere_density_plot(star_rs, R)

