import numpy as np
import math, time
import random
import matplotlib.pyplot as plt
from numpy import linalg as LA
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib.animation import FuncAnimation
import matplotlib.patches as mpatches
from scipy import stats

import os
import sqlite3

k_B = 1.380648e-23 #Boltzmann constant (J/K)

current_milli_time = lambda: int(round(time.time() * 1000))

dir_principal = os.getcwd()
carpeta_data = dir_principal + '\Data'

###########################################################
#														  #
#						FUNCTIONS						  #
#													  	  #
###########################################################

def mod(v):
    """
        Computes the squared sum over the last axis of the numpy.ndarray v
    """
    return np.sum(v * v, axis=-1)

def calculate_T(v, m, N):
	suma = np.sum(v * v)

	temperature = suma * m / (3*N * k_B)

	return temperature

def calculate_P(N, x, T, dt, p):
	average_p = np.sum(p) / len(p)
	Pressure = average_p / (dt * 6 * x**2)

	theory = N * k_B * T / (x**3)

	return (Pressure, theory)

def vel_hist(x, v, v_max, dv, m, T):

	def fmb(m, T, v):
		return (m/(2.*np.pi*k_B*T))**1.5 * 4.*np.pi*v*v * np.exp((-m*v*v)/(2.*k_B*T))  #Maxwell-Boltzman distribution
 
	v = np.arange(0.0, v_max, dv)

	plt.plot(v, fmb(m, T, v), 'r-')
	n, bins, patches = plt.hist(x, 40, normed=1, facecolor='g', alpha=0.5)
	plt.xlabel('Speed -- v (m/s)')
	plt.ylabel('Probability -- f(v)')
	plt.title('Histogram of Maxwell-Boltzmann speed distribution')
	plt.grid(True)
	
	os.chdir(carpeta_data)
	plt.savefig('Velocity_hist.png')
	plt.gcf().clear()

def create_database():
	database = sqlite3.connect("gas_database.db")
	cursor = database.cursor()

	tablas = [
			"""
				CREATE TABLE IF NOT EXISTS ideal_gas(
					temperature REAL NOT NULL,
					theo_pressure REAL NOT NULL,
					exp_pressure REAL NOT NULL,
					rel_error REAL NOT NULL
				);
			"""
		]
	for tabla in tablas:
		cursor.execute(tabla);

	tablas = [
			"""
				CREATE TABLE IF NOT EXISTS real_gas(
					temperature REAL NOT NULL,
					theo_pressure REAL NOT NULL,
					exp_pressure REAL NOT NULL,
					rel_error_pressure REAL NOT NULL,
					theo_avg_vel REAL NOT NULL,
					exp_avg_vel REAL NOT NULL,
					rel_error_vel REAL NOT NULL
				);
			"""
		]
	for tabla in tablas:
		cursor.execute(tabla);

	database.commit()

def fill_database(name, table, interrogants, values):
	database = sqlite3.connect('%s.db' % name)
	cursor = database.cursor()

	sentencia = "INSERT INTO %s VALUES %s" % (table, interrogants) 
	cursor.execute(sentencia, values)
	database.commit()
	
def read_database(name, table, name_2, order):
	database = sqlite3.connect('%s.db' % name)
	cursor = database.cursor()
	sentencia = "SELECT * FROM '%s' ORDER BY %s %s;" % (table, name_2, order)
 
	cursor.execute(sentencia)
	
	lines = cursor.fetchall()

	cursor.close()
	
	return lines

###########################################################
#														  #
#						CLASSES						  	  #
#													  	  #
###########################################################

class IdealGas(object):
	def __init__(self, N, L, T, m, R, dt, T_MAX):
		self.N = N #Number of particles
		self.L = L #Length of the box
		self.T = T #Temperature of the system
		self.m = m #Particle mass (Oxygen  2.66e-23) (He 6.64e-24)
		self.R = R #Particle radius
		self.dt = dt #Step size
		self.T_MAX = T_MAX #Maximum simulation time

		self.Nt = int(self.T_MAX / self.dt) #Number of simulation steps/iterations

		self.r = np.zeros((self.N, 3)) #Array of vector positions
		self.v = np.zeros((self.N, 3)) #Array of vector velocities
		self.p = np.zeros(int(self.Nt)) # Array of momentum exchanged with the wall
		self.CM = (0,0,0) #Center of masses

		self.CM_t = []
		
	def init_particles(self):

		self.r = np.random.rand(self.N, 3) * 2 * (self.L/2 - self.R) - (self.L/2 - self.R)

		#vmed = np.sqrt(8 * k_B * T / (np.pi * m))
		vrms = np.sqrt(3*k_B * self.T / self.m)

		v_polar = np.random.random((self.N, 2))

		self.v[:,0] = (np.sin(v_polar[:,0] * np.pi) * np.cos(v_polar[:,1] * 2 * np.pi)) 
		self.v[:,1] = (np.sin(v_polar[:,0] * np.pi) * np.sin(v_polar[:,1] * 2 * np.pi)) 
		self.v[:,2] = np.cos(v_polar[:,0] * np.pi)

		self.v *= vrms

	def check_walls(self, i, x): #Checks for particle-wall collisions and return the exchanged momentum
		pex = 0
		for k in range(3):
			if abs(self.r[i][k]) + self.R > x/2:
				self.v[i][k] *= -1
				if self.r[i][k] < 0:
					self.r[i][k] = -x/2 + self.R
				else :
					self.r[i][k] = x/2 - self.R

				pex += 2 * self.m * abs(self.v[i][k])
		return pex

	def simulation(self):
		self.init_particles()

		for k in range(self.Nt):
			for i in range(self.N):
				self.p[k] += self.check_walls(i, self.L)
				self.r[i] += self.dt * self.v[i]

			self.CM_t.append(np.sum(self.r, axis=0) / self.N)

		return self.p

	def isothermal_expansion(self, pull_period, L0, dL):
		self.init_particles()

		Vs = []
		Vs_discrete = []
		theo_pr_discrete = []
		time = []
		Exp_pr = []
		Theo_pr = []


		Vs_discrete.append(L0**3)	
		theo_pr_discrete.append(self.N * k_B * self.T / (L0**3))

		self.L = L0

		for k in range(self.Nt):
			if k % (pull_period/self.dt) == 0 and k > 0:
				self.L = self.L + dL
				Vs_discrete.append(self.L**3)
				theo_pr_discrete.append(self.N * k_B* self.T/(self.L**3))

			for i in range(self.N):
				self.p[k] += self.check_walls(i, self.L)
				self.r[i] += self.dt * self.v[i]

			if k > 1:

				Exp_P = self.p[k] / (self.dt * 6 * self.L**2)
				Theo_P = self.N * k_B * self.T / (self.L**3)

				Exp_pr.append(Exp_P)
				Theo_pr.append(Theo_P)
				Vs.append(self.L**3)
				time.append(k*self.dt)

		if pull_period > self.dt :		
			a = int(self.T_MAX/pull_period) #Times the pressure will remain constant (as V is doesn't change)

			average_t = [] #average of the pressure in this times

			for h in range(a):
				avrg = 0
				realtime = int(pull_period/self.dt)
				for s in range(h*realtime, (h+1)*realtime-2):
					avrg = avrg + Exp_pr[s]
				average_t.append(avrg/realtime)

		plt.plot(Vs, Theo_pr, ':o', label = 'Theoretical pressure')
		plt.plot(Vs, Exp_pr, 'x', label = 'Experimental pressure')

		if pull_period > self.dt :
			plt.plot(Vs_discrete, average_t, ls = 'none', marker = 'o', label = 'Average experimental pressure')

		plt.title('Isothermal compression (T=%.2f K)' % self.T)
		plt.xlabel('Volume ($m^3$)')
		plt.ylabel('Pressure (Pa)')
		plt.legend()
		
		os.chdir(carpeta_data)
		plt.savefig('Expansion_volume.png')
		plt.gcf().clear()

		plt.plot(time, Theo_pr, '-', label = 'Theoretical pressure')
		plt.plot(time, Exp_pr, 'x', label = 'Experimental pressure')

		if pull_period > self.dt :
			plt.plot(np.linspace(0+ pull_period/2, self.T_MAX - pull_period/2, (a)), average_t, ls = 'none', marker='o', label = 'Average experimental pressure')

		plt.title('Isothermal compression (T=%.2f K)' % self.T)
		plt.xlabel('time ($s$)')
		plt.ylabel('Pressure (Pa)')
		plt.legend()
		
		os.chdir(carpeta_data)
		plt.savefig('Expansion_time.png')
		plt.gcf().clear()			

	def isochoric_heating(self, heat_period, T0, dT):	
		T_discrete = []
		T_discrete.append(T0)
		counter =0

		Exp_pr = []
		Theo_pr = []
		temperatures = []
		time = []

		self.init_particles()

		for k in range(self.Nt):

			if k % (heat_period / self.dt) == 0 and k > 0:
				self.T = self.T + dT
				T_discrete.append(self.T)

				vmed = np.sqrt(8 * k_B * self.T / (np.pi * self.m)) #Re-init particles velocities for each T
				v_polar = np.random.random((self.N, 2))

				self.v[:,0] = (np.sin(v_polar[:,0] * np.pi) * np.cos(v_polar[:,1] * 2 * np.pi)) 
				self.v[:,1] = (np.sin(v_polar[:,0] * np.pi) * np.sin(v_polar[:,1] * 2 * np.pi)) 
				self.v[:,2] = np.cos(v_polar[:,0] * np.pi)

				self.v *= vmed

			for i in range(self.N):
				self.p[k] += self.check_walls(i, self.L)
				self.r[i] += self.dt * self.v[i]

			Exp_P = self.p[k]/(self.dt * 6 * self.L**2)
			Theo_P = self.N * k_B * self.T / (self.L**3)

			Exp_pr.append(Exp_P)
			Theo_pr.append(Theo_P)
			temperatures.append(self.T)
			time.append(k * self.dt)

		if heat_period > self.dt:
			a = int(self.T_MAX / heat_period) #Times the pressure will remain constant (as T doesn't change)

			average_t = [] #average of the pressure in this times

			for h in range(a):
				avrg = 0
				realtime = int(heat_period / self.dt)
				for s in range(h*realtime, (h+1)*realtime-2):
					avrg = avrg + Exp_pr[s]
				average_t.append(avrg/realtime)

		plt.plot(temperatures, Theo_pr, ':.', label = 'Theoretical pressure')
		plt.plot(temperatures, Exp_pr, 'x', label = 'Experimental pressure')

		if heat_period > self.dt:
			plt.plot(T_discrete, average_t, '-o', label = 'Average experimental pressure')

		plt.title('Isochoric process (V=%.2f $m^3$)'%(self.L**3))
		plt.xlabel('Temperature ($K$)')
		plt.ylabel('Pressure (Pa)')
		plt.legend()

		os.chdir(carpeta_data)
		plt.savefig('Heating_temperature.png')
		plt.gcf().clear()

		plt.plot(time, Theo_pr, '-', label = 'Theoretical pressure')
		plt.plot(time, Exp_pr, 'x', label = 'Experimental pressure')

		if heat_period > self.dt:
			plt.plot(np.linspace(0+ heat_period/2, self.T_MAX - heat_period/2, (a)), average_t, ls = 'none', marker='o', label = 'Average experimental pressure')

		plt.title('Isochoric process (V=%.2f $m^3$)'%(self.L**3))
		plt.xlabel('time ($s$)')
		plt.ylabel('Pressure (Pa)')
		plt.legend()

		os.chdir(carpeta_data)
		plt.savefig('Heating_time.png')
		plt.gcf().clear()

	def animation(self):

		self.init_particles()

		def update(t, lines):
			k = int(t / self.dt)

			for i in range(self.N):
				self.p[k] += self.check_walls(i, self.L)
				self.r[i] += self.dt * self.v[i]

			self.CM = np.sum(self.r, axis=0) / self.N
			self.CM_t.append(self.CM)

			lines[0].set_data(self.r[:,0], self.r[:,1])
			lines[0].set_3d_properties(self.r[:,2])

			lines[1].set_data([self.CM[0]], [self.CM[1]])
			lines[1].set_3d_properties([self.CM[2]])

			return lines

		# Attaching 3D axis to the figure
		fig = plt.figure()
		ax = p3.Axes3D(fig)


		# Setting the axes properties
		ax.set_xlim3d([-self.L/2, self.L/2])
		ax.set_xlabel('X')

		ax.set_ylim3d([-self.L/2, self.L/2])
		ax.set_ylabel('Y')

		ax.set_zlim3d([-self.L/2, self.L/2])
		ax.set_zlabel('Z')

		T1 = calculate_T(self.v, self.m, self.N)

		lines = []
		lines.append(ax.plot(self.r[:,0], self.r[:,1], self.r[:,2], ls='None', marker='.', label = 'Particles at $T_0=%.2f\,K$'%T1)[0])
		lines.append(ax.plot([self.CM[0]], [self.CM[1]], [self.CM[2]], marker='o', color='r', label='Center of masses')[0])

		ani = FuncAnimation(fig, update, fargs=(lines,), frames=np.linspace(0, self.T_MAX-self.dt, int(self.T_MAX / self.dt)),
		                    blit=True, interval=1, repeat=False)
		#ani.save('animation.mp4', fps=20, writer="ffmpeg", codec="libx264")
		plt.legend(loc = 'upper left')
		plt.show()

	def experiments(self):

		Exp_P, Theo_P = calculate_P(self.N, self.L, self.T, self.dt, self.p)
		T1 = calculate_T(self.v, self.m, self.N)

		A = 6*self.L**2

		theo_P_2 = self.N * k_B * T1 / (self.L**3)

		rel_error = abs(Theo_P - Exp_P) / Theo_P * 100
		rel_error2 = abs(theo_P_2 - Exp_P) / theo_P_2 * 100

		os.chdir(carpeta_data)
		fill_database('gas_database', 'ideal_gas', '(?,?,?,?)', [T1, Theo_P, Exp_P, rel_error])

		#Pressure plot
		theo_array = []
		exp_array = []

		for i in range(self.Nt):
			theo_array.append(Theo_P)
			exp_array.append(Exp_P)

		plt.plot(np.linspace(0, self.T_MAX - self.dt, self.Nt), self.p / (self.dt * A), ':', label = 'Experimental pressure')
		plt.plot(np.linspace(0, self.T_MAX - self.dt, self.Nt), theo_array, label = 'Theoretical pressure')
		plt.plot(np.linspace(0, self.T_MAX - self.dt, self.Nt), exp_array, label = 'Average experimental pressure')
		plt.title('Pressure of the system')
		plt.xlabel('t(s)')
		plt.ylabel('P(Pa)')
		plt.legend()

		os.chdir(carpeta_data)
		plt.savefig('Pressure.png')
		plt.gcf().clear()

		#CM PLOTS
		transposed = np.transpose(self.CM_t)
		t = np.linspace(0, self.T_MAX - self.dt, len(self.CM_t))

		std_x = np.std(transposed[:][0])
		std_y = np.std(transposed[:][1])
		std_z = np.std(transposed[:][2])

		mean_x = np.mean(transposed[:][0])
		mean_y = np.mean(transposed[:][1])
		mean_z = np.mean(transposed[:][2])

		plt.subplot(2, 2, 1)
		plt.plot(t, transposed[:][0]) #CM_x
		meanx = plt.plot(t, np.linspace(mean_x, mean_x, len(self.CM_t)), color='red', label='Mean') #mean
		sdx = plt.plot(t, np.linspace(std_x, std_x, len(self.CM_t)), ls = ':', color='green', label=r'$\pm$ SD') #+ standard deviation
		plt.plot(t, np.linspace(-std_x, -std_x, len(self.CM_t)), ls = ':', color='green') #- standard deviation
		plt.xlabel('t')
		plt.ylabel('x')

		plt.subplot(2, 2, 2)
		plt.plot(t, transposed[:][1]) #CM_y
		plt.plot(t, np.linspace(mean_y, mean_y, len(self.CM_t)), color='red', label='Mean') #mean
		plt.plot(t, np.linspace(std_y, std_y, len(self.CM_t)), ls = ':', color='green', label=r'$\pm$ SD') #+ standard deviation
		plt.plot(t, np.linspace(-std_y, -std_y, len(self.CM_t)), ls = ':', color='green') #- standard deviation
		plt.xlabel('t')
		plt.ylabel('y')

		plt.subplot(2, 2, 3)
		plt.plot(t, transposed[:][2]) #CM_z
		plt.plot(t, np.linspace(mean_z, mean_z, len(self.CM_t)), color='red', label='Mean') #mean
		plt.plot(t, np.linspace(std_z, std_z, len(self.CM_t)), ls = ':', color='green', label=r'$\pm$ SD') #+ standard deviation
		plt.plot(t, np.linspace(-std_z, -std_z, len(self.CM_t)), ls = ':', color='green') #- standard deviation
		plt.xlabel('t')
		plt.ylabel('z')

		plt.suptitle(r'$\vec{r}_{CM}$ vs t')
		plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
		
		os.chdir(carpeta_data)
		plt.savefig('CM_plots.png')
		plt.gcf().clear()

class RealGas(object):
	def __init__(self, N, L, T, m, R, dt, T_MAX):
		self.N = N #Number of particles
		self.L = L #Length of the box
		self.T = T #Temperature of the system
		self.m = m #Particle mass (Oxygen  2.66e-23) (He 6.64e-24)
		self.R = R #Particle radius
		self.dt = dt #Step size
		self.T_MAX = T_MAX #Maximum simulation time

		self.Nt = int(self.T_MAX / self.dt) #Number of simulation steps/iterations

		self.r = np.zeros((self.N, 3)) #Array of vector positions
		self.v = np.zeros((self.N, 3)) #Array of vector velocities
		self.p = np.zeros(int(self.Nt)) # Array of momentum exchanged with the wall
		self.CM = (0,0,0) #Center of masses

		self.hist = []
		self.CM_t = []

	def init_particles(self):

		self.r = np.random.rand(self.N, 3) * 2 * (self.L/2 - self.R) - (self.L/2 - self.R)

		#vmed = np.sqrt(8 * k_B * T / (np.pi * m))
		vrms = np.sqrt(3*k_B * self.T / self.m)

		v_polar = np.random.random((self.N, 2))

		self.v[:,0] = (np.sin(v_polar[:,0] * np.pi) * np.cos(v_polar[:,1] * 2 * np.pi)) 
		self.v[:,1] = (np.sin(v_polar[:,0] * np.pi) * np.sin(v_polar[:,1] * 2 * np.pi)) 
		self.v[:,2] = np.cos(v_polar[:,0] * np.pi)

		self.v *= vrms

	def check_walls(self, i, x): #Checks for particle-wall collisions and return the exchanged momentum
		pex = 0
		for k in range(3):
			if abs(self.r[i][k]) + self.R > x/2:
				self.v[i][k] *= -1
				if self.r[i][k] < 0:
					self.r[i][k] = -x/2 + self.R
				else :
					self.r[i][k] = x/2 - self.R

				pex += 2 * self.m * abs(self.v[i][k])
		return pex

	def check_collision(self):  #Check for collisions with other particles of te same gas

		dists = np.sqrt(mod(self.r - self.r[:,np.newaxis])) 
		cols2 = (0 < dists) & (dists < 2*self.R)
		idx_i, idx_j = np.nonzero(cols2)
		# ***possibility to simplify this *** #
		for i, j in zip(idx_i, idx_j):
		    if j < i:
		        # skip duplications and same particle
		        continue 

		    rij = self.r[i] - self.r[j]
		    d = mod(rij)
		    vij = self.v[i] - self.v[j]
		    dv = np.dot(vij, rij) * rij / d
		    self.v[i] -= dv
		    self.v[j] += dv
			
		    # update the positions so they are no longer in contact
		    self.r[i] += self.dt * self.v[i]
		    self.r[j] += self.dt * self.v[j]

	def simulation(self):
		self.init_particles()

		for k in range(self.Nt):
			self.check_collision()

			modv = 0
			for i in range(self.N):
				self.p[k] += self.check_walls(i, self.L)
				self.r[i] += self.dt * self.v[i]
				modv = modv + LA.norm(self.v[i])

			self.CM_t.append(np.sum(self.r, axis=0) / self.N)

			average_v = modv / self.N
			self.hist.append(average_v)

		return self.p

	def isothermal_expansion(self, pull_period, L0, dL):
		self.init_particles()

		Vs = []
		Vs_discrete = []
		theo_pr_discrete = []
		time = []
		Exp_pr = []
		Theo_pr = []


		Vs_discrete.append(L0**3)	
		theo_pr_discrete.append(self.N * k_B * self.T / (L0**3))

		self.L = L0

		for k in range(self.Nt):
			self.check_collision()

			if k % (pull_period/self.dt) == 0 and k > 0:
				self.L = self.L + dL
				Vs_discrete.append(self.L**3)
				theo_pr_discrete.append(self.N * k_B* self.T/(self.L**3))

			for i in range(self.N):
				self.p[k] += self.check_walls(i, self.L)
				self.r[i] += self.dt * self.v[i]

			if k > 1:

				Exp_P = self.p[k] / (self.dt * 6 * self.L**2)
				Theo_P = self.N * k_B * self.T / (self.L**3)

				Exp_pr.append(Exp_P)
				Theo_pr.append(Theo_P)
				Vs.append(self.L**3)
				time.append(k*self.dt)

		if pull_period > self.dt :		
			a = int(self.T_MAX/pull_period) #Times the pressure will remain constant (as V is doesn't change)

			average_t = [] #average of the pressure in this times

			for h in range(a):
				avrg = 0
				realtime = int(pull_period/self.dt)
				for s in range(h*realtime, (h+1)*realtime-2):
					avrg = avrg + Exp_pr[s]
				average_t.append(avrg/realtime)

		plt.plot(Vs, Theo_pr, ':o', label = 'Theoretical pressure')
		plt.plot(Vs, Exp_pr, 'x', label = 'Experimental pressure')

		if pull_period > self.dt :
			plt.plot(Vs_discrete, average_t, ls = 'none', marker = 'o', label = 'Average experimental pressure')

		plt.title('Isothermal compression (T=%.2f K)' % self.T)
		plt.xlabel('Volume ($m^3$)')
		plt.ylabel('Pressure (Pa)')
		plt.legend()
		
		os.chdir(carpeta_data)
		plt.savefig('Expansion_volume.png')
		plt.gcf().clear()

		plt.plot(time, Theo_pr, '-', label = 'Theoretical pressure')
		plt.plot(time, Exp_pr, 'x', label = 'Experimental pressure')

		if pull_period > self.dt :
			plt.plot(np.linspace(0+ pull_period/2, self.T_MAX - pull_period/2, (a)), average_t, ls = 'none', marker='o', label = 'Average experimental pressure')

		plt.title('Isothermal compression (T=%.2f K)' % self.T)
		plt.xlabel('time ($s$)')
		plt.ylabel('Pressure (Pa)')
		plt.legend()
		
		os.chdir(carpeta_data)
		plt.savefig('Expansion_time.png')
		plt.gcf().clear()			

	def isochoric_heating(self, heat_period, T0, dT):	
		T_discrete = []
		T_discrete.append(T0)
		counter =0

		Exp_pr = []
		Theo_pr = []
		temperatures = []
		time = []

		self.init_particles()

		for k in range(self.Nt):
			self.check_collision()

			if k % (heat_period / self.dt) == 0 and k > 0:
				self.T = self.T + dT
				T_discrete.append(self.T)

				vmed = np.sqrt(8 * k_B * self.T / (np.pi * self.m)) #Re-init particles velocities for each T
				v_polar = np.random.random((self.N, 2))

				self.v[:,0] = (np.sin(v_polar[:,0] * np.pi) * np.cos(v_polar[:,1] * 2 * np.pi)) 
				self.v[:,1] = (np.sin(v_polar[:,0] * np.pi) * np.sin(v_polar[:,1] * 2 * np.pi)) 
				self.v[:,2] = np.cos(v_polar[:,0] * np.pi)

				self.v *= vmed

			for i in range(self.N):
				self.p[k] += self.check_walls(i, self.L)
				self.r[i] += self.dt * self.v[i]

			Exp_P = self.p[k]/(self.dt * 6 * self.L**2)
			Theo_P = self.N * k_B * self.T / (self.L**3)

			Exp_pr.append(Exp_P)
			Theo_pr.append(Theo_P)
			temperatures.append(self.T)
			time.append(k * self.dt)

		if heat_period > self.dt:
			a = int(self.T_MAX / heat_period) #Times the pressure will remain constant (as T doesn't change)

			average_t = [] #average of the pressure in this times

			for h in range(a):
				avrg = 0
				realtime = int(heat_period / self.dt)
				for s in range(h*realtime, (h+1)*realtime-2):
					avrg = avrg + Exp_pr[s]
				average_t.append(avrg/realtime)

		plt.plot(temperatures, Theo_pr, ':.', label = 'Theoretical pressure')
		plt.plot(temperatures, Exp_pr, 'x', label = 'Experimental pressure')

		if heat_period > self.dt:
			plt.plot(T_discrete, average_t, '-o', label = 'Average experimental pressure')

		plt.title('Isochoric process (V=%.2f $m^3$)'%(self.L**3))
		plt.xlabel('Temperature ($K$)')
		plt.ylabel('Pressure (Pa)')
		plt.legend()

		os.chdir(carpeta_data)
		plt.savefig('Heating_temperature.png')
		plt.gcf().clear()

		plt.plot(time, Theo_pr, '-', label = 'Theoretical pressure')
		plt.plot(time, Exp_pr, 'x', label = 'Experimental pressure')

		if heat_period > self.dt:
			plt.plot(np.linspace(0+ heat_period/2, self.T_MAX - heat_period/2, (a)), average_t, ls = 'none', marker='o', label = 'Average experimental pressure')

		plt.title('Isochoric process (V=%.2f $m^3$)'%(self.L**3))
		plt.xlabel('time ($s$)')
		plt.ylabel('Pressure (Pa)')
		plt.legend()

		os.chdir(carpeta_data)
		plt.savefig('Heating_time.png')
		plt.gcf().clear()


	def animation(self):

		self.init_particles()

		def update(t, lines):
			k = int(t / self.dt)
			self.check_collision(self.r, self.v)

			modv = 0
			for i in range(self.N):
				self.p[k] += self.check_walls(i, self.L)
				self.r[i] += self.dt * self.v[i]
				modv = modv + LA.norm(self.v[i])

			self.CM = np.sum(self.r, axis=0) / self.N
			self.CM_t.append(self.CM)

			average_v = modv / self.N
			self.hist.append(average_v)

			lines[0].set_data(self.r[:,0], self.r[:,1])
			lines[0].set_3d_properties(self.r[:,2])

			lines[1].set_data([self.CM[0]], [self.CM[1]])
			lines[1].set_3d_properties([self.CM[2]])

			return lines

		# Attaching 3D axis to the figure
		fig = plt.figure()
		ax = p3.Axes3D(fig)


		# Setting the axes properties
		ax.set_xlim3d([-self.L/2, self.L/2])
		ax.set_xlabel('X')

		ax.set_ylim3d([-self.L/2, self.L/2])
		ax.set_ylabel('Y')

		ax.set_zlim3d([-self.L/2, self.L/2])
		ax.set_zlabel('Z')

		T1 = calculate_T(self.v, self.m, self.N)

		lines = []
		lines.append(ax.plot(self.r[:,0], self.r[:,1], self.r[:,2], ls='None', marker='.', label = 'Particles at $T_0=%.2f\,K$'%T1)[0])
		lines.append(ax.plot([self.CM[0]], [self.CM[1]], [self.CM[2]], marker='o', color='r', label='Center of masses')[0])

		ani = FuncAnimation(fig, update, fargs=(lines,), frames=np.linspace(0, self.T_MAX-self.dt, self.Nt),
		                    blit=True, interval=1, repeat=False)
		#ani.save('animation.mp4', fps=20, writer="ffmpeg", codec="libx264")
		plt.legend(loc = 'upper left')
		plt.show()

	def experiments(self):

		Exp_P, Theo_P = calculate_P(self.N, self.L, self.T, self.dt, self.p)
		T1 = calculate_T(self.v, self.m, self.N)

		A = 6*self.L**2

		theo_P_2 = self.N * k_B * T1 / (self.L**3)

		rel_error = abs(Theo_P - Exp_P) / Theo_P * 100
		rel_error2 = abs(theo_P_2 - Exp_P) / theo_P_2 * 100

		#PRESSURE PLOT
		theo_array = []
		exp_array = []

		for i in range(self.Nt):
			theo_array.append(Theo_P)
			exp_array.append(Exp_P)

		plt.plot(np.linspace(0, self.T_MAX - self.dt, self.Nt), self.p / (self.dt * A), ':', label = 'Experimental pressure')
		plt.plot(np.linspace(0, self.T_MAX - self.dt, self.Nt), theo_array, label = 'Theoretical pressure')
		plt.plot(np.linspace(0, self.T_MAX - self.dt, self.Nt), exp_array, label = 'Average experimental pressure')
		plt.title('Pressure of the system')
		plt.xlabel('t(s)')
		plt.ylabel('P(Pa)')
		plt.legend()
		
		os.chdir(carpeta_data)
		plt.savefig('Pressure.png')
		plt.gcf().clear()

		#VELOCITY PLOT
		theo_array = []
		theo_array_rms = []

		vmed = np.sqrt(8 * k_B * self.T / (np.pi * self.m))
		vrms = np.sqrt(3 * k_B * self.T / self.m)

		exp_avg_vel = np.sum(self.hist) / len(self.hist)

		rel_error_vel = abs(vmed-exp_avg_vel)/vmed * 100

		fill_database('gas_database', 'real_gas', '(?,?,?,?,?,?,?)', [T1, Theo_P, Exp_P, rel_error, vmed, exp_avg_vel, rel_error_vel])

		plt.plot(np.linspace(0, self.T_MAX - self.dt, len(self.hist)), self.hist, '-', label = 'Experimental data')
		plt.plot(np.linspace(0, self.T_MAX - self.dt, self.Nt), np.linspace(vmed, vmed, self.Nt), label = r'Maxwell-Boltzman $v_{med}$')
		plt.plot(np.linspace(0, self.T_MAX - self.dt, self.Nt), np.linspace(vrms, vrms, self.Nt), label = r'Maxwell-Boltzman $v_{rms}$')
		plt.title('Average velocity ($T=%.2f\,K$)' % T1)
		plt.xlabel('t(s)')
		plt.ylabel('$v_{med}$(m/s)')
		plt.legend()
		
		os.chdir(carpeta_data)
		plt.savefig('Velocity.png')
		plt.gcf().clear()

		#VELOCITY HISTOGRAM

		modv = []
		for i in range(len(self.v)):
			mod = self.v[i][0]**2 + self.v[i][1]**2 + self.v[i][2]**2
			modv.append(np.sqrt(mod))

		v_max = 3*np.sqrt((3*k_B* self.T)/ self.m)
		dv = 0.01

		vel_hist(modv, self.v, v_max, dv, self.m, self.T)

		#CENTER OF MASSES MOVEMENT PLOTS

		transposed = np.transpose(self.CM_t)
		t = np.linspace(0, self.T_MAX - self.dt, len(self.CM_t))

		std_x = np.std(transposed[:][0])
		std_y = np.std(transposed[:][1])
		std_z = np.std(transposed[:][2])

		mean_x = np.mean(transposed[:][0])
		mean_y = np.mean(transposed[:][1])
		mean_z = np.mean(transposed[:][2])

		plt.subplot(2, 2, 1)
		plt.plot(t, transposed[:][0]) #CM_x
		meanx = plt.plot(t, np.linspace(mean_x, mean_x, len(self.CM_t)), color='red', label='Mean') #mean
		sdx = plt.plot(t, np.linspace(std_x, std_x, len(self.CM_t)), ls = ':', color='green', label=r'$\pm$ SD') #+ standard deviation
		plt.plot(t, np.linspace(-std_x, -std_x, len(self.CM_t)), ls = ':', color='green') #- standard deviation
		plt.xlabel('t')
		plt.ylabel('x')

		plt.subplot(2, 2, 2)
		plt.plot(t, transposed[:][1]) #CM_y
		plt.plot(t, np.linspace(mean_y, mean_y, len(self.CM_t)), color='red', label='Mean') #mean
		plt.plot(t, np.linspace(std_y, std_y, len(self.CM_t)), ls = ':', color='green', label=r'$\pm$ SD') #+ standard deviation
		plt.plot(t, np.linspace(-std_y, -std_y, len(self.CM_t)), ls = ':', color='green') #- standard deviation
		plt.xlabel('t')
		plt.ylabel('y')

		plt.subplot(2, 2, 3)
		plt.plot(t, transposed[:][2]) #CM_z
		plt.plot(t, np.linspace(mean_z, mean_z, len(self.CM_t)), color='red', label='Mean') #mean
		plt.plot(t, np.linspace(std_z, std_z, len(self.CM_t)), ls = ':', color='green', label=r'$\pm$ SD') #+ standard deviation
		plt.plot(t, np.linspace(-std_z, -std_z, len(self.CM_t)), ls = ':', color='green') #- standard deviation
		plt.xlabel('t')
		plt.ylabel('z')

		plt.suptitle(r'$\vec{r}_{CM}$ vs t')
		plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

		os.chdir(carpeta_data)
		plt.savefig('CM_plots.png')
		plt.gcf().clear()


#Fer isothermal and isochoric processes a la classe de la interfase per utilitzar la barra de carrega
