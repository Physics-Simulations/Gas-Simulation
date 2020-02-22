import sys, re
import PyQt5
from PyQt5.QtWidgets import *
from PyQt5 import uic
from PyQt5.QtCore import pyqtSlot, QDate, Qt
from PyQt5.QtGui import QIcon, QPixmap, QFont, QImage

import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.animation import FuncAnimation
import matplotlib.colors
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.patches as mpatches
from scipy import stats

import time
import os
from shutil import copyfile

from Gas_classes import *

dir_principal = os.getcwd()
carpeta_data = dir_principal + '\Data'

if not os.path.exists(carpeta_data): os.mkdir(carpeta_data)

os.chdir(carpeta_data)
create_database()

class Window(QMainWindow): 
	def __init__(self):
		QMainWindow.__init__(self)
		os.chdir(carpeta_data)
		uic.loadUi('main.ui', self)
		os.chdir(dir_principal)

		self.showMaximized()

		#Animation and simulations
		self.simulation.clicked.connect(self.fer_simu)
		self.animation.clicked.connect(self.fer_anim)

		#Save files
		self.save_files = Save_files()

		self.save_pressure.clicked.connect(self.obrir_save_files_pressure)
		self.save_cm.clicked.connect(self.obrir_save_files_cm)
		self.save_velocity_plot.clicked.connect(self.obrir_save_files_velocity_plot)
		self.save_velocity_hist.clicked.connect(self.obrir_save_files_velocity_hist)

		self.save_pv_plot.clicked.connect(self.obrir_save_files_pv_plot)
		self.save_expansion_V.clicked.connect(self.obrir_save_files_expansion_V)
		self.save_expansion_t.clicked.connect(self.obrir_save_files_expansion_time)

		self.save_pt_plot.clicked.connect(self.obrir_save_files_pt_plot)
		self.save_heating_T.clicked.connect(self.obrir_save_files_heating_T)
		self.save_heating_t.clicked.connect(self.obrir_save_files_heating_t)

		self.save_kt.clicked.connect(self.obrir_save_files_kt)

		#Plots
		self.pressure.clicked.connect(self.veure_grafic_pressure)
		self.cm.clicked.connect(self.veure_grafic_CM)
		self.velocity_plot.clicked.connect(self.veure_grafic_velocity)
		self.velocity_hist.clicked.connect(self.veure_grafic_velocity_hist)

		self.pv_plot.clicked.connect(self.veure_grafic_pv_diagram)
		self.expansion_V.clicked.connect(self.veure_grafic_expansion_volume)
		self.expansion_t.clicked.connect(self.veure_grafic_expansion_time)

		self.pt_plot.clicked.connect(self.veure_grafic_pt_diagram)
		self.heating_T.clicked.connect(self.veure_grafic_heating_temperature)
		self.heating_t.clicked.connect(self.veure_grafic_heating_time)

		self.k_T_plot.clicked.connect(self.veure_grafic_kT)

		#L evolution
		self.pv_diagram.clicked.connect(self.P_V_diagram)
		self.expansion.clicked.connect(self.isothermal_expansion)

		#T evolution
		self.pt_diagram.clicked.connect(self.P_T_diagram)
		self.heating.clicked.connect(self.isochoric_heating)

		#kT(T)
		self.calculate_kt.clicked.connect(self.kT_T)


	def obrir_save_files_pressure(self):
		self.save_files.saveFileDialog('Pressure')

	def obrir_save_files_cm(self):
		self.save_files.saveFileDialog('CM_plots')

	def obrir_save_files_velocity_plot(self):
		self.save_files.saveFileDialog('Velocity')

	def obrir_save_files_velocity_hist(self):
		self.save_files.saveFileDialog('Velocity_hist')

	def obrir_save_files_pv_plot(self):
		self.save_files.saveFileDialog('PV_diagram')

	def obrir_save_files_expansion_V(self):
		self.save_files.saveFileDialog('Expansion_volume')

	def obrir_save_files_expansion_time(self):
		self.save_files.saveFileDialog('Expansion_time')

	def obrir_save_files_pt_plot(self):
		self.save_files.saveFileDialog('PT_diagram')

	def obrir_save_files_heating_T(self):
		self.save_files.saveFileDialog('Heating_temperature')

	def obrir_save_files_heating_t(self):
		self.save_files.saveFileDialog('Heating_time')

	def obrir_save_files_kt(self):
		self.save_files.saveFileDialog('k_T(T)')


	def veure_grafic_pressure(self):
		self.imageLabel.clear()
		filename = 'Pressure.png'

		os.chdir(carpeta_data)

		if os.path.exists(filename):
			image = QImage(filename)

			self.imageLabel.setPixmap(QPixmap.fromImage(image))
		else:
			QMessageBox.warning(self, 'Warning!', 'You must first simulate the system!')

	def veure_grafic_CM(self):
		self.imageLabel.clear()
		filename = 'CM_plots.png'

		os.chdir(carpeta_data)

		if os.path.exists(filename):
			image = QImage(filename)

			self.imageLabel.setPixmap(QPixmap.fromImage(image))
		else:
			QMessageBox.warning(self, 'Warning!', 'You must first simulate the system!')

	def veure_grafic_velocity(self):
		self.imageLabel.clear()
		filename = 'Velocity.png'

		os.chdir(carpeta_data)

		if os.path.exists(filename):
			image = QImage(filename)

			self.imageLabel.setPixmap(QPixmap.fromImage(image))
		else:
			QMessageBox.warning(self, 'Warning!', 'You must first simulate the real gas system!')

	def veure_grafic_velocity_hist(self):
		self.imageLabel.clear()
		filename = 'Velocity_hist.png'

		os.chdir(carpeta_data)

		if os.path.exists(filename):
			image = QImage(filename)

			self.imageLabel.setPixmap(QPixmap.fromImage(image))
		else:
			QMessageBox.warning(self, 'Warning!', 'You must first simulate the real gas system!')

	def veure_grafic_pv_diagram(self):
		self.imageLabel.clear()
		filename = 'PV_diagram.png'

		os.chdir(carpeta_data)

		if os.path.exists(filename):
			image = QImage(filename)

			self.imageLabel.setPixmap(QPixmap.fromImage(image))
		else:
			QMessageBox.warning(self, 'Warning!', 'You must first simulate the system!')

	def veure_grafic_expansion_volume(self):
		self.imageLabel.clear()
		filename = 'Expansion_volume.png'

		os.chdir(carpeta_data)

		if os.path.exists(filename):
			image = QImage(filename)

			self.imageLabel.setPixmap(QPixmap.fromImage(image))
		else:
			QMessageBox.warning(self, 'Warning!', 'You must first simulate the system!')

	def veure_grafic_expansion_time(self):
		self.imageLabel.clear()
		filename = 'Expansion_time.png'

		os.chdir(carpeta_data)

		if os.path.exists(filename):
			image = QImage(filename)

			self.imageLabel.setPixmap(QPixmap.fromImage(image))
		else:
			QMessageBox.warning(self, 'Warning!', 'You must first simulate the system!')

	def veure_grafic_pt_diagram(self):
		self.imageLabel.clear()
		filename = 'PT_diagram.png'

		os.chdir(carpeta_data)

		if os.path.exists(filename):
			image = QImage(filename)

			self.imageLabel.setPixmap(QPixmap.fromImage(image))
		else:
			QMessageBox.warning(self, 'Warning!', 'You must first simulate the system!')

	def veure_grafic_heating_temperature(self):
		self.imageLabel.clear()
		filename = 'Heating_temperature.png'

		os.chdir(carpeta_data)

		if os.path.exists(filename):
			image = QImage(filename)

			self.imageLabel.setPixmap(QPixmap.fromImage(image))
		else:
			QMessageBox.warning(self, 'Warning!', 'You must first simulate the system!')

	def veure_grafic_heating_time(self):
		self.imageLabel.clear()
		filename = 'Heating_time.png'

		os.chdir(carpeta_data)

		if os.path.exists(filename):
			image = QImage(filename)

			self.imageLabel.setPixmap(QPixmap.fromImage(image))
		else:
			QMessageBox.warning(self, 'Warning!', 'You must first simulate the system!')

	def veure_grafic_kT(self):
		self.imageLabel.clear()
		filename = 'K_T(T).png'

		os.chdir(carpeta_data)

		if os.path.exists(filename):
			image = QImage(filename)

			self.imageLabel.setPixmap(QPixmap.fromImage(image))
		else:
			QMessageBox.warning(self, 'Warning!', 'You must first simulate the system!')


	def fer_simu(self):
		t_0 = current_milli_time()

		N = self.N.value()
		L = self.L.value() 
		T = self.T.value()
		m = float(self.m.text())
		R = float(self.R.text())
		dt = float(self.dt.text())
		T_MAX = self.t_max.value()

		if self.gas_type_sim.currentText() == "Ideal gas":
			gas = IdealGas(N, L, T, m, R, dt, T_MAX)
			gas.simulation()
			gas.experiments()
	
		else:
			gas = RealGas(N, L, T, m, R, dt, T_MAX)
			gas.simulation()
			gas.experiments()

		time_used = (current_milli_time() - t_0) / 1000

		self.execution_time.setText(str(time_used) + ' s')
		QMessageBox.information(self, 'Information', 'Simulation finished!')

	def fer_anim(self):
		t_0 = current_milli_time()

		N = self.N.value()
		L = self.L.value() 
		T = self.T.value()
		m = float(self.m.text())
		R = float(self.R.text())
		dt = float(self.dt.text())
		T_MAX = self.t_max.value()

		if self.gas_type_sim.currentText() == "Ideal gas":
			gas = IdealGas(N, L, T, m, R, dt, T_MAX)
			gas.animation()
			gas.experiments()
	
		else:
			gas = RealGas(N, L, T, m, R, dt, T_MAX)
			gas.animation()
			gas.experiments()

		time_used = (current_milli_time() - t_0) / 1000

		self.execution_time.setText(str(time_used) + ' s')

	def P_V_diagram(self):
		t_0 = current_milli_time()

		N = self.N.value()
		L = self.L.value() 
		T = self.T.value()
		m = float(self.m.text())
		R = float(self.R.text())
		dt = float(self.dt.text())
		T_MAX = self.t_max.value()

		initial_L = self.initial_L.value()
		final_L = self.final_L.value()
		n = self.points_L.value()
		pull_period = self.pull_period.value()

		Ls = np.linspace(initial_L, final_L, n)

		exp_pressures = []
		theo_pressures = []

		if self.isothermal_process.currentText() == "Ideal gas":

			i = 0
			total = len(Ls)

			for L in Ls:
				i += 1

				gas = IdealGas(N, L, T, m, R, dt, T_MAX)
				p = gas.simulation()

				(exp, theo) = calculate_P(N, L, T, dt, p)

				exp_pressures.append(exp)
				theo_pressures.append(theo)

				self.pv_bar.setValue(int(i / total * 100))

			Vs = []
			for i in range(len(Ls)):
				Vs.append(Ls[i]**3)

			plt.plot(Vs, theo_pressures, '-o', label = 'Ideal gas pressure')
			plt.plot(Vs, exp_pressures, '-o', label = 'Experimental pressure')
			plt.title('Isothermal compression (T=%.2f K)'%T)
			plt.xlabel('Volume ($m^3$)')
			plt.ylabel('Pressure (Pa)')
			plt.legend()
			
			os.chdir(carpeta_data)
			plt.savefig('PV_diagram.png')
			plt.gcf().clear()

			time_used = (current_milli_time() - t_0) / 1000
			self.execution_time.setText(str(time_used) + ' s')

			QMessageBox.information(self, 'Information', 'Simulation finished!')

			self.pv_bar.setValue(0)
	
		else:

			i = 0
			total = len(Ls)

			for L in Ls:
				i += 1

				gas = RealGas(N, L, T, m, R, dt, T_MAX)
				p = gas.simulation()

				(exp, theo) = calculate_P(N, L, T, dt, p)

				exp_pressures.append(exp)
				theo_pressures.append(theo)

				self.pv_bar.setValue(int(i / total * 100))

			Vs = []
			for i in range(len(Ls)):
				Vs.append(Ls[i]**3)

			plt.plot(Vs, theo_pressures, '-o', label = 'Ideal gas pressure')
			plt.plot(Vs, exp_pressures, '-o', label = 'Experimental pressure')
			plt.title('Isothermal compression (T=%.2f K)'%T)
			plt.xlabel('Volume ($m^3$)')
			plt.ylabel('Pressure (Pa)')
			plt.legend()
			
			os.chdir(carpeta_data)
			plt.savefig('PV_diagram.png')
			plt.gcf().clear()

			time_used = (current_milli_time() - t_0) / 1000
			self.execution_time.setText(str(time_used) + ' s')
		
			QMessageBox.information(self, 'Information', 'Simulation finished!')
			self.pv_bar.setValue(0)

	def isothermal_expansion(self):
		t_0 = current_milli_time()

		N = self.N.value()
		L = self.L.value() 
		T = self.T.value()
		m = float(self.m.text())
		R = float(self.R.text())
		dt = float(self.dt.text())
		T_MAX = self.t_max.value()

		Lf = self.final_L.value()
		L0 = self.initial_L.value()
		pull_period = self.pull_period.value()

		dL = pull_period * (Lf-L0)/T_MAX

		if pull_period > T_MAX :
			QMessageBox.warning(self, 'Warning', 'Pull period cant\'t be greater than simulation time!')

		else:
			if self.isothermal_process.currentText() == "Ideal gas":
				gas = IdealGas(N, L, T, m, R, dt, T_MAX)
				gas.isothermal_expansion(pull_period, L0, dL)
			else:
				gas = RealGas(N, L, T, m, R, dt, T_MAX)
				gas.isothermal_expansion(pull_period, L0, dL)

			time_used = (current_milli_time() - t_0) / 1000
			self.execution_time.setText(str(time_used) + ' s')

			QMessageBox.information(self, 'Information', 'Simulation finished!')

	def P_T_diagram(self):
		t_0 = current_milli_time()

		N = self.N.value()
		L = self.L.value() 
		T = self.T.value()
		m = float(self.m.text())
		R = float(self.R.text())
		dt = float(self.dt.text())
		T_MAX = self.t_max.value()

		initial_T = self.initial_T.value()
		final_T = self.final_T.value()
		n = self.points_T.value()
		heat_period = self.heat_period.value()

		Ts = np.linspace(initial_T, final_T, n)

		exp_pressures = []
		theo_pressures = []

		if self.isochoric_process.currentText() == "Ideal gas":
			
			i = 0
			total = len(Ts)

			for T in Ts:
				i += 1

				gas = IdealGas(N, L, T, m, R, dt, T_MAX)
				p = gas.simulation()

				(exp, theo) = calculate_P(N, L, T, dt, p)

				exp_pressures.append(exp)
				theo_pressures.append(theo)

				self.pt_bar.setValue(int(i / total * 100))

			plt.plot(Ts, theo_pressures, '-o', label = 'Ideal gas pressure')
			plt.plot(Ts, exp_pressures, '-o', label = 'Experimental pressure')
			plt.title(r'Isochoric heating (V=%.2f $m^3$)'% L**3)
			plt.xlabel(r'Volume ($m^3$)')
			plt.ylabel('Pressure (Pa)')
			plt.legend()
			
			os.chdir(carpeta_data)
			plt.savefig('PT_diagram.png')
			plt.gcf().clear()

			time_used = (current_milli_time() - t_0) / 1000
			self.execution_time.setText(str(time_used) + ' s')

			QMessageBox.information(self, 'Information', 'Simulation finished!')

			self.pt_bar.setValue(0)
	
		else:

			i = 0
			total = len(Ts)

			for T in Ts:
				i += 1

				gas = RealGas(N, L, T, m, R, dt, T_MAX)
				p = gas.simulation()

				(exp, theo) = calculate_P(N, L, T, dt, p)

				exp_pressures.append(exp)
				theo_pressures.append(theo)

				self.pt_bar.setValue(int(i / total * 100))

			plt.plot(Ts, theo_pressures, '-o', label = 'Ideal gas pressure')
			plt.plot(Ts, exp_pressures, '-o', label = 'Experimental pressure')
			plt.title(r'Isothermal compression (V=%.2f $m^3$)' % L**3)
			plt.xlabel(r'Volume ($m^3$)')
			plt.ylabel('Pressure (Pa)')
			plt.legend()
			
			os.chdir(carpeta_data)
			plt.savefig('PT_diagram.png')
			plt.gcf().clear()

			time_used = (current_milli_time() - t_0) / 1000
			self.execution_time.setText(str(time_used) + ' s')
		
			QMessageBox.information(self, 'Information', 'Simulation finished!')
			self.pt_bar.setValue(0)

	def isochoric_heating(self):
		t_0 = current_milli_time()

		N = self.N.value()
		L = self.L.value() 
		T = self.T.value()
		m = float(self.m.text())
		R = float(self.R.text())
		dt = float(self.dt.text())
		T_MAX = self.t_max.value()

		Tf = self.final_T.value()
		T0 = self.initial_T.value()
		heat_period = self.heat_period.value()

		dT = heat_period * (Tf-T0)/T_MAX

		if heat_period > T_MAX :
			QMessageBox.warning(self, 'Warning', 'Pull period cant\'t be greater than simulation time!')

		else:
			if self.isochoric_process.currentText() == "Ideal gas":
				gas = IdealGas(N, L, T, m, R, dt, T_MAX)
				gas.isochoric_heating(heat_period, T0, dT)
			else:
				gas = RealGas(N, L, T, m, R, dt, T_MAX)
				gas.isochoric_heating(heat_period, T0, dT)

			time_used = (current_milli_time() - t_0) / 1000
			self.execution_time.setText(str(time_used) + ' s')

			QMessageBox.information(self, 'Information', 'Simulation finished!')

	def kT_T(self):

		t_0 = current_milli_time()

		N = self.N.value()
		L = self.L.value() 
		T = self.T.value()
		m = float(self.m.text())
		R = float(self.R.text())
		dt = float(self.dt.text())
		T_MAX = self.t_max.value()

		initial_L = self.initial_L.value()
		final_L = self.final_L.value()
		n_L = self.points_L.value()

		initial_T = self.initial_T.value()
		final_T = self.final_T.value()
		n_T = self.points_T.value()

		Ts = np.linspace(initial_T, final_T, n_T)
		Ls = np.linspace(initial_L, final_L, n_L)

		COLORS = ['r', 'b', 'g', 'y', 'm']

		k_Ts = []
		pressures = []

		total = len(Ts) * len(Ls)
		counter = 0

		for T in Ts:

			Exp_pr = []
			Theo_pr = []

			for L in Ls:
				gas = IdealGas(N, L, T, m, R, dt, T_MAX)
				gas.init_particles()

				p = gas.simulation()

				(exp, theo) = calculate_P(N, L, T, dt, p)

				Exp_pr.append(exp)
				Theo_pr.append(theo)

				counter += 1
				self.kt_bar.setValue(int(counter / total * 100))

			Vs = []
			for i in range(len(Ls)):
				Vs.append(Ls[i]**3)

			#Calculate k_T
			k_T = []
			k_T_theo = []
			real_pr = []

			for i in range(1, len(Ls)):
				k_T.append((-1/Vs[i])*((Vs[i]-Vs[i-1])/(Exp_pr[i]-Exp_pr[i-1])))
				real_pr.append(Exp_pr[i-1])
				k_T_theo.append(1/Theo_pr[i-1])

			#print('k_T experimental at T=%i:' % T, k_T)
			#print('k_T teòrica:', k_T_theo)

			k_Ts.append(k_T)
			pressures.append(real_pr)

		for i in range(len(Ts)):
			plt.plot(k_Ts[i], pressures[i], color = COLORS[i % len(COLORS)], marker = 's', label = 'T=%i' % Ts[i])

		plt.title(r'$k_T(T)$')
		plt.xlabel(r'P ($Pa$)')
		plt.ylabel(r'$k_T$ ($Pa^{-1}$)')
		plt.legend()
		
		os.chdir(carpeta_data)
		plt.savefig('K_T(T).png')

		time_used = (current_milli_time() - t_0) / 1000
		self.execution_time.setText(str(time_used) + ' s')

		QMessageBox.information(self, 'Information', 'Simulation finished!')
		self.kt_bar.setValue(0)

	def results(self):
		pass
		
	def closeEvent(self, event):
		os.chdir(carpeta_data)
		result = QMessageBox.question(self, 'Leaving...','Do you want to exit?', QMessageBox.Yes | QMessageBox.No)
		if result == QMessageBox.Yes:
			event.accept()

			if os.path.exists('Pressure.png'): os.remove('Pressure.png')
			if os.path.exists('CM_plots.png'): os.remove('CM_plots.png')
			if os.path.exists('Velocity.png'): os.remove('Velocity.png')
			if os.path.exists('Velocity_hist.png'): os.remove('Velocity_hist.png')
			if os.path.exists('PT_diagram.png'): os.remove('PT_diagram.png')
			if os.path.exists('PV_diagram.png'): os.remove('PV_diagram.png')
			if os.path.exists('Expansion_time.png'): os.remove('Expansion_time.png')
			if os.path.exists('Expansion_volume.png'): os.remove('Expansion_volume.png')
			if os.path.exists('Heating_temperature.png'): os.remove('Heating_temperature.png')
			if os.path.exists('Heating_time.png'): os.remove('Heating_time.png')
			if os.path.exists('k_T(T).png'): os.remove('k_T(T).png')
			
		else:event.ignore()

class Save_files(QFileDialog):
	def __init__(self):
		QFileDialog.__init__(self)

		self.title = 'Save files'
		self.left = 10
		self.top = 10
		self.width = 640
		self.height = 400 

		self.initUI()

	def initUI(self):
		self.setWindowTitle(self.title)
		self.setGeometry(self.left, self.top, self.width, self.height)

	def saveFileDialog(self, name):
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog

		fileName, _ = QFileDialog.getSaveFileName(self, 'Save files') 

		if fileName:
			os.chdir(carpeta_data)
			if os.path.exists('%s.png' % name): copyfile('%s.png' % name, fileName + '.png')
			else: QMessageBox.warning(self, 'Warning!', 'The plot doesn\'t exist!') 

app = QApplication(sys.argv)
_window=Window()
_window.show()
app.exec_()