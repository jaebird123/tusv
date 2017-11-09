#     file: jplot.py
#   author: Jesse Eaton
#  created: 11/5/2017
# modified: 11/5/2017
#  purpose: plots for tusv paper


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import numpy as np

import matplotlib.pyplot as plt
import matplotlib
import scipy
import scipy.stats
from scipy.optimize import curve_fit
from matplotlib import pyplot

import pylab as pl


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

# here


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)

	font = {'family' : 'normal',
					'size'   : 16}
	pl.rc('font', **font)

	# plot_param_tune_lambda1()
	# plot_param_tune_lambda2()
	plot_m_diff()
	

def plot_param_tune_lambda1():
	lmb1 = [0.0, 0.01, 0.05, 0.25, 1.25, 6.25]
	c_bp = [19.2, 28.0, 22.8, 6.0, 22.4, 17.4]
	c_sg = [20.0, 5.4, 5.0, 1.8, 3.4, 8.0]
	c    = [39.2, 33.4, 27.8, 7.8, 25.8, 25.4]

	fig, ax = pl.subplots()
	pl.plot(lmb1, c_bp, '-s', label = '$C_{bp}$')
	pl.plot(lmb1, c_sg, '-^', label = '$C_{sg}$')
	pl.plot(lmb1, c, '-o', label = '$C_{tot}$')

	ax.set_xscale("log", nonposx='clip')

	pl.legend(loc = 'upper right')
	pl.xlabel('$\lambda_1$')
	pl.ylabel('$|C_{tru}-C_{obs}|$')

def plot_param_tune_lambda2():
	lmb2 = [0.0, 0.01, 0.05, 0.25, 1.25, 6.25]
	c_bp = [6.0, 21.2, 12.6, 11.2, 8.0, 6.0]
	c_sg = [1.8, 4.4, 2.2, 2.2, 2.0, 1.5]
	c    = [7.8, 25.6, 14.8, 13.4, 10.0, 7.5]

	fig, ax = pl.subplots()
	pl.plot(lmb2, c_bp, '-s', label = '$C_{bp}$')
	pl.plot(lmb2, c_sg, '-^', label = '$C_{sg}$')
	pl.plot(lmb2, c, '-o', label = '$C_{tot}$')

	ax.set_xscale("log", nonposx='clip')

	pl.legend(loc = 'upper right')
	pl.xlabel('$\lambda_2$')
	pl.ylabel('$|C_{tru}-C_{obs}|$')

	pl.show()

def plot_m_diff():
	m = [1, 5, 10]
	#  m == 1, m == 5, m == 10
	c_bps = [[4.0, 0.0, 8.0], [10.0, 14, 0, 0], [10, 14]]
	c_sgs = [[2.0, 0.0, 5.0], [3.0, 3, 0, 0], [3, 3]]
	cs    = [[6.0, 0.0, 13.0], [13.0, 17, 0, 0], [13, 17]]

	c_bp, c_sg, c = [], [], []
	for i in xrange(0, len(m)):
		c_bp.append(np.mean(c_bps[i]))
		c_sg.append(np.mean(c_sgs[i]))
		c.append(np.mean(cs[i]))

	fig, ax = pl.subplots()
	pl.plot(m, c_bp, '-s', label = '$C_{bp}$')
	pl.plot(m, c_sg, '-^', label = '$C_{sg}$')
	pl.plot(m, c, '-o', label = '$C_{tot}$')

	pl.legend(loc = 'upper right')
	pl.xlabel('number of samples $m$')
	pl.ylabel('$|C_{tru}-C_{obs}|$')

	pl.show()


#
#   PLOT SV LENS
#

def plot_sv_lens():
	tum_lens = np.genfromtxt('sv_lens_59_samples.txt', delimiter = ',')
	sim_lens = np.genfromtxt('sv_lens_sim_samples.txt', delimiter = ',')

	# tum_lens = tum_lens / float(len(tum_lens))
	# sim_lens = sim_lens / float(len(sim_lens))
	# print len(tum_lens)
	# exit()

	# y_tum, x_ = plt.hist(tum_lens, 100)
	# y_sim, bins, patches = plt.hist(sim_lens, 100)
	# print n
	# print bins
	# plt.show()
	# exit()

	# popt, pcov = curve_fit(exponenial_func, bins, n, p0=(1, 1e-6, 1))
	# print popt
	# print pconv
	# exit()

	# dist = getattr(scipy.stats, 'expon')
	# param = dist.fit(n)
	# pdf_fitted = dist.pdf(bins, *param[:-2], loc=param[-2], scale=param[-1])
	# plt.plot(pdf_fitted)
	# plt.show()

	# plt.hist(sim_lens, 100)
	# plt.show()

	bins = np.linspace(200000, 6 * 10 ** 7, 100)
	bins = np.linspace(0, 200000, 100)
	pyplot.hist(tum_lens, bins, alpha=0.5, label='TCGA', normed=True)
	pyplot.hist(sim_lens, bins, alpha=0.5, label='Simulated', normed=True)
	pyplot.legend(loc='upper right')
	pyplot.xlabel('SV length')
	pyplot.ylabel('Percent occurance')
	pyplot.show()

	# print np.mean(tum_lens)
	# print np.std(tum_lens)
	# print ''
	# print np.mean(sim_lens)
	# print np.std(sim_lens)
	# print tum_lens

def exponenial_func(x, a, b, c):
	return a*np.exp(-b*x)+c


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'jplot.py', description = "plots for tusv paper")
	return vars(parser.parse_args(argv))

def is_valid_file(parser, arg):
	if not os.path.exists(arg):
		parser.error('The file \"' + str(arg) + '\" could not be found.')
	else:
		return open(arg, 'r')


# # # # # # # # # # # # # # # # # # # # # # # # #
#   C A L L   T O   M A I N   F U N C T I O N   #
# # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
	main(sys.argv[1:])
