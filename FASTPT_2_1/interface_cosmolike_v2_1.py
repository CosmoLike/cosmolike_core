#####################
#from pylab import *
import numpy as np
#import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize as opt
import copy
#from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
from scipy.integrate import *
from scipy.special import jn
#import scipy.signal
from time import time
import subprocess as S
import argparse
import FASTPT

def call_FPT_bias(k,plin):
  #power-law bias and zero padding
  nu=-2; 
  n_pad=len(k)
  res = np.zeros((8,len(k)))
  # set the parameters for the power spectrum window and
  # Fourier coefficient window 
  P_window=np.array([.2,.2])  
  C_window=.65
  fastpt=FASTPT.FASTPT(k,nu,n_pad=n_pad, to_do=['dd_bias']) 
  bias_fpt=fastpt.one_loop_dd_bias(plin,P_window=P_window,C_window=C_window)
  res[0:7,:]  = np.asarray(bias_fpt[0:7])
  res[7,0] =bias_fpt[7]
  return np.asarray(res)
