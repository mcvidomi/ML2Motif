'''
Created on 16.06.2015

@author: marinavidovic
'''
from __future__ import division
import numpy as np
import math
import pdb
import matplotlib
matplotlib.use('Agg')
import pylab as py
import itertools

def test():

    plt.plot([1,2,3])
    plt.savefig('testBB.png')

size = 44
def plot_gPOIM(poim, savefile="gPOIM.pdf", show=False):
    py.close("all")
    matplotlib.use('Agg')
    py.figure(figsize=(14, 12))
    py.ioff()
    motivelen = int(np.log(len(poim)) / np.log(4))
    ylabel = []
    for i in range(int(math.pow(4, motivelen))):
        label = []
        index = i
        for j in range(motivelen):
            label.append(index % 4)
            index = int(index / 4)
        label.reverse()
        ylabel.append(veclisttodna(label))
        
    py.pcolor(poim)
    
    cb=py.colorbar()
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(size)  
    py.axis([0, len(poim[0]), 0, len(poim)])
    diff = int((len(poim[0]) / 5)) - 1
   # x_places = py.arange(4.5, len(poim[0]), diff)
   # py.xticks(x_places, np.arange(5, len(poim[0]) + 1, diff),fontsize=40) 
    py.xlabel("Position", fontsize=size)
    py.ylabel("Motif", fontsize=size) 
    py.xticks(fontsize=size)
    py.yticks(np.arange(math.pow(4, motivelen)) + 0.5, (ylabel), fontsize=size)   
    if savefile != "":
#        F = py.gcf()
#        # Now check everything with the defaults:
#        DefaultSize = F.get_size_inches()
#        # Now make the image twice as big, while keeping the fonts and all the
#        # same size
#        F.set_size_inches((DefaultSize[0]*2, DefaultSize[1]*2))
        py.savefig(savefile,  dpi=100) 
    
    if show:
        py.show()


def figurepoimsimple_base(poim, savefile, bases, small_k):
    py.close("all")
   
    matplotlib.use('Agg')
    #py.figure(figsize=(14, 12))
    bl = len(bases)
    py.ioff()
    #ylabel = [''.join(p) for p in itertools.product(bases, repeat=small_k)]
    ylabel = [list(p) for p in itertools.product(bases, repeat=small_k)]
    
    py.pcolor(poim)
    
    cb=py.colorbar()
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(40)  
    py.axis([0, len(poim[0]), 0, len(poim)])
    diff = int((len(poim[0]) / 5)) - 1
   # x_places = py.arange(4.5, len(poim[0]), diff)
   # py.xticks(x_places, np.arange(5, len(poim[0]) + 1, diff),fontsize=40) 
    py.xlabel("Position", fontsize=46)
    py.ylabel("Motif", fontsize=16)
    py.yticks(np.arange(len(ylabel)) + 0.5, (ylabel), fontsize=12)   
    if savefile != "":
#        F = py.gcf()
#        # Now check everything with the defaults:
#        DefaultSize = F.get_size_inches()
#        # Now make the image twice as big, while keeping the fonts and all the
#        # same size
#        F.set_size_inches((DefaultSize[0]*2, DefaultSize[1]*2))
        py.savefig(savefile) 
    py.close("all")
    
    
     
def figurepoimsimple_small(poim, l, start, savefile, show):
    R = poim

    py.figure(figsize=(14, 12))
    motivelen = int(np.log(len(poim)) / np.log(4))
    ylabel = []
    for i in range(int(math.pow(4, motivelen))):
        label = []
        index = i
        for j in range(motivelen):
            label.append(index % 4)
            index = int(index / 4)
        label.reverse()
        ylabel.append(veclisttodna(label))
    py.pcolor(R[:, start:start + l])
    cb=py.colorbar()
    for t in cb.ax.get_yticklabels():
     t.set_fontsize(40)   
    diff = int((l / 5)) - 1
    x_places = py.arange(0.5, l, diff)
    xa = np.arange(start, start + l, diff)
    diff = int((l / 4)) 
    x_places = py.arange(0.5, l , diff)
    xa = np.arange(start + 1, start + 1 + l, diff)
    py.xlabel("Position", fontsize=46)
    py.ylabel("Motif", fontsize=46)
    py.yticks(np.arange(math.pow(4, motivelen)) + 0.5, (ylabel),fontsize=40)   
    py.xticks(x_places, (xa.tolist()),fontsize=40)
    if savefile != "":
        py.savefig(savefile)  
    print "the poim should show up here"
    if show:
        py.show()


def figuremaxPOIM_small(poim, savepath, l , start, show):  
    R = poim

    py.figure(figsize=(14, 12))
    py.pcolor(R[:, start:start + l])
    cb=py.colorbar()
    for t in cb.ax.get_yticklabels():
     t.set_fontsize(40)
    diff = int((l / 4)) 
    x_places = py.arange(0.5, l , diff)
    xa = np.arange(start + 1, start + 1 + l, diff)
    py.xlabel("Position", fontsize=46)
    py.ylabel("Motif length", fontsize=46)
    py.xticks(x_places, (xa.tolist()),fontsize=40)
    py.yticks(np.arange(0, len(poim)) + 0.5, np.arange(1, len(poim) + 1),fontsize=40)    
    py.savefig(savepath)
    if show:
        py.show()

  
def figuremaxPOIM(R, savepath, show):  

    py.figure(figsize=(14, 12))
    py.pcolor(R)
    cb=py.colorbar()
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(size)
    py.axis([0, len(R[0]), 0, len(R)])
    diff = 5#int((len(R[0]) / 5)) - 1
    x_places = py.arange(4.5, len(R[0]), diff)
    py.xticks(x_places, np.arange(5, len(R[0]) + 1, diff),fontsize=size) 
    py.xlabel("Position", fontsize=size)
    py.ylabel("Motif length", fontsize=size)
    py.yticks(np.arange(0, len(R)) + 0.5, np.arange(1, len(R) + 1),fontsize=size)    
    py.savefig(savepath, dpi=100)
    if show:
        py.show()


def figure_column_sum_poim(poim,savepath):
    py.close("all")
    maxpoim = np.sum(poim, axis=0)
    py.plot(maxpoim)
    py.plot([np.mean(maxpoim)]*len(maxpoim))
    py.xlabel("Position")
    py.ylabel("Column-wise sum")
    py.savefig(savepath)
    
def veclisttodna(vec):
    """
    converts the digits to the dna alphabet and returns dna string
    """
    str = ""
    for i in range(len(vec)):
        if vec[i] == 0:
            str += 'A'
        elif vec[i] == 1:
            str += ('C')
        elif vec[i] == 2:
            str += ('G')
        else:
            str += ('T')
    return str


# 
# from matplotlib.patches import Patch
import matplotlib.pyplot as plt
# import pylab as py
# import toolbox.addmotivepoims as add

# import numpy as np
# import math
def vectodna(vec):
    """
    converts the digits to the dna alphabet and returns dna string
    """
    str = []
    for i in range(len(vec)):
        if vec[i] == 0:
            str.append('A')
        elif vec[i] == 1:
            str.append('C')
        elif vec[i] == 2:
            str.append('G')
        else:
            str.append('T')
    return str

def veclisttodna(vec):
    """
    converts the digits to the dna alphabet and returns dna string
    """
    str = ""
    for i in range(len(vec)):
        if vec[i] == 0:
            str += 'A'
        elif vec[i] == 1:
            str += ('C')
        elif vec[i] == 2:
            str += ('G')
        else:
            str += ('T')
    return str

def figurepoimcompare(R, Q, filename):
    motivelen = int(np.log(len(Q)) / np.log(4))
    plt.clf()
    ylabeling = []
    for i in range(int(math.pow(4, motivelen))):
        label = [] 
        index = i
        for j in range(motivelen):
            label.append(index % 4)
            index = int(index / 4)
        label.reverse()
        ylabeling.append(veclisttodna(label)) 
    
    fignumber = 0
    if len(np.shape(Q)) > 2:
        figquantity = np.shape(Q)[0] 
        for i in range(figquantity):
            fignumber += 1
            py.subplot(2, float(figquantity) / 2 + figquantity % 2 + 1, fignumber)
            py.pcolor(Q[i])
            diff = int((len(Q[0]) / 5)) - 1
            x_places = py.arange(0.5, len(Q[0]), diff)
            py.xticks(x_places, np.arange(0, len(Q[0]) + 1, diff)) 
            py.yticks(np.arange(math.pow(4, motivelen)) + 0.5, (ylabeling))
            py.title("Q")
            py.axis([0, len(Q[i][0]), 0, len(Q[i])])
            fignumber += 1
            py.subplot(2, float(figquantity) / 2 + figquantity % 2 + 1, fignumber)
            py.pcolor(R[i])
            py.colorbar()
            py.title("R")
            py.yticks(visible=False)
            py.xticks(x_places, np.arange(0, len(Q[0]) + 1, diff)) 
        if filename != "":    
            py.savefig(filename + ".png")

    else:
        fignumber += 1
        py.subplot(1, 2, fignumber)
        py.pcolor(Q)
        diff = int((len(Q[0]) / 5)) - 1
        x_places = py.arange(0.5, len(Q[0]), diff)
        py.xticks(x_places, np.arange(0, len(Q[0]) + 1, diff)) 
        py.yticks(np.arange(math.pow(4, motivelen)) + 0.5, (ylabeling))
        py.title("Q")
        py.subplot(1, 2, fignumber + 1)
        py.pcolor(R)
        py.title("R")
        fignumber += 1
        py.colorbar()
        py.yticks(visible=False)
        py.xticks(x_places, np.arange(0, len(Q[0]) + 1, diff)) 
        
        
        if filename != "": 
            py.savefig(filename + ".png")
        

    py.close()

def show_values(sequences, pc, fmt="%.2f", **kw):
    '''
    Heatmap with text in each cell with matplotlib's pyplot
    Source: http://stackoverflow.com/a/25074150/395857 
    By HYRY
    '''
    from itertools import izip
    pc.update_scalarmappable()
    ax = pc.get_axes()
    for p, color, value in izip(pc.get_paths(), pc.get_facecolors(), pc.get_array()):
        x, y = p.vertices[:-2, :].mean(0)
        if np.all(color[:3] > 0.5):
            color = (0.0, 0.0, 0.0)
        else:
            color = (1.0, 1.0, 1.0)
        ax.text(x, y, sequences[int(round(y)-1)][int(round(x)-1)], ha="center", va="center", color=color, **kw)

def cm2inch(*tupl):
    '''
    Specify figure size in centimeter in matplotlib
    Source: http://stackoverflow.com/a/22787457/395857
    By gns-ank
    '''
    inch = 2.54
    if type(tupl[0]) == tuple:
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

def heatmap(test_sequences,AUC, title, xlabel, ylabel, xticklabels, yticklabels,xinch=50,yinch=15):
    '''
    Inspired by:
    - http://stackoverflow.com/a/16124677/395857 
    - http://stackoverflow.com/a/25074150/395857
    '''

    # Plot it out
    fig, ax = py.subplots()    
    c = ax.pcolor(AUC, edgecolors='k', linestyle= 'solid', linewidths=0.5, cmap='coolwarm')

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(AUC.shape[0]) + 0.5, minor=False)

    ax.set_xticks(np.arange(AUC.shape[1]) + 0.5, minor=False)

    # set tick labels
    #ax.set_xticklabels(np.arange(1,AUC.shape[1]+1), minor=False)
    #ax.set_xticklabels(xticklabels, minor=False,fontsize=22)
    ax.set_yticklabels(yticklabels, minor=False,fontsize=22)

    # set title and x/y labels
  #  py.title(title)
    py.xlabel(xlabel,fontsize=22)
    py.ylabel(ylabel,fontsize=22)      

    # Remove last blank column
    py.xlim( (0, AUC.shape[1]) )

    # Turn off all the ticks
    
    ax = py.gca()
    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    # Add color bar
   # py.colorbar(c)

    # Add text in each cell 
    show_values(test_sequences,c)
    # resize 
    fig = py.gcf()
    fig.set_size_inches(cm2inch(xinch, yinch))
    
