import tools
import numpy as np
import pickle
import pdb
import os
import view
from SVM import *
from POIM import *
from Motif import *

# CREATE RESULT FOLDER
if not os.path.exists("results"):
    os.makedirs("results")


# READ DATA
print("LOAD/READ DATA")
xtrain,xtest,y_train,ytest = tools.read_data(job="read", lines=10000)
print("LOAD/READ DATA --- DONE!")


# TRAIN SVM
print("TRAIN SVM")

Cobj = SVM(xtrain,y_train)
Cobj.train(C=1.)
Cobj.svm_save("results/svm_trained.pkl")


# COMPUTE gPOIM
print("COMPUTE gPOIM")
small_k = 2
Pobj = gPOIM()
Pobj.set_up(Cobj, samplesize=100, small_k = small_k)
Pobj.save("results/gPOIM.pkl")


# PLOT gPOIM
view.plot_gPOIM(Pobj.gPOIM,"results/gPOIM.pdf")


# EXTRACT Motif
print("EXTRACT MOTIF")
motif_start,motif_len = get_motif_pos_len_similar_diffPOIM(Pobj.gPOIM)
motif_len=[20]
Mobj = Motif()
motif = Mobj.find(2,motif_len,Pobj.gPOIM,motif_start, base = ['A','C','G','T'])
np.savetxt("results/pwm.pkl", motif[0])
seqLogo("results/pwm.pkl")
