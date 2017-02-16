import sys
sys.path.append("/home/mvidovic/POIMall/MarinaVidovic/Exploration/")
import copy
import pdb
import math
import numpy as np
import pickle

from poim import *
import tools
import view
import itertools



class POIM:
    '''
    --- use shogun toolbox ---
    '''

    def __init__(self, SVMobj=None,maxdegree=7):
        if not(SVMobj is None):
            self.SVMobj = SVMobj
            self.poims = self.SVMobj.ma[0]
            self.maxPOIM = self.SVMobj.ma[1]
            self.diffPOIM = self.SVMobj.ma[2]
        self.maxdegree = maxdegree

    def save(self,file):
        fobj = open(file, 'wb')
        pickle.dump(self.__dict__, fobj)
        fobj.close()

    def load(self,file):
        fobj = open(file, 'rb')
        self.__dict__ = pickle.load(fobj)
        fobj.close()



    def makePOIM(self, SVMobj=None, tofile = None):
        if SVMobj is None:
            SVMobj = self.SVMobj
        from poim import *
        ma = compute_poims(SVMobj.svm, SVMobj.kernel, self.maxdegree, SVMobj.L)
        if tofile != None: self.save(tofile)
        self.poims = ma[0]
        self.maxPOIM = ma[1]
        self.diffPOIM = ma[2]


    def normalize(self,Q):
        self.poim_norm = tools.normalizemotif(Q)

    def plot_poim(self,mode,Q,tofile):
        '''
        plot poim to file
        :param mode:
            1 - plot poim
            2 - plot diffPOIM
        :param Q: poim
        :param tofile: path to file
        :return: -
        '''
        if mode == 1:
            view.figurepoimsimple(Q, tofile, False)
        if mode == 2:
            view.figuremaxPOIM(Q, tofile, False)

    def poim2file(self, poim, file):
        fobj = open(file, 'wb')
        pickle.dump(poim, fobj)
        fobj.close()

    def loadpoim(self,file):
        fobj = open(file, 'rb')
        self.gPOIM = pickle.load(fobj)
        fobj.close()
        return self.gPOIM


class gPOIM:

    def get_scores(self,k_mer, position, SVMobj, xte):
        xte2 = []
        for x in xte:
            xte2.append(x[0:position] + k_mer + x[(position + len(k_mer)):])

        scores, _ = SVMobj.apply(xte2)

        return scores

    #def diffPOIM(self, max_k=3):
    #
    #    for i in range(1,k+1):



    def set_up(self, SVMobj, samplesize=100, small_k = 2, subPOIM=None):

        self.small_k = small_k
        self.samplesize = samplesize
        self.L = SVMobj.L
        # poim = conditioned_sampling_POIM(clf,samplesize,prior_probabilities,tally,x,y)
        xte = tools.gensequences_rnd(SVMobj.L, self.samplesize)
        base = ['A', 'C', 'G', 'T']

        if subPOIM is None:
            start = 0
            end = SVMobj.L - self.small_k + 1
            k = end - start
        else:
            start = subPOIM[0]
            end = subPOIM[1]


        kmers = [''.join(p) for p in itertools.product(base, repeat=self.small_k)]

        gPOIM_sub = np.zeros((int(math.pow(4, self.small_k)), end-start))
        c = -1
        for position in range(start, end):
            c += 1
            for km in range(len(kmers)):
                scores = self.get_scores(kmers[km], position, SVMobj, xte)
                gPOIM_sub[km][c] = math.pow(np.mean(scores), 1)

        self.gPOIM = np.zeros((int(math.pow(4, self.small_k)), SVMobj.L - self.small_k + 1))
        self.gPOIM[:, start:end] = gPOIM_sub

    def save(self,file=None):
        if file is None:
            file = "gPOIM.pkl"
        fobj = open(file, 'wb')
        pickle.dump(self.__dict__, fobj)
        fobj.close()

    def load(self,file=None):
        if file is None:
            file = "gPOIM.pkl"
        fobj = open(file, 'rb')
        self.__dict__ = pickle.load(fobj)
        fobj.close()

    def plot_poim(self,tofile):
        view.figurepoimsimple(self.gPOIM,tofile,show=False)









class Consensus:

    def __init__(self, exp=""):
        self.seq_list= []
        self.score_list = []
        self.category_list = []
        self.consensus_list = []
        self.exp = exp

    def add_data(self, sequences, scores, category, SVMobj, start, end, overlap, samplesize = 100):
        self.seq_list = sequences
        self.score_list = scores
        self.category_list = category
        self.start = start
        self.end = end
        self.overlap = overlap
        self.samplesize = samplesize
        self.SVMobj = SVMobj
        self.L = len(sequences[0])

        for seq in self.seq_list:
            self.consensus_list.append(self.get_consensus(seq))


    def get_consensus(self,x):
        x_rnd, y_rnd = tools.gensequences_repeatable(self.L, 0, self.samplesize, 0.0, "", 0,
                                                            np.random.randint(1, 10000)) # here seek maybe better!

        A = np.zeros(self.end-self.start-self.overlap)
        for i in range(self.start, self.end-self.overlap):
            x_rnd_copy = np.array(copy.deepcopy(x_rnd))
            # place nucleotide
            for j in range(len(x_rnd_copy)):
                x_rnd_copy[j] = x_rnd_copy[j][:i] + x[i:i+1+self.overlap] + x_rnd_copy[j][i+1+self.overlap:]
            vals,blub = self.SVMobj.apply(x_rnd_copy)
            A[i-self.start] = np.mean(vals)
        return A

    def normalize(self):

        self.consensus_list = np.array(self.consensus_list)
        self.consensus_list += np.min(self.consensus_list)
        self.consensus_list = self.consensus_list.transpose()/np.sum(self.consensus_list, axis=1)
        self.consensus_list = self.consensus_list.transpose()


    def save(self):
        f = open(self.exp + "_consensus.pkl",'wb')
        pickle.dump(self.__dict__,f)
        f.close()

    def load(self):
        f = open(self.exp + "_consensus.pkl",'rb')
        self.__dict__ = pickle.load(f)
        f.close()


    def cut_seq(self,start,stop):
        for i in range(len(self.seq_list)):
            self.seq_list[i]=self.seq_list[i][start:stop]

    def plot(self, fout = "",xinch=50,yinch=15):
        import matplotlib
        matplotlib.use('Agg')
        import pylab as py
        if type(self.consensus_list!= numpy.ndarray):
            data = np.array(self.consensus_list)
        ylabels_score = [str(np.round(p,2)) + self.category_list[i] for i,p in enumerate(self.score_list)]
        py.close("all")
        x_axis_size = self.end - self.start-self.overlap +1
        y_axis_size = len(self.seq_list)
        title = "Pixel-wise explanation"
        xlabel = "Sequence positions"
        ylabel = "Scores"

        xticklabels = range(0, x_axis_size)  # could be text
        yticklabels = ylabels_score  # could be text
        title = "test"
        view.heatmap(self.seq_list, data, title, xlabel, ylabel, xticklabels, yticklabels,xinch=xinch,yinch=yinch)
        py.xticks(range(0, self.end - self.start - self.overlap, 10), range(self.start, self.end - self.overlap, 10), fontsize=22)
        py.xlabel("Position", fontsize=22)
        py.ylabel("s(x)", fontsize=22)
        py.yticks(fontsize=22)
        if fout =="":
            py.savefig('results/consensus1.pdf',format='pdf',bbox_inches='tight')  # use format='svg' or 'pdf' for vectorial pictures
        else:
            py.savefig(fout, format='pdf', bbox_inches='tight')