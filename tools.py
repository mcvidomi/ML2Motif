import numpy as np
import pickle
import copy
import pdb
from sklearn.metrics.pairwise import pairwise_kernels
from sklearn.svm.libsvm import decision_function
from sklearn.model_selection import train_test_split
import numpy.matlib
import os
import pdb


def extractRealData(datapath,savepath,lines):
#path = "/home/mvidovic/POIMall/data/real/"
#filename = "human_acceptor_splice_data.txt"
    data =file(datapath).readlines()[:lines]
    labels = []
    x = []
    cn=0
    for i in range(len(data)):
        labels.append(int(data[i][0:2]))
        x.append(data[i][3:-1])
        if int(data[i][0:2])==1:
            cn=cn+1
    print "numper of positive labels: " , cn
    if savepath !="":
        fobj = open(savepath,'wb')
        pickle.dump([x,labels],fobj)
        fobj.close()
    return x,labels


def read_data(job="read",lines=100000):
    if not os.path.exists("data"):
        os.makedirs("data")
    if job=="read" or not os.path.exists("data/human_splice_data.pkl"):
        # DOWNLOAD DATA
        from StringIO import StringIO
        import urllib2
        import bz2
        # Get zip file from website
        url = 'http://www.kyb.tuebingen.mpg.de/bs/people/raetsch/LSMKL/human_acceptor_splice_data.txt_00.bz2' # just an example bz2 file
        archive = StringIO()
        url_data = urllib2.urlopen(url)
        archive.write(url_data.read())
        archive.seek(0)

        # Extract the training data
        data = bz2.decompress(archive.read())
        data = data.split("\n")
        y = []
        x = []
        for i in range(np.min([lines,len(data)])):
            try:
                y.append(int(data[i][0:2]))
                x.append(data[i][3:-1])
            except:
                pass

        f = open('data/human_splice_data.pkl', 'wb')
        pickle.dump([x,y],f)
        f.close()

    else:
        f = open('data/human_splice_data.pkl', 'rb')
        x,y = pickle.load(f)
        f.close()


    xtrain, xtest, ytrain, ytest = train_test_split(x,y,test_size=0.2, train_size=0.8, random_state=15, stratify=None)

    return xtrain, xtest, ytrain, ytest


def mutate_motif(motiv,probmut):
    dna = ['A', 'C', 'G', 'T']
    mutmot = ""
    for i in range(len(motiv)):
        rnd = random.random()
        if (rnd <= probmut):
            dnashort =['A', 'C', 'G', 'T']
            dnashort.pop(dnashort.index(motiv[i]))
            mutmot += random.choice(dnashort)
        else:
            mutmot +=motiv[i]
    return mutmot


def ini_pwm(L, N, B):
    '''
    generates a random start pwm for a motive with length motivelength
    :param L: length of motifs
    :param N: number of motifs
    :param B: lenght of base, for example, if DNA is the base, than the length is 4
    :return:
    '''

    pwm = np.zeros((N, B, L))

    for m in range(N):
        for j in range(L):
            sumcol = 0
            for i in range(B):
                pwm[m][i][j] = np.random.uniform(1, 100)
                sumcol += pwm[m][i][j]

            for i in range(B):
                pwm[m][i][j] = np.round(pwm[m][i][j] / sumcol, 2)
    return pwm

def simulate_sequence(L):
    dna = ['A', 'C', 'G', 'T']
    sequence = ''
    for i in range(L):
        sequence += np.random.choice(dna)  #zufaellig Element aus dna anhaengen
    return sequence

def gensequences_rnd(L,N):
    '''
    Generate random DNA sequences
    :param L: length of sequences
    :param N: number of sequences
    :return: randomly generated sequences
    '''
    seq_list = []
    for i in range(N):
        seq = simulate_sequence(L)
        seq_list.append(seq)
    return seq_list

def gensequences_repeatable(tally,positives,sequenceno,prob,motif,mu,seed_no):
    sequences = []
    dna = ['A', 'C', 'G', 'T']
    y=np.ones(sequenceno)*(-1)
    ml = len(motif)
    np.random.seed(seed_no)
    for i in range(sequenceno):
        A = ''
        for j in range(tally):
            A += np.random.choice(dna)  #zufaellig Element aus dna anhaengen
        if i < positives:
            y[i]=1
            mut=mutate_motif(motif,prob)
            A = A.replace(A[mu:mu + ml], mut)
        sequences.append(A)
    return sequences,y


def dna2digits_binary(dna_sequence):
    seq = []
    for i in range(len(dna_sequence)):
        if (dna_sequence[i] == 'A') or (dna_sequence[i] == 'a'):
            seq.append([0, 0, 0, 1])
        elif (dna_sequence[i] == 'C') or (dna_sequence[i] == 'c'):
            seq.append([0, 0, 1, 0])
        elif (dna_sequence[i] == 'G') or (dna_sequence[i] == 'g'):
            seq.append([0, 1, 0, 0])
        else:
            seq.append([1, 0, 0, 0])

    seq = np.array(seq).flatten()
    return seq

def dna2digit(x):
    xbinary = []
    for i in range(len(x)):
        xbinary.append(dna2digits_binary(x[i]))
    return xbinary


def dna2digits(dna_sequence):
    """
    converts the digits to the dna alphabet and returns dna string
    """
    seq = []
    for i in range(len(dna_sequence)):
        if dna_sequence[i] == 'A':
            seq.append(1)
        elif dna_sequence[i] == 'C':
            seq.append(2)
        elif dna_sequence[i] == 'G':
            seq.append(3)
        else:
            seq.append(4)
    return seq

def dnalist2digits(dna_list):
    seq = []
    for i in range(len(dna_list)):
        seq.append(dna2digits(dna_list[i]))
    return seq


def normalizemotif(poim):
    min_poim = np.min(poim)
    if min_poim < 0:
        poim = (poim + abs(min_poim) + 0.0)
    poim_norm = poim / np.sum(poim, axis=0)

    return poim_norm


def startpwmrnd_base2(P, baselen):
    '''
    generates a random start pwm for a motive with length motivelength
    '''
    pwm = []
    for l in range(len(P)):
        pwm.append(np.zeros((baselen, P[l])))
    for l in range(len(P)):
        for j in range(int(P[l])):
            sumcol = 0
            for i in range(baselen):
                pwm[l][i][j] = np.random.uniform(1, 100)
                sumcol += pwm[l][i][j]

            for i in range(baselen):
                pwm[l][i][j] = np.round(pwm[l][i][j] / sumcol, 2)

    return pwm


def converter(par, P):
    if type(par[0]) == np.float64:
        r = []
        count = 0
        for p in range(len(P)):
            pwm = np.zeros((4, P[p]))
            for i in range(4):
                for j in range(P[p]):
                    pwm[i][j] = par[count]
                    count += 1

            r.append(pwm)
        return r
    return par