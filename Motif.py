import numpy as np
import math
import tools
from scipy.optimize import minimize
from scipy.optimize.lbfgsb import fmin_l_bfgs_b
import pdb
import copy

class Motif():

    def f_L2(self, x, Q, L, mu, small_k, k):
        n_items = int(len(x) / k)
        r = np.reshape(x, ((int(len(x) / k), k)))
        R = np.zeros((int(math.pow(n_items, small_k)), L))

        D = k - small_k + 1
        po = int(math.pow(n_items, small_k))
        zlcut = int(math.pow(2, n_items) - 1)  # higher small_k than 2 may change this statment
        r1 = range(small_k - 1)
        f = 0
        for y in range(po):
            for d in range(D):
                j = mu + d
                for i in r1:
                    mov_bits = n_items / 2  # (i + 1) * 2

                    # change if small_k is bigger than 1
                    zl = (y << mov_bits & zlcut) >> mov_bits
                    zr = y >> mov_bits

                    R[y][j] += r[zr][d] * r[zl][d + 1]

                f += math.pow(R[y][j] - Q[y][j], 2)

        return f


    def f_L2_2(self,x, Q, L,  mu, small_k, P):
        r = tools.converter(x, P)
        n_items = 4  # int(len(r[0]))
        # r = np.reshape(x,((int(len(x)/k),k)))

        R = np.zeros((int(math.pow(n_items, small_k)), L))

        po = int(math.pow(n_items, small_k))
        # zlcut = int(math.pow(2, small_k * 2) - 1)
        zlcut = int(math.pow(2, n_items) - 1)  # higher small_k than 2 may change this statement
        r1 = range(small_k - 1)
        f = 0
        c = 0

        # zlcut = int(math.pow(2, small_k * 2) - 1)
        # zlcut = int(math.pow(2, n_items) - 1) # higher small_k than 2 may change this statement
        r1 = np.arange(small_k)
        r1 = r1[::-1]

        for l in range(len(P)):
            D = P[l] - small_k + 1
            for y in range(po):

                for d in range(D):
                    j = mu[c] + d
                    val = 1.
                    y_copy = copy.deepcopy(y)
                    for i in r1:
                        if i > 0:
                            mov_bits = int(math.pow(2, i))
                            z = y_copy >> mov_bits
                            # print z
                            val *= r[l][z][d + small_k - 1 - i]
                            zrcut = int(math.pow(4, i) - 1)
                            y_copy = y_copy & zrcut
                        else:
                            z = y_copy
                            # print z
                            val *= r[l][z][d + small_k - 1 - i]

                    R[y][j] += val

            c += 1

        for y in range(po):
            for j in range(len(Q[0])):
                f += math.pow(R[y][j] - Q[y][j], 2)

        return f


    def find(self,small_k, P, poim, mu, base=['A','C','G','T']):
        xtmp = tools.startpwmrnd_base2(P, len(base))
        x0 = np.array([])
        for l in range(len(P)):
            x0 = np.append(x0, xtmp[l].flatten())

        L = len(poim[0])
        lb = np.ones(x0.shape) * 0.001
        ub = np.ones(x0.shape) * 0.99

        lenA = int(len(x0))
        lenk = int(len(x0)) / len(base)
        Aeq = np.zeros((lenk, lenA))
        beq = np.ones(lenk)
        c = 0
        c2 = 0
        for l in range(len(P)):
            for i in range(P[l]):
                for pk in range(c2 + i, c2 + int(P[l] * len(base)), P[l]):
                    Aeq[c + i, pk] = 1
            c += P[l]
            c2 += P[l] * len(base)
        cons = {'type': 'eq', 'fun': lambda x: np.dot(Aeq, x) - beq}
        bnds = []
        for i in range(len(x0)):
            bnds.append((lb[i], ub[i]))
        result = minimize(self.f_L2_2, x0, args=(poim, L, mu, small_k, P), method='SLSQP', jac=False,
                          bounds=bnds, constraints=cons)  # ,options={'ftol': 0.01})#

        self.motif_pwm = tools.converter(result.x, P)
        fopt = result.fun


        return self.motif_pwm


    def find2(self, POIMobj, motif_len, motif_start, base, path2pwm=None,solver="NLP"):
        self.motif_start = motif_start
        self.motif_len = motif_len
        x0 = tools.ini_pwm(motif_len, 1, len(base))[0]

        x0 = x0.flatten()

        lb = np.ones(x0.shape) * 0.001
        ub = np.ones(x0.shape) * 0.999
        iprint = 0
        maxIter = 1000
        ftol = 1e-04
        gradtol = 1e-03
        diffInt = 1e-05
        contol = 1e-02
        maxFunEvals = 1e04
        maxTime = 100

        lenA = int(len(x0))
        lenk = int(len(x0)) / len(base)
        Aeq = np.zeros((lenk, lenA))
        beq = np.ones(lenk)
        for i in range(lenk):
            for pk in range(i, lenA, lenk):
                Aeq[i, pk] = 1

                # ,Aeq=Aeq,beq=beq,
        cons = {'type': 'eq', 'fun': lambda x: np.dot(Aeq, x) - beq}
        bnds = []
        for i in range(len(x0)):
            bnds.append((lb[i], ub[i]))
        # bnds = np.vstack((lb,ub))


        if solver == "ralg":
            from openopt import NLP
            p = NLP(self.f_L2, x0,lb=lb, ub=ub, Aeq=Aeq,beq=beq, args=(POIMobj.gPOIM,POIMobj.L,motif_start,POIMobj.small_k,motif_len),  diffInt=diffInt, ftol=ftol, plot=0, iprint=iprint,maxIter = maxIter, maxFunEvals = maxFunEvals, show=False, contol=contol)
            result = p._solve(solver)
            x = result.xf
            f = result.ff
        elif solver == "LBFGSB":
            x, f, d = fmin_l_bfgs_b(self.f_L2, x0,
                                    args=(POIMobj.gPOIM, POIMobj.L, motif_start, POIMobj.small_k, motif_len),
                                    approx_grad=True)#constraints=cons)#
        elif solver == "SLSQP":
            result = minimize(self.f_L2, x0,args=(POIMobj.gPOIM, POIMobj.L, motif_start, POIMobj.small_k, motif_len),method='SLSQP',bounds=bnds,constraints=cons)
            x = result.x
            f = result.fun
        self.motif_pwm = np.reshape(x, (4, motif_len))
        fopt = f
        self.normalize()
        if not(path2pwm is None):
            np.savetxt(path2pwm, self.poim_norm)

        return self.motif_pwm

    def normalize(self):
        motif_pos = self.motif_pwm - np.min(self.motif_pwm)
        self.poim_norm = motif_pos / np.sum(motif_pos, axis=0)
        return self.poim_norm


    def proofmotif(self,x,y):
        self.prob = np.ones(len(x))
        for i,xi in enumerate(x):
            for p in range(self.motif_len):
                if (xi[self.motif_start+p] == 'A'):
                    self.prob[i] *=self.poim_norm[0][p]
                elif (xi[self.motif_start+p] == 'C'):
                    self.prob[i] *= self.poim_norm[1][p]
                elif (xi[self.motif_start+p] == 'G'):
                    self.prob[i] *= self.poim_norm[2][p]
                else:
                    self.prob[i] *= self.poim_norm[3][p]
        import scipy.stats

        group1 = self.prob[np.where(y==1)[0]]
        group2 = self.prob[np.where(y==-1)[0]]
        f,p_value = scipy.stats.f_oneway(group1, group2)

        from sklearn import metrics
        fpr, tpr, thresholds = metrics.roc_curve(y, self.prob)
        auc = metrics.auc(fpr, tpr)
        return p_value,auc


    def probmotifs(self,x,y,motif_start,motif_len,motifs):

        prob = np.ones(len(x))
        for m in range(len(motif_start)):
            for i, xi in enumerate(x):
                for p in range(motif_len[m]):
                    if (xi[motif_start[m] + p] == 'A'):
                        prob[i] *= motifs[m][0][p]
                    elif (xi[motif_start[m] + p] == 'C'):
                        prob[i] *= motifs[m][1][p]
                    elif (xi[motif_start[m] + p] == 'G'):
                        prob[i] *= motifs[m][2][p]
                    else:
                        prob[i] *= motifs[m][3][p]
        import scipy.stats

        group1 = prob[np.where(y == 1)[0]]
        group2 = prob[np.where(y == -1)[0]]
        f, p_value = scipy.stats.f_oneway(group1, group2)

        from sklearn import metrics
        fpr, tpr, thresholds = metrics.roc_curve(y, prob)
        auc = metrics.auc(fpr, tpr)
        return p_value, auc


from bisect import bisect_left

def takeClosest(myList, myNumber,L):
    """
    Assumes myList is sorted. Returns closest value to myNumber before myNumber and after myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return [0, myList[0]]
    if pos == len(myList):
        return  [myList[-1], L] # L is len of maxpoim
    before = myList[pos - 1]
    after = myList[pos]
    return [before,after]


def get_motif_pos_len(poim,N=1):
    motif_start=[]
    motif_len=[]
    maxpoim = np.sum(poim, axis=0)
    baseline = np.mean(maxpoim)
    minval = np.min(maxpoim)
    for i in range(N):

        peak = np.argmax(maxpoim)
        nonpeak = np.where(maxpoim < baseline)[0]
        #print peak
        #print nonpeak
        before,after =  takeClosest(nonpeak,peak,len(maxpoim))
        #print before
        #print " after " + str(after)
        if before == 0:
            motif_start.append(before)
            motif_len.append(after - before)
        else:
            motif_start.append(before+1)
            motif_len.append(after-before)

        maxpoim[motif_start[i]:motif_start[i]+motif_len[i]] = minval
    print motif_start
    print motif_len
    return [motif_start,motif_len]

def get_motif_pos_len_similar_diffPOIM(poim,N=1):
    motif_start=[]
    motif_len=[]
    maxpoim = np.max(poim, axis=0)
    baseline = np.mean(maxpoim) + np.std(maxpoim)
    minval = np.min(maxpoim)
    for i in range(N):

        peak = np.argmax(maxpoim)
        nonpeak = np.where(maxpoim < baseline)[0]
        #print peak
        #print nonpeak
        before,after =  takeClosest(nonpeak,peak,len(maxpoim))
        #print before
        #print " after " + str(after)
        if before == 0:
            motif_start.append(before)
            motif_len.append(after - before)
        else:
            motif_start.append(before+1)
            motif_len.append(after-before)

        maxpoim[motif_start[i]:motif_start[i]+motif_len[i]] = minval
    print motif_start
    print motif_len
    return [motif_start,motif_len]


def get_pval_auc(max_num_motif,gPOIMobj,x,y):
    pvals = []
    aucvals = []
    Mobj = Motif()
    for num_motif in range(1,max_num_motif+1):
        motif_start, motif_len = get_motif_pos_len(gPOIMobj.diffPOIM,num_motif)
        motifs = []
        for i in range(len(motif_start)):
            motifs.append(Mobj.find2(gPOIMobj, motif_len[i], motif_start[i], ['A', 'C', 'G', 'T'], path2pwm=None , solver="SLSQP")) # "SLSQP"
        p,auc= Mobj.probmotifs(x,y, motif_start, motif_len, motifs)
        pvals.append(p)
        aucvals.append(auc)
    return pvals,aucvals

def seqLogo(path2pwm):
    from rpy2.robjects.packages import importr
    import rpy2.robjects as robjects
    importr("seqLogo")
    importr("grid")
    m = robjects.r('''m = read.table( ( "%(path2pwm)s"))''' % locals())
    pwm = robjects.r.makePWM(m)
    s = path2pwm[:-4] + "_seqLogo.pdf"
    robjects.r.pdf(file=path2pwm[:-4] + "_seqLogo.pdf")
    robjects.r.seqLogo(pwm, robjects.r("ic=FALSE"), xfontsize=34, yfontsize=34)
    robjects.r("dev.off()")

