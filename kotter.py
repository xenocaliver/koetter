#!/usr/local/sage/default/sage -python
# coding: UTF-8

r"""
This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

decode a code word with BMS algorithm

AUTHORS:

- Akiyoshi Hashimoto (2014-06): initial implementation

"""

# import command line arguments functions
import sys
# for option parser
import optparse
# import all sage function
from sage.all import *
#
F = GF(8, 'a')
a = F.gen()
K = F['z']
z = K.gen()
gamma = 3
BASISINDEX = FiniteEnumeratedSet([3, 5, 7])
GAP = FiniteEnumeratedSet([1, 2, 4])
l = []
asigma = []
alambda = []
delta = []
Delta = []
xsigma = []
xlambda = []
j = []
MATRIXSIZE = 128

def initSyndrome():
    global F
    global a
    global MATRIXSIZE

    S = matrix(F, MATRIXSIZE, MATRIXSIZE)      # define syndrome matrix
    S[0, 0] = a*a*a*a
    S[3, 0] = a
    S[5, 0] = a*a*a
    S[6, 0] = a**6
    S[7, 0] = 1
    S[8, 0] = a**6
    S[10, 0] = a**5
    S[11, 0] = a*a
    S[12, 0] = a**4
    S[13, 0] = a**5
    S[0, 3] = a
    S[3, 3] = a**6
    S[5, 3] = a**6
    S[7, 3] = a**5
    S[8, 3] = a*a
    S[9, 3] = a*a*a*a
    S[10, 3] = a**5
    S[0, 5] = a*a*a
    S[3, 5] = a**6
    S[5, 5] = a**6
    S[6, 5] = a*a
    S[7, 5] = a**4
    S[8, 5] = a
    S[0, 6] = a**6
    S[5, 6] = a*a
    S[6, 6] = a*a*a*a
    S[7, 6] = a**5
    S[0, 7] = 1
    S[3, 7] = a**5
    S[5, 7] = a*a*a*a
    S[6, 7] = a**5
    S[0, 8] = a**6
    S[3, 8] = a*a
    S[5, 8] = a
    S[3, 9] = a*a*a*a
    S[0, 10] = a**5
    S[3, 10] = a**5
    S[0, 11] = a*a
    S[0, 12] = a*a*a*a
    S[0, 13] = a**5
    return S

def inductionSyndrome(i, j):
    global F
    global K
    global a
    global z
    global gamma
    global BASISINDEX
    global GAP
    global S

    if i in GAP or j in GAP:
        return 0
    elif (i not in GAP) and (j not in GAP) and (i%gamma != 0) and (j%gamma != 0):
        return S[i + j - 7, 0]
    else:
        return 0

def initialize():
    global F
    global K
    global a
    global z
    global gamma
    global BASISINDEX
    global GAP
    global l
    global asigma
    global alambda
    global xsigma
    global xlambda
    global j
    global Delta
    global delta

    # initialize Syndrome
    initSyndrome()
    # initialize l[i]
    basislist = BASISINDEX.list()
    l = [0]*len(basislist)
    for i in range(0, len(l)):
        x = 0
        while True:
            ingap = (x in GAP)
            k = x%gamma
            if ingap == False:
                if k == i:
                    l[i] = x
                    break
            x += 1
    # initialize a_sigma
    asigma = [0]*len(l)
    alambda = [0]*len(l)
    xsigma  = [0]*len(l)
    xlambda = [0]*len(l)
    Delta = [0]*len(l)
    delta = [0]*len(l)
    j = [0]*len(l)
    for i in range(0, len(l)):
        xsigma[i] = 1 + 0*z
        asigma[i] = -l[i]
        xlambda[i] = 0 + 0*z
        alambda[i] = l[i] - gamma
    print "{0} {1} {2} {3}".format(asigma, alambda, xsigma, xlambda)

def update(r, vote):
    global F
    global K
    global a
    global z
    global gamma
    global BASISINDEX
    global GAP
    global l
    global asigma
    global alambda
    global xsigma
    global xlambda
    global j
    global Delta
    global delta
    global S

    for i in range(0, len(asigma)):
        j[i] = asigma[i]%gamma

    if vote == False:
        for i in range(0, len(asigma)):
            csigma = xsigma[i].coeffs()
            if asigma[i] in GAP:
                Delta[i] = 0
            else:
                Delta[i] = 0
                for h in range(0, r + 1 - asigma[i] + 1):
                    for k in range(0, r + 1 - asigma[i] + 1):
                        csigma.append(0)
                    Delta[i] += csigma[h]*S[asigma[i], r + 1 - asigma[i] - h]

    for i in range(0, len(asigma)):
#        print "Delta = {0}".format(Delta[i])
        if Delta[i] != 0 and alambda[j[i]] < asigma[i]:
            delta[i] = 1
        else:
            delta[i] = 0
    sigmatmp = [0]*len(asigma)
    lambdatmp = [0]*len(asigma)
    asigmatmp = [0]*len(asigma)
    alambdatmp = [0]*len(alambda)
    for i in range(0, len(asigma)):
        sigmatmp[i] = xsigma[i]
        lambdatmp[i] = xlambda[i]
        asigmatmp[i] = asigma[i]
        alambdatmp[i] = alambda[i]
    for i in range(0, len(asigma)):
        xsigma[i] = sigmatmp[i] - Delta[i]*lambdatmp[j[i]]
        if Delta[i] == 0:
            xlambda[j[i]] = (1 - delta[i])*z*lambdatmp[j[i]]
        else:
            xlambda[j[i]] = delta[i]*z*sigmatmp[i]/Delta[i] + (1 - delta[i])*z*lambdatmp[j[i]]
    for i in range(0, len(asigma)):
        asigma[i] = (1 - delta[i])*asigmatmp[i] + delta[i]*alambdatmp[j[i]] + 1
#        print "delta = {0} asigmatmp = {1} alambdatmp = {2}".format(delta[i], asigmatmp[i], alambdatmp[j[i]])
        alambda[j[i]] = (1 - delta[i])*alambdatmp[j[i]] + delta[i]*asigmatmp[i]

    print "{0} {1} {2} {3} {4}".format(asigma, alambda, xsigma, xlambda, Delta, delta)
def estimatehatS(m):
    global S
    global F
    global MATRIXSIZE

    hatS = matrix(F, MATRIXSIZE, MATRIXSIZE)
    nrows = S.nrows()
    ncols = S.ncols()
    for i in range(0, nrows):
        for j in range(0, ncols):
            hatS[i, j] = S[i, j]
    hatS[m, 0] = 0
    for j in range(1, m + 1):
        hatS[m - j, j] = inductionSyndrome(m - j, j)
#    for j in range(0, m + 1):
#        print "hatS[{0}, {1}] = {2}".format(m - j, j, hatS[m - j, j])
    return hatS

def majorityVoting(m):
    global F
    global K
    global a
    global z
    global gamma
    global BASISINDEX
    global GAP
    global l
    global asigma
    global alambda
    global xsigma
    global xlambda
    global j
    global Delta
    global delta
    global S
    global hatS

    hatDelta = [0]*len(asigma)
    # estimate temporal estimate syndromes
    hatS = estimatehatS(m)
    # estimate temporal estimate discrepancy
    for i in range(0, len(asigma)):
        csigma = xsigma[i].coeffs()
        if asigma[i] in GAP:
            hatDelta[i] = 0
        else:
            Delta[i] = 0
            for h in range(0, m - asigma[i] + 1):
                for k in range(0, m - asigma[i] + 1):
                    csigma.append(0)
                hatDelta[i] += csigma[h]*S[asigma[i], m  - asigma[i] - h]
            hatDelta[i] += hatS[asigma[i], m - asigma[i]]
#        print "hatDelta[{0}] = {1}".format(i, hatDelta[i])
    # count ballot
    asigmatmp = [0]*len(asigma)
    alambdatmp = [0]*len(alambda)
    ballot = [0]*len(asigma)
    for i in range(0, len(asigma)):
        j[i] = asigma[i]%gamma
    for i in range(0, len(asigma)):
        asigmatmp[i] = asigma[i]
        alambdatmp[i] = alambda[i]
    for i in range(0, len(asigma)):
        asigma[i] = (1 - delta[i])*asigmatmp[i] + delta[i]*alambdatmp[j[i]] + 1
        alambda[j[i]] = (1 - delta[i])*alambdatmp[j[i]] + delta[i]*asigmatmp[i]
    for i in range(0, len(asigma)):
        ballot[i] += max(0, asigma[i] - alambda[j[i]])
#        print "ballot = {0}".format(ballot[i])
    for i in range(0, len(asigma)):
        asigma[i] = asigmatmp[i]
        alambda[j[i]] = alambdatmp[j[i]]
    # voting
    maxvote = 0
    maxindex = 0
    for i in range(0, len(asigma)):
        if maxvote < ballot[i]:
            maxvote = ballot[i]
            maxindex = i
    S[m, 0] = -hatDelta[maxindex]
    for i in range(0, len(asigma)):
        Delta[i] = hatDelta[i] + S[m, 0]
        print "Delta[{0}] = {1}".format(i, Delta[i])
    print "S[{0}, 0] = {1}".format(m, S[m, 0])
    return S[m, 0]

def updateSyndrome(ss, m):
    global S
    global gamma
    global GAP
    S[m, 0] = ss
    for j in range(1, m + 1):
        if ((m - j) in GAP) or (j in GAP):
            S[m - j, j] = 0
        elif ((m - j) not in GAP) and (j not in GAP) and ((m - j)%gamma != 0) and (j%gamma != 0):
            S[m - j, j] = ss + S[m - 7, 0]
        else:
            S[m - j, j] = ss
#        print "S[{0}, {1}] = {2}".format(m - j, j, S[m - j, j])
if __name__ == '__main__':
    S = initSyndrome()
    initialize()
    for r in range(-1, 13):
        update(r, False)
    ss = majorityVoting(14)
    updateSyndrome(ss, 14)
    update(14, True)
    ss = majorityVoting(15)
    updateSyndrome(ss, 15)
    update(15, True)
