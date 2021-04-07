# -*- coding: utf-8 -*-
"""

"""

from MySeq import MySeq
from MyMotifs import MyMotifs

class MotifFinding:
    
    def __init__(self, size = 8, seqs = None):
        self.motifSize = size
        if (seqs != None):
            self.seqs = seqs
            self.alphabet = seqs[0].alfabeto()
        else:
            self.seqs = []
                    
    def __len__ (self):
        return len(self.seqs)
    
    def __getitem__(self, n):
        return self.seqs[n]
    
    def seqSize (self, i):
        return len(self.seqs[i])
    
    def readFile(self, fic, t):
        for s in open(fic, "r"):
            self.seqs.append(MySeq(s.strip().upper(),t))
        self.alphabet = self.seqs[0].alfabeto()
        
        
    def createMotifFromIndexes(self, indexes):
        pseqs = []
        for i,ind in enumerate(indexes):
            pseqs.append( MySeq(self.seqs[i][ind:(ind+self.motifSize)], self.seqs[i].tipo) )
        return MyMotifs(pseqs)
        
        
    # SCORES
        
    def score(self, s):
        score = 0
        motif = self.createMotifFromIndexes(s)
        motif.doCounts()
        mat = motif.counts
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for  i in range(1, len(mat)):
                if mat[i][j] > maxcol: 
                    maxcol = mat[i][j]
            score += maxcol
        return score
   
    def scoreMult(self, s):
        score = 1.0
        motif = self.createMotifFromIndexes(s)
        motif.createPWM()
        mat = motif.pwm
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for  i in range(1, len(mat)):
                if mat[i][j] > maxcol: 
                    maxcol = mat[i][j]
            score *= maxcol
        return score     
       
    # EXHAUSTIVE SEARCH
       
    def nextSol (self, s):
        nextS = [0]*len(s)
        pos = len(s) - 1     
        while pos >=0 and s[pos] == self.seqSize(pos) - self.motifSize:
            pos -= 1
        if (pos < 0): 
            nextS = None
        else:
            for i in range(pos): 
                nextS[i] = s[i]
            nextS[pos] = s[pos]+1;
            for i in range(pos+1, len(s)):
                nextS[i] = 0
        return nextS
        
    def exhaustiveSearch(self):
        melhorScore = -1
        res = []
        s = [0]* len(self.seqs)
        while (s!= None):
            sc = self.score(s)
            if (sc > melhorScore):
                melhorScore = sc
                res = s
            s = self.nextSol(s)
        return res
     
    # BRANCH AND BOUND     
     
    def nextVertex (self, s):
        res =  []
        if len(s) < len(self.seqs): # internal node -> down one level
            for i in range(len(s)): 
                res.append(s[i])
            res.append(0)
        else: # bypass
            pos = len(s)-1 
            while pos >=0 and s[pos] == self.seqSize(pos) - self.motifSize:
                pos -= 1
            if pos < 0: res = None # last solution
            else:
                for i in range(pos): res.append(s[i])
                res.append(s[pos]+1)
        return res
    
    
    def bypass (self, s):
        res =  []
        pos = len(s) -1
        while pos >=0 and s[pos] == self.seqSize(pos) - self.motifSize:
            pos -= 1
        if pos < 0: res = None 
        else:
            for i in range(pos): res.append(s[i])
            res.append(s[pos]+1)
        return res
        
    def branchAndBound (self):
        melhorScore = -1
        melhorMotif = None
        size = len(self.seqs)
        s = [0]*size
        while s != None:
            if len(s) < size:
                optimScore = self.score(s) + (size-len(s)) * self.motifSize
                if optimScore < melhorScore: s = self.bypass(s)
                else: s = self.nextVertex(s)
            else:
                sc = self.score(s)
                if sc > melhorScore:
                    melhorScore = sc
                    melhorMotif  = s
                s = self.nextVertex(s)
        return melhorMotif

    # Consensus (heuristic)
  
    def heuristicConsensus(self):
        # procuras as posições para o motif nas 2 primeiras sequências
        # procura exaustiva nas duas primeiras sequências
        mf = MotifFinding(self.motifSize, self.seqs[:2])
        s = mf.exhaustiveSearch() #(1,3)
        #avalia a melhor posição para cada uma das sequencias
        #seguintes uma a uma, guardando a melhor posição (maximiza o score)
        for i in range(2,len(self.seqs)):
            s.append(0)
            score = -1
            pos = 0
            for j in range(self.seqSize(i)-self.motifSize+1):
                s[i] = j
                score_atual = self.score(s)
                if score_atual > score:
                    score = score_atual
                    pos = j
                s[i] = pos
        return s

    # Consensus (heuristic)

    def heuristicStochastic (self):
        from random import randint
        s = [0] * len(self.seqs)
        #passo1 
        #Iniciar todas as posições com valores aleatórios
        for i in range(len(self.seqs)):
            #randint(A,B) =>   A<=x<=B
            s[i] = randint(0, self.seqSize(i)-self.motifSize)
        #passo2
        best_score = self.score(s)
        improve = True
        while improve:
            #constroi o perfil com base nas posições iniciais s
            motif = self.createMotifFromIndexes(s)
            motif.createPWM()
            #avalia a melhor posição inicial para cada sequência com base no perfil
            for i in range(len(self.seqs)):
                s[i] = motif.mostProbableSeq(self.seqs[i])
            #verifica se houve melhoria
            scr = self.score(s)
            if scr > best_score:
                best_score = scr
            else:
                improve = False
        return s

    # Gibbs sampling 

    def gibbs (self, iterations):
        from random import randint
        s = []
        for i in range(len(self.seqs)):
            s.append(randint(0, len(self.seqs[0]) - self.motifSize - 1))
        best_s = list(s)
        best_score = self.scoreMult(s)
        for it in range(iterations):
            # randomly pick a sequence
            seq_idx = randint(0, len(self.seqs) - 1)
            seq_sel = self.seqs[seq_idx]
            s.pop(seq_idx)
            removed = self.seqs.pop(seq_idx)
            motif = self.createMotifFromIndexes(s)
            motif.createPWM()
            self.seqs.insert(seq_idx, removed)
            r = motif.probAllPositions(seq_sel)
            pos = self.roulette(r)
            s.insert(seq_idx, pos)
            score = self.scoreMult(s)
            if score > best_score:
                best_score = score
                best_s = list(s)
        return best_s

    def roulette(self, f):
        from random import random
        tot = 0.0
        for x in f: tot += (0.01+x)
        val = random()* tot
        acum = 0.0
        ind = 0
        while acum < val:
            acum += (f[ind] + 0.01)
            ind += 1
        return ind-1



# Exercicios aula 5

# Alteração de algumas funções de forma a ser possível analisar através do algoritmo heurístico estocástico para pseudocontagens

    def pseudscore(self,s):
        score = 0
        motif = self.createMotifFromIndexes(s)
        motif.doCounts()
        mat = []
        for l in range(len(motif.counts)):
            row = []
            for k in range(len(motif.counts[l])):
                row.append(motif.counts[l][k] + 1)
            mat.append(row)
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for  i in range(1, len(mat)):
                if mat[i][j] > maxcol: 
                    maxcol = mat[i][j]
            score += maxcol
        return score
    
    def pseudprobabSeq (self, seq, pwm):
        res = 1.0
        for i in range(self.motifSize):
            lin = self.alphabet.index(seq[i])
            res *= pwm[lin][i]
        return res

    def pseudmostProbableSeq(self, seq, pwm):
        maximo = -1.0
        maxind = -1
        for k in range(len(seq)-self.motifSize):
            p = self.pseudprobabSeq(seq[k:k+ self.motifSize], pwm)
            if(p > maximo):
                maximo = p
                maxind = k
        return maxind
    
    def pseudpwm(self, pwm):
        pseudopwm = []
        for i in range(len(pwm)):
            row = []
            for j in range(len(pwm[0])):
                row.append(pwm[i][j] + 0.1)
            pseudopwm.append(row)
        return pseudopwm
    
    def pseudprobAllPositions(self, seq, pwm):
        res = []
        for k in range(len(seq)-self.motifSizes+1):
            res.append(self.pseudprobabSeq(seq, pwm))
        return res
    
    def pseudscoreMult(self, s):
        score = 1.0
        motif = self.createMotifFromIndexes(s)
        motif.createPWM()
        mat = motif.pwm
        pseudopwm = self.pseudpwm(mat)
        for j in range(len(pseudopwm[0])):
            maxcol = pseudopwm[0][j]
            for  i in range(1, len(pseudopwm)):
                if pseudopwm[i][j] > maxcol: 
                    maxcol = pseudopwm[i][j]
            score *= maxcol
        return score
    
# Analise através algoritmo heurístico estocástico para pseudocontagens:
    def pseudheuristicStochastic (self):
        from random import randint
        s = [0] * len(self.seqs)
        #passo1 
        #Iniciar todas as posições com valores aleatórios
        for i in range(len(self.seqs)):
            #randint(A,B) =>   A<=x<=B
            s[i] = randint(0, self.seqSize(i)-self.motifSize)
        #passo2
        best_score = self.pseudscore(s)
        improve = True
        while improve:
            #constroi o perfil com base nas posições iniciais s
            motif = self.createMotifFromIndexes(s)
            motif.createPWM()
            pwm = motif.pwm
            pseudopwm = self.pseudpwm(pwm)
            #avalia a melhor posição inicial para cada sequência com base no perfil
            for i in range(len(self.seqs)):
                s[i] = self.pseudmostProbableSeq(self.seqs[i], pseudopwm)
            #verifica se houve melhoria
            scr = self.pseudscore(s)
            if scr > best_score:
                best_score = scr
            else:
                improve = False
        return s
    
    def pseudgibbs (self, iterations):
        from random import randint
        s = []
        for i in range(len(self.seqs)):
            s.append(randint(0, len(self.seqs[0]) - self.motifSize - 1))
        best_s = list(s)
        best_score = self.pseudscore(s)
        for it in range(iterations):
            # randomly pick a sequence
            seq_idx = randint(0, len(self.seqs) - 1)
            seq_sel = self.seqs[seq_idx]
            s.pop(seq_idx)
            removed = self.seqs.pop(seq_idx)
            motif = self.createMotifFromIndexes(s)
            pwm = motif.createPWM()
            pseudopwm = self.pseudpwm(pwm)
            self.seqs.insert(seq_idx, removed)
            r = motif.pseudprobAllPositions(seq_sel,pseudopwm)
            pos = self.roulette(r)
            s.insert(seq_idx, pos)
            score = self.pseudscore(s)
            if score > best_score:
                best_score = score
                best_s = list(s)
        return best_s
        
# tests

def test1():  
    sm = MotifFinding()
    sm.readFile("exemploMotifs.txt","dna")
    sol = [25,20,2,55,59]
    sa = sm.score(sol)
    print(sa)
    scm = sm.scoreMult(sol)
    print(scm)

def test2():
    print ("Test exhaustive:")
    seq1 = MySeq("ATAGAGCTGA","dna")
    seq2 = MySeq("ACGTAGATGA","dna")
    seq3 = MySeq("AAGATAGGGG","dna")
    mf = MotifFinding(3, [seq1,seq2,seq3])
    sol = mf.exhaustiveSearch()
    print ("Solution", sol)
    print ("Score: ", mf.score(sol))
    print("Consensus:", mf.createMotifFromIndexes(sol).consensus())

    print ("Branch and Bound:")
    sol2 = mf.branchAndBound()
    print ("Solution: " , sol2)
    print ("Score:" , mf.score(sol2))
    print("Consensus:", mf.createMotifFromIndexes(sol2).consensus())
    
    print ("Heuristic consensus: ")
    sol1 = mf.heuristicConsensus()
    print ("Solution: " , sol1)
    print ("Score:" , mf.score(sol1))

def test3():
    mf = MotifFinding()
    mf.readFile("exemploMotifs.txt","dna")
    print ("Branch and Bound:")
    sol = mf.branchAndBound()
    print ("Solution: " , sol)
    print ("Score:" , mf.score(sol))
    print("Consensus:", mf.createMotifFromIndexes(sol).consensus())

def test4():
    mf = MotifFinding()
    mf.readFile("exemploMotifs.txt","dna")
    print("Heuristic stochastic")
    sol = mf.heuristicStochastic()
    print ("Solution: " , sol)
    print ("Score:" , mf.score(sol))
    print ("Score mult:" , mf.scoreMult(sol))
    print("Consensus:", mf.createMotifFromIndexes(sol).consensus())
    
    sol2 = mf.gibbs(1000)
    print ("Score:" , mf.score(sol2))
    print ("Score mult:" , mf.scoreMult(sol2))

def test5():
    mf = MotifFinding()
    mf.readFile("exemploMotifs.txt","dna")
    print("Heuristic stochastic")
    sol = mf.pseudheuristicStochastic()
    print ("Solution: " , sol)
    print ("Score:" , mf.pseudscore(sol))
    print ("Score mult:" , mf.scoreMult(sol))
    print("Consensus:", mf.createMotifFromIndexes(sol).consensus())
    
    sol2 = mf.gibbs(1000)
    print ("Score:" , mf.pseudscore(sol2))
    print ("Score mult:" , mf.scoreMult(sol2))
    
test1()
print('-=' * 30)
test2()
print('-=' * 30)
test3()
print('-=' * 30)
test4()
print('-=' * 30)
test5()
