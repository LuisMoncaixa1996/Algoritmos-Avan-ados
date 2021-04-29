from EvolAlgorithm import EvolAlgorithm
from Popul import PopulInt, PopulReal
from MotifFinding import MotifFinding
from MyMotifs import MyMotifs



def createMatZeros(nl, nc):
    res = []
    for _ in range(0, nl):
        res.append([0]*nc)
    return res


def printMat(mat):
    for i in range(0, len(mat)):
        for j in range(len(mat[i])):
            print(f"{mat[i][j]:.3f}", end=' ')
        print()


class EAMotifsInt (EvolAlgorithm):
    def __init__(self, popsize, numits, noffspring, filename):
        self.motifs = MotifFinding()
        self.motifs.readFile(filename, "dna")
        indsize = len(self.motifs)
        EvolAlgorithm.__init__(self, popsize, numits, noffspring, indsize)

    def initPopul(self, indsize):
        maxvalue = self.motifs.seqSize(0) - self.motifs.motifSize
        self.popul = PopulInt(self.popsize, indsize,
                              maxvalue, [])

    def evaluate(self, indivs):
        for i in range(len(indivs)):
            ind = indivs[i]
            sol = ind.getGenes()
            fit = self.motifs.score(sol)
            fit1 = self.motifs.scoreMult(sol)
            ind.setFitness(fit)
            ind.setmultFitness(fit1)
        
        


class EAMotifsReal (EvolAlgorithm):
    def __init__(self, popsize, numits, noffspring, filename):
        self.motifs = MotifFinding()
        self.motifs.readFile(filename, "dna")
        indsize = self.motifs.motifSize * len(self.motifs.alphabet)
        EvolAlgorithm.__init__(self, popsize, numits, noffspring, indsize)
        

    def initPopul(self, indsize):
        maxvalue = self.motifs.seqSize(0) - self.motifs.motifSize
        self.popul = PopulReal(self.popsize, indsize, 0,
                              maxvalue, [])


    def vec_to_pwm(self,v):
        n_alph = len(self.motifs.alphabet)
        n_motif = self.motifs.motifSize
        pwm = createMatZeros(n_alph, n_motif)
        for i in range(0,len(v),n_alph):
            col_idx = int(i / n_alph)
            col = v[i:i+n_alph]
            soma = sum(col)
            for j in range(n_alph):
                pwm[j][col_idx] = col[j] / soma
        return pwm
    
    def probabSeq(self, seq):
        res = 1.0
        for i in range(self.motifs.motifSize):
            lin = self.motifs.alphabet.index(seq[i])
            res *= self.motifs.pwm[lin][i]
        return res
    
    def mostProbableSeq(self, seq):
        maximo = -1.0
        maxind = -1
        for k in range(len(seq)-self.motifs.motifSize):
            p = self.probabSeq(seq[k:k + self.motifs.motifSize])
            if(p > maximo):
                maximo = p
                maxind = k
        return maxind
    
    
    def evaluate(self, indivs):
        
        for i in range(len(indivs)):
            ind = indivs[i]
            sol = ind.getGenes()
            pwm = self.vec_to_pwm(sol)
            mtf = MyMotifs(pwm=pwm, alphabet=self.motifs_finding.alphabet)
            s = []
            for j in range(len(self.motifs_finding.seqs)):
                seq = self.motifs_finding.seqs[j].seq
                p = mtf.mostProbableSeq(seq)
                s.append(p)  
            fit1 = self.motifs.scoreMult(s,self.motifs.pwm) # com scoreMult
            fit = self.motifs.score(s) 
            ind.setFitness(fit)
            ind.setmultFitness(fit1)


def test1():
    ea = EAMotifsInt(100, 4, 50, "exemploMotifs.txt")
    ea.run()
    ea.printBestSolution()


def test2():
    ea = EAMotifsReal(100, 10, 50, "exemploMotifs.txt")
    ea.run()
    ea.printBestSolution()


#test1()
test2()
