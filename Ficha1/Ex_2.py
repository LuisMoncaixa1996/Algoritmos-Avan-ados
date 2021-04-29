class SuffixTree:
    
    def __init__(self):
        self.nodes = { 0:(-1,{}) } # root node
        self.num = 0
        self.seq1 = ''
        self.seq2 = ''
    
    
    def print_tree(self):
        for k in self.nodes.keys():
            if self.nodes[k][0] != -1:
                iseq,c = self.nodes[k][0]
            else:
                iseq = self.nodes[k][0]
                c = ''
            if iseq < 0:
                print (k, "->", self.nodes[k][1]) 
            else:
                print (k, ":", iseq,c)
                
    def add_node(self, origin, symbol, leafnum = -1):
        self.num += 1
        self.nodes[origin][1][symbol] = self.num
        self.nodes[self.num] = (leafnum,{})
        
    def add_suffix(self, p, sufnum):
        pos = 0
        node = 0
        while pos < len(p):
            if p[pos] not in self.nodes[node][1].keys():
                if pos == len(p) - 1:
                    self.add_node(node, p[pos], sufnum)
                else:
                    self.add_node(node, p[pos])
            node = self.nodes[node][1][p[pos]]
            pos += 1
    
    def suffix_tree_from_seq(self, seq1, seq2): 
        seq1 += "$"
        seq2 += "#"
        self.seq1 = seq1 #Adiciona as sequências ao init 
        self.seq2 = seq2
        for i in range(len(seq1)):
            self.add_suffix(seq1[i:], (0,i)) #Criar túpulos correspondentes às leafs onde o 0 e 1 corresponde à sequência e o i corresponde ao arco 
        for i in range(len(seq2)):
            self.add_suffix(seq2[i:], (1,i))
            
    def find_pattern(self, pattern):
        pos = 0
        node = 0
        for pos in range(len(pattern)):
            if pattern[pos] in self.nodes[node][1].keys():
                node = self.nodes[node][1][pattern[pos]]
                pos += 1
            else:
                return None
        return self.get_leafes_below(node)
        

    def get_leafes_below(self, node):
        seq1 = []
        seq2 = []
        if self.nodes[node][0] != -1:
            iseq,c = self.nodes[node][0]
        else:
            iseq = self.nodes[node][0]
            c = ''
        if iseq >=0:
            if iseq == 0:
                seq1.append(c)
            else:
                seq2.append(c)
        else:
            for k in self.nodes[node][1].keys():
                newnode = self.nodes[node][1][k]
                s0,s1 = self.get_leafes_below(newnode)
                seq1.extend(s0)
                seq2.extend(s1)
        return (seq1, seq2)
    
    def largestCommonSubstring(self):
        padrao1 = ''
        count1 = 0
        for i in range(len(self.seq1)): # Vai iterar sobre o tamanho da primeira sequência
            padrao = ''
            contagem = 0
            for j in range(len(self.seq2)): # Iteração sobre o tamanho da segunda sequência
                if self.seq1[i] == self.seq2[j]: # Caso ocorra match
                    padrao += self.seq1[i] #Adiciona ao padrão a letra correspondente e aumenta a contagem e avança uma posição em relação à primeira sequência
                    contagem +=1
                    i +=1
                else:
                    if contagem > count1: #Se não ocorrer match vai verificar se a contagem de posições é superior á contagem inicial, caso se observe, adiciona o padrão obtido e determina define a contagem atual
                        padrao1 = padrao
                        count1 = contagem
                    contagem = 0 #Se a contagem não for superior dá-se reset à contagem e ao padrão
                    padrao = ''
        print(padrao1)
                        

def test():
    seq1 = "TACTA"
    seq2 = "TCATCAAT"
    st = SuffixTree()
    st.suffix_tree_from_seq(seq1,seq2)
    st.print_tree()
    print (st.find_pattern("TA"))
    print (st.find_pattern("ACG"))
    print(st.get_leafes_below(2))

def test2():
    seq1 = "GADGGFGGGGGGGGGHLDHOFOIGJKCTA"
    seq2 = "TAAGADGGFGGGGGGGGGHLDHOFOIGJKCTA"
    st = SuffixTree()
    st.suffix_tree_from_seq(seq1,seq2)
    print (st.find_pattern("TA"))
    st.largestCommonSubstring()
    


test()
#print()
#test2()