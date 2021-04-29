# -*- coding: utf-8 -*-

class SuffixTree:
    
    def __init__(self):
        self.nodes = { 0:(-1,{}) } # root node
        self.num = 0
        self.seq = ''
    
    def print_tree(self):
        for k in self.nodes.keys():
            if self.nodes[k][0] < 0:
                print (k, "->", self.nodes[k][1]) 
            else:
                print (k, ":", self.nodes[k][0])
                
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
    
    def suffix_tree_from_seq(self, text):
        self.seq = text
        t = text + "$"
        for i in range(len(t)):
            self.add_suffix(t[i:], i)
            
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
        res = []
        if self.nodes[node][0] >=0: 
            res.append(self.nodes[node][0])            
        else:
            for k in self.nodes[node][1].keys():
                newnode = self.nodes[node][1][k]
                leafes = self.get_leafes_below(newnode)
                res.extend(leafes)
        return res
    
    def nodes_below(self, node): #exercicio 1 a)
        lb = list(self.nodes[node][1].values()) # Colocar em uma lista os nº do nó de destino em cada arco correspondente ao nó selecionado
        for i in lb:
            lb.extend(list(self.nodes[i][1].values())) # Aumentar a lista com os nº do nó de destino dos arcos já existentes na lista
        return lb
    
    def matches_prefix(self, prefix): # exercicio 1 b)
        pre = SuffixTree.find_pattern(self, prefix) #Definir as posições iniciais onde se encontra o prefixo
        if pre == None or []:
            return None
        else:
            match = []
            for i in pre:
                match.append(self.seq[i:]) #Adicionar à lista os padrões cujo, as posições iniciais são as correspondentes ao prefixo
            for j in range(len(match)):
                dp = len(match[j]) #Tamanho do padrão encontrado
                c = 1
                while dp > len(prefix): 
                    match.append(match[j][:-c]) #Adicionar à lista os padrões que contêm o prefixo, retirando um caracter de cada vez até o tamanho do padrão ser 0
                    dp -= 1
                    c += 1
            distintpattern = []
            for l in match:      #Colocar numa nova lista os padrões encontrados, retirando os duplicados
                if l not in distintpattern:
                    distintpattern.append(l)
            return distintpattern
                        

def test():
    seq = "TACTA"
    st = SuffixTree()
    st.suffix_tree_from_seq(seq)
    st.print_tree()
    #print (st.find_pattern("TA"))
    #print (st.find_pattern("ACG"))
    print(st.nodes_below(0))
    #print(st.get_leafes_below(7))

def test2():
    seq = "TACTAGHF"
    st = SuffixTree()
    st.suffix_tree_from_seq(seq)
    print (st.find_pattern("TA"))
    print(st.matches_prefix('TA'))


#test()
#print()
test2()

      
            
    
    
