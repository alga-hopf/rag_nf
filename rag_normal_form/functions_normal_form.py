

# # Normal form for right angled groups
# 
# These functions allow to compute explicitly the normal form for elements in a right-angled group (RAG)


import networkx as nx
from matplotlib import pyplot as plt
from itertools import combinations
from collections import Counter
import numpy as np

# Function for setting an element for which we want to compute the normal form. l is a list containing the representation of the element that
# we want to study. More precisely, each element of the list is the index of the generator appearing in that position

def element(l,GG):
    return [l+['nan'],GG]

# Function that defines the bicolored graph. It takes as input some data of the graph and the list of indices of black vertices.
# The data can be a dictionary of connections, the adjacency matrix, the incidence matrix or whatever object you can feed the NetworkX function
# "Graph" with
def bic_graph(vert_dict,Vb):
    G=nx.Graph(vert_dict)
    Vw=[v for v in range(1,G.number_of_nodes()+1) if v not in Vb]
    return [G,[Vb,Vw]]

# Function that draws the bicolored graph from which we want to build the group
def draw_graph(GG):
    G=GG[0]
    Vb=GG[1][0]
    Vw=GG[1][1]
    vertex_colors=['grey' if i+1 in Vb else 'white' for i in range(G.number_of_nodes())]
    nx.draw(G,node_color=vertex_colors,with_labels=True,linewidths=1,edgecolors='black')
    plt.show()

# Some functions used to handle graphs
    
# List containing the cliques of the graph G
def cliques(GG):
    G=GG[0]
    cliques=nx.enumerate_all_cliques(G)
    return cliques

# List containing the "extended" cliques of the graph, i.e., the cliques of the graph having as vertices also the
# inverses of the vertices in the original graph and where the edges are defined accordingly
def cliquesExt(GG):
    cliquesExt=[]
    for s in cliques(GG):
        S=s+[-s[i] for i in range(len(s))]
        for j in range(len(s),len(S)+1):
            tempComb=list(combinations(S,j))
            for k in range(len(tempComb)):
                cliquesExt.append(tempComb[k])
    cliquesExt=map(set,cliquesExt)
    return list(cliquesExt)


# Some functions used to compute the normal form for RAGs

# First decomposition of the element into conductors, without moving anything
# X=[Xword,G] and G=X[1] is the graph
def decomposition(X):
    cliques_ext=cliquesExt(X[1])
    dec=[]
    start=0
    for i in range(1,len(X[0])):
        if set(X[0][start:i]) in cliques_ext and set(X[0][start:i]).union(set([X[0][i]])) not in cliques_ext:
            dec.append(X[0][start:i])
            start=i
        else:
            continue
    return dec

# Left conductors of the elements X
def left_conductor(X):
    dec=decomposition(X)
    for i in range(len(dec)-1):
        cX=dec[len(dec)-1-i]
        decitemp=[cX[k] for k in range(len(cX))] 
        delete=[]
        for j in range(len(cX)):
            if set([cX[j]]).union(set(dec[len(dec)-2-i])) in cliquesExt(X[1]):
                dec[len(dec)-2-i].append(cX[j])
                delete.append(j)
        dec[len(dec)-1-i]=[decitemp[k] for k in range(len(decitemp)) if k not in delete]
    return dec

# Normal form for an element in a RAG, expressed as a list of dictionaries (the conductors are the indices of the
# list). Each dictionary has as keys the vertices of the graph and as values the exponents of the 
# corresponding generators in the corresponding conductor.
def normal_form_dict(X):
    XX=left_conductor(X)
    conductors=[]
    for i in range(len(XX)):
        expCond={}
        exponents=Counter(XX[i])
        for j in range(1,X[1][0].number_of_nodes()+1):
            if j in X[1][1][0]:
                expCond[j]=(exponents[j]-exponents[-j])%2
            else:
                expCond[j]=exponents[j]-exponents[-j]
        conductors.append(expCond)
    return conductors

# Normal form for an element in a RAG, expressed as a list of lists (the conductors are the indices of the
# external list). Each list has length equal to the numbert of vertices and in each entry there are the exponents
# of the corresponding generator in the corresponding conductor
def normal_form(X):
    XX=left_conductor(X)
    conductors=[]
    for i in range(len(XX)):
        expCond=[]
        exponents=Counter(XX[i])
        for j in range(1,X[1][0].number_of_nodes()+1):
            if j in X[1][1][0]:
                expCond.append((exponents[j]-exponents[-j])%2)
            else:
                expCond.append(exponents[j]-exponents[-j])
        conductors.append(expCond)
    conductors=np.array(conductors)
    p=[]
    for i in range(len(conductors)):
        if (conductors[i]==np.zeros(X[1][0].number_of_nodes())).all():
            p.append(i)
    conductors=np.delete(conductors,p,axis=0)
            
    flag=1
    while flag==1:
        positions=[]
        flag=0
        for i in range(len(conductors)-1):
            if (conductors[i]==conductors[i+1]).all():
                flag=1
                positions.append(i)
                positions.append(i+1)
        if positions!=[]:
            conductors=np.delete(conductors,positions,axis=0)
    
    return conductors


# Normal form of an element X returned as a Tietze list
def normal_form_tietze(X):
    x_nf=normal_form(X)
    tietze_list=[]
    for i in range(len(x_nf)):
        for j in range(len(x_nf[i])):
            if x_nf[i][j]!=0:
                tietze_list.append(j+1)
    return [tietze_list+['nan'],X[1]]

# Tests whether two words are equal, i.e., the function compares the normal forms computed as above of X and Y
def are_equal(X,Y):
    if nx.is_isomorphic(X[1][0],Y[1][0])==False:
        return 'The elements do not belong to isomorphic groups!'
    else:
        XX=normal_form(X)
        YY=normal_form(Y)
        if len(XX)!=len(YY):
            return 'No'
        else:
            for i in range(min(len(XX),len(YY))):
                if XX[i]==YY[i]:
                    continue
                else:
                    return 'No'
                    break
        return 'Yes'

# Print the normal form of an element X with respect to the generators y_i. y has to be defined as a variable before calling the function
def print_normal_form(y,X):
    nf=normal_form(X)
    s=''
    
    x=str(y)
    for i in range(len(nf)):
        for j in range(len(nf[i])):
            if nf[i][j]!=0:
                if nf[i][j]!=1:
                    s+=x+str(j+1)+"^"+str(nf[i][j])+str(" ")
                else:
                    s+=x+str(j+1)+str(" ")
    return s

# Return the product of elements in the list x. The list x is made of elements of the group
def dot(x):
    X=x[0][0][:-1]
    for i in range(1,len(x)):
        if nx.is_isomorphic(x[0][1][0],x[i][1][0])==False:
            return 'The elements belong to non-isomorphic groups!'
            break
        X=X+x[i][0][:-1]
    return [X+['nan'],x[0][1]]
