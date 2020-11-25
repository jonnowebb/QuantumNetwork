# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 12:05:26 2020

@author: jonno
"""

import networkx as nx
import numpy as np
import pylab as plt
import scipy

#X=Qubit/vertice
#Z=Entanglement/Edge
#Y=Phase change
#I=Identity

HilbertSpace=2
#----Setting up Pauli matrices----
I=np.identity(HilbertSpace)
X=np.array([[0,1],[1,0]])
Y=np.array([[0,-1j],[1j,0]])
Z=np.array([[1,0],[0,-1]])

#----Testing Pauli's----
PauliSquaredX=X**2
PauliSquaredY=Y**2
PauliSquaredZ=Z**2

PauliXY=np.dot(X,Y)
PauliYZ=np.dot(Y,Z)
PauliZX=np.dot(Z,X)
PauliYX=np.dot(Y,X)
PauliZY=np.dot(Z,Y)
PauliXZ=np.dot(X,Z)

#----Setting up check matrix----
#Goal: Input any graph from Bristol paper
#and get a tableau representation out
#Let's start with a 3 qubit set up
#The check matrix will be 6 long by 3 tall
Qubits=[1,2,3,4,5]
#CheckMatrix=np.zeros(shape=(len(Qubits),len(Qubits)*2))
#First half is the X matrix, second is the Z matrix
#Each row is a generator, so for n qubits there is n generators,
#and each column is what the generator is made of

#Lets initilise our qubits,
G=nx.Graph()
for i in Qubits:
    G.add_node(i)
#And now lets create the entanglement we want to start with
#Entangle=[[1,2],[1,3]]
#G.add_edges_from(Entangle)
#nx.draw(G,with_labels=True)
#This is the numerics approach, one we need to move away from
#as we cannot use any of this information for good other than
#making a graph we want plotted.

Edges=G.edges()
OurEdges=[]
for u,v in Edges: OurEdges.append([u,v])


#All our qubits are initilised as |+>
h=np.array([[1],[0]])
v=np.array([[0],[1]])
d=1/np.sqrt(2)*(h+v)
a=1/np.sqrt(2)*(h-v)
r=1/np.sqrt(2)*(h+1j*v)
l=1/np.sqrt(2)*(h-1j*v)
InitialQubits=np.ones(len(Qubits))


Generators=np.array((['X','X','X','X','X'],['I','Z','Z','Z','Z'],['Z','I','Z','Z','Z'],['Z','Z','I','Z','Z'],['Z','Z','Z','I','Z']))
if len(Generators)!=len(Qubits):
    raise TypeError("Generator list mismatches qubit list")
XMatrix = Generators == 'X'
#This is for helping sort out the Standard Form further down
XMatrix2 = Generators == 'X'
ZMatrix = Generators == 'Z'
CheckMatrix=np.concatenate((XMatrix,ZMatrix),axis=1)
#We need to check a valid input generator
Zeros=np.zeros(shape=(len(Qubits),len(Qubits)))
Identity=np.identity(len(Qubits))
Tri1=np.concatenate((Zeros, Identity), axis=1)
Tri2=np.concatenate((Identity, Zeros), axis=1)
Tri=np.concatenate((Tri1,Tri2))
#Check if input is valid
IsValid=np.dot(CheckMatrix,np.dot(Tri,np.transpose(CheckMatrix))) %2
if np.count_nonzero(IsValid)==0:
    #Apply Hadamards
    for i in range(0,len(Qubits)):
        if np.any(XMatrix[i])==False:
            XMatrix[:,i]=ZMatrix[:,i]
            ZMatrix[:,i]=XMatrix2[:,i]
        if XMatrix[-1][-1]==False:
            XMatrix[:,-1]=ZMatrix[:,-1]
            ZMatrix[:,-1]=XMatrix2[:,-1]
    CheckMatrix=np.concatenate((XMatrix,ZMatrix),axis=1)
    #Put into standard form
    StandardForm=scipy.linalg.solve(XMatrix,CheckMatrix) %2
else:
    #Try another generator sequence if this happens
    print('Invalid generator')
#We can read off this Standard Form matrix to produce the graph
StdFormGraph=nx.Graph()
if sum(np.diag(XMatrix))==sum(np.diag(Identity)):
    for i in range(0,len(np.diag(XMatrix))):
        StdFormGraph.add_node(i)
        if ZMatrix[i,:].any()==1:
            for j in range(0,len(np.where(ZMatrix[i,:]==True)[0][:])):
                StdFormGraph.add_edge(np.where(ZMatrix[i,:]==True)[0][j],np.where(XMatrix[i]==True)[0][0])
nx.draw(StdFormGraph,with_labels=True)
