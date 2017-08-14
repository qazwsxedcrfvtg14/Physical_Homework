from __future__ import division, print_function
import visual
import math
from visual.graph import *
from visual.factorial import *
def combin(a,b):
    ans=1
    for i in range(a,a-b,-1):
        ans*=i
    for i in range(1,b+1):
        ans/=i
    return ans
def microstates(q,N):
   """ returns the number of microstates Omega for a system with 
       q energy quanta and N oscillators 
       INPUT: 
          q: number of energy quanta
          N: number of oscilators
   """
   return math.factorial(q+N-1)//math.factorial(q)//math.factorial(N-1)

def entropy(q,N):
   """ returns the entropy S for a system with 
       q energy quanta and N oscillators 
       INPUT: 
          q: number of energy quanta
          N: number of oscilators
   """
   
   return 1.38064852e-23*math.log(microstates(q,N))

def temperature(q,N,w):
   """ returns the temperature T for a system with 
       q energy quanta and N oscillators of energy quantum \Delta E=\hbar w
       INPUT: 
          q: number of energy quanta
          N: number of oscilators
          w: oscillator angular frequency
   """
   dE=1.054571800e-34*w
   dS=entropy(q+1,N)-entropy(q,N)
   return dE/dS

def H(a,b):
    return visual.factorial.combin(a+b-1,b)
#print(combin(5,3))

from visual import *
Ntotal = 700 ## total number of oscillators 
N1 = 400 ## number of oscillators in object 1 
N2 = Ntotal-N1 ## number of oscillators in object 2 
qtotal = 100 ## total quanta of energy shared among all the oscillators

waygraph = gvbars(delta=0.7, color=color.red) # to make vertical bar graph 

q1 = 0 ## start with no quanta of energy in object 1
max_pos=0
max_dis=0
while q1 <= qtotal: ## for each possible value of energy in object 1 
    q2 = qtotal-q1 ## number of quanta of energy in object 2 
    ways1 = H(N1,q1) ## Calculate number of ways to arrange q1 quanta in object 1 
    ways2 = H(N2,q2) ## Calculate number of ways to arrange q2 quanta in object 2: 
    waygraph.plot( pos=(q1,ways1*ways2) ) ## Plot number of ways to arrange energy in both objects 
    if ways1*ways2 > max_pos:
        max_pos=ways1*ways2
        max_dis=q1
    q1 = q1+1
print(max_dis,qtotal-max_dis)
max_pos/=2
q1 = 0 ## start with no quanta of energy in object 1
half_pos=float("Inf")
max_dis=0
while q1 <= qtotal: ## for each possible value of energy in object 1 
    q2 = qtotal-q1 ## number of quanta of energy in object 2 
    ways1 = H(N1,q1) ## Calculate number of ways to arrange q1 quanta in object 1 
    ways2 = H(N2,q2) ## Calculate number of ways to arrange q2 quanta in object 2: 
    if math.fabs(ways1*ways2-max_pos) < half_pos:
        half_pos=math.fabs(ways1*ways2-max_pos)
        max_dis=q1
    q1 = q1+1
print(max_dis,qtotal-max_dis)

graph = gdisplay(x=0, y=400, height=300)
waygraph1 = gcurve(color=color.red)
waygraph2 = gcurve(color=color.blue)
waygraph3 = gcurve(color=color.cyan)

q1 = 0 ## start with no quanta of energy in object 1
max_q1=0
max_omega=0
while q1 <= qtotal: ## for each possible value of energy in object 1 
    q2 = qtotal-q1 ## number of quanta of energy in object 2 
    waygraph1.plot( pos=(q1,math.log(microstates(q1,N1))) ) ## Plot number of ways to arrange energy in both objects 
    waygraph2.plot( pos=(q1,math.log(microstates(q2,N2))) ) ## Plot number of ways to arrange energy in both objects 
    omega=math.log(microstates(q1,N1)*microstates(q2,N2))
    waygraph3.plot( pos=(q1,omega) ) ## Plot number of ways to arrange energy in both objects 
    if omega > max_omega:
        max_omega=omega
        max_q1=q1
    q1 = q1+1
print(max_omega,max_q1)

k=16*4
N=105
q=150
m=27/1000/(6.02e23)
w=sqrt(k/m)
#print(q,N,w,entropy(q,N),microstates(q,N))
print(temperature(q,N,w))