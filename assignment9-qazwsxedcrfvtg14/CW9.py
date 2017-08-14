from __future__ import division, print_function
from visual import *
scene.width = 1000
scene.height = 700
R = 0.08
Q = 50e-9
Nq = 25   
dtheta = 2*pi/Nq
dQ = Q/Nq

e0=8.854187817e-12
pi=3.14159

scalefactor = 3e-7 # you may need to change this

## source charges
ring(pos=vector(0,0,0), radius=R, color=color.red, thickness=0.005)
sources = []
angles=arange(0,2*pi,dtheta)
for theta in angles:
    a = sphere(pos=vector(0, R*cos(theta),R*sin(theta)), radius=0.007, color=color.cyan, q=dQ)
    sources.append(a)

pointer=[]

for i in range(17):
    pos=(vector(-2*R,0,0)*(16-i)+vector(2*R,0,0)*(i))/16
    #axis=vector(sin(pi*2/12*i)*1e-11,cos(pi*2/12*i)*1e-11,0)
    a=vector(0,0,0)
    for j in range(Nq):
        a+=1/(4*pi*e0)*dQ/(pos-sources[j].pos).mag2*(pos-sources[j].pos).norm()
    axis=a*scalefactor
    pointer.append(arrow(pos=pos, axis=axis))

for i in range(7):
    pos=(vector(R/2,0,-2*R)*(6-i)+vector(R/2,0,2*R)*(i))/6
    #axis=vector(sin(pi*2/12*i)*1e-11,cos(pi*2/12*i)*1e-11,0)
    a=vector(0,0,0)
    for j in range(Nq):
        a+=1/(4*pi*e0)*dQ/(pos-sources[j].pos).mag2*(pos-sources[j].pos).norm()
    axis=a*scalefactor
    pointer.append(arrow(pos=pos, axis=axis))

#p3 =sphere(pos=vector(2*R,0,0), radius=0.007, color=color.green, make_trail=True )

mouseposition = scene.waitfor('click').pos 
#p3.pos=mouseposition
#p3 =sphere(pos=vector(2*R,0,0), radius=0.007, color=color.green, make_trail=True )
p3 =sphere(pos=mouseposition, radius=0.007, color=color.green, make_trail=True )


p3.v=vector(0,0,0)
p3.m=9.1e-31
q3=-1.6e-19
k=8.9875517873681764e9
dt=1e-11
t=0
while True:
    t+=dt
    rate(300)
    f=vector(0,0,0)
    for j in range(Nq):
        f+=k*dQ*q3/(p3.pos-sources[j].pos).mag2*(p3.pos-sources[j].pos).norm()
    p3.v+=f*dt/p3.m
    p3.pos+=p3.v*dt
    



## outer loop picks observation location 
#for Ea in obs: 
#    Enet = vector(0,0,0) 
    ## inner loop goes through all source charges
    #for scharge in sources: 
        ## add code to calculate contribution of this source charge 
        ## add this to net field at this location 
    #Ea.axis = Enet*scalefactor
