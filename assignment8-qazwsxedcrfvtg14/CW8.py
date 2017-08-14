from __future__ import division, print_function
from visual import *
from visual.graph import *
scene.width = scene.height = 800

## constants
oofpez = 9e9 # stands for One Over Four Pi Epsilon-Zero
qe = 1.6e-19  
s = 4e-11   
scalefactor = 1e-23 # you may need to change this

p1 = sphere(pos=vector(0,s/2,0), radius=1e-11, color=color.red)
q1 = qe
p2 = sphere(pos=vector(0,-s/2,0), radius=1e-11, color=color.blue)
q2 = -qe

p3 =sphere(pos=vector(-3*s,0,0), radius=3e-12, color=color.green, make_trail=True )
p3.v=vector(0,0,0)
p3.m=1.672621898e-27
q3=qe

## for plotting

egraphs=gdisplay(x=600, width=600)   ## move graph so it's not on top of scene
#scene.width = scene.height = 600
Ug = gcurve(color=color.yellow)
Kg = gcurve(color=color.cyan)
KUg = gcurve(color=color.magenta)
pointer=[]
e0=8.854187817e-12
pi=3.14159
for i in range(12):
    pos=(p1.pos+p2.pos)/2+vector(sin(pi*2/12*i)*s,cos(pi*2/12*i)*s,0)*1.2
    #axis=vector(sin(pi*2/12*i)*1e-11,cos(pi*2/12*i)*1e-11,0)
    a1=1/(4*pi*e0)*q1/(pos-p1.pos).mag2*(pos-p1.pos).norm()
    a2=1/(4*pi*e0)*q2/(pos-p2.pos).mag2*(pos-p2.pos).norm()
    axis=(a1+a2)*scalefactor
    pointer.append(arrow(pos=pos, axis=axis))

for i in range(12):
    pos=(p1.pos+p2.pos)/2+vector(0,sin(pi*2/12*i)*s,cos(pi*2/12*i)*s)*1.2
    #axis=vector(sin(pi*2/12*i)*1e-11,0,cos(pi*2/12*i)*1e-11)
    a1=1/(4*pi*e0)*q1/(pos-p1.pos).mag2*(pos-p1.pos).norm()
    a2=1/(4*pi*e0)*q2/(pos-p2.pos).mag2*(pos-p2.pos).norm()
    axis=(a1+a2)*scalefactor
    pointer.append(arrow(pos=pos, axis=axis))

k=8.9875517873681764e9
dt=1e-17
t=0
while True:
    t+=dt
    rate(300)
    f1=k*q1*q3/(p3.pos-p1.pos).mag2*(p3.pos-p1.pos).norm()
    f2=k*q2*q3/(p3.pos-p2.pos).mag2*(p3.pos-p2.pos).norm()
    p3.v+=(f1+f2)*dt/p3.m
    p3.pos+=p3.v*dt
    K=1./2*p3.m*p3.v.mag2
    Kg.plot(pos=(t,K))
    U=k*q1*q3/(p3.pos-p1.pos).mag+k*q2*q3/(p3.pos-p2.pos).mag
    Ug.plot(pos=(t,U))
    KUg.plot(pos=(t,K+U))
