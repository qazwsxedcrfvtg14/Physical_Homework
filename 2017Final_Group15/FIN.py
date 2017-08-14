from random import random

import numpy as np
from visual import *
from visual.graph import *

R=0.3

scene.width = 800
scene.height = 800
#cir = shapes.circle(radius=R+R*0.1)
#circ = shapes.circle(radius=R)
#wash = extrusion(pos=[(0, 0, R), (0, 0, -R)], shape=cir-circ, color=color.white)
#bottom = extrusion(pos=[(0,0,-R), (0, 0, -R)], shape=circ, color=color.white)


N=11
M=5
W=0.2
L=0.1
size=0.01
bal=[]
nxt=[]
stri=[]
stor=[]
K=100000
KQQ=0.00000
KQQ_D=10.03
rot_speed=100
water=15
air=30
chg_T=0.25
g=vector(0,-9.8,0)*100
str_max=1.2
water_level=-0.7*R


cir = shapes.circle(radius=R+R*0.1)
circ = shapes.circle(radius=R)
#semicircle = paths.arc(radius=R, angle1=5*pi/4, angle2=7*pi/4)
boxbox = shapes.rectangle(pos=(0,0), width=3*R)
bb=extrusion(pos=[(0, 0, R), (0, 0, -R)], shape=boxbox-cir, color=color.white)
wash = extrusion(pos=[(0, 0, R), (0, 0, -R)], shape=cir-circ, color=color.white)
bottom = extrusion(pos=[(0,0,-R), (0, 0, -R)], shape=circ, color=color.white)
#water = extrusion(pos=[(0, 0, R), (0, 0, -R)], shape=semicircle, color=color.cyan)
water_box = box(pos=(0, (water_level-R)/2, 0), length=R*2, height=R+water_level, width=2*R-0.01, color=color.cyan, opacity=0.2)
pointer_rad=0.
pointer_scaler=0.0015
pointer = arrow(pos=(cos(pointer_rad)*R,sin(pointer_rad)*R,R),axis=(cos(pointer_rad)*R/3,sin(pointer_rad)*R/3,0),color=color.black)
RR=R
R-=size

"""
for i in range(N):
    for j in range(M):
        pos=vector( (-W*(N-1-i)+W*(i))/(N-1),0, (-L*(M-1-j)+L*(j))/(N-1) )
        #bal.append(sphere(pos=pos, radius = size, color=(random(), random(), random())))
        bal.append(sphere(pos=pos, radius = size, color=(0, 0, 255)))
        nxt.append([])
        if i>0:
            nxt[-1].append((i-1)*M+j)
        if j>0:
            nxt[-1].append(i*M+j-1)
        if i<N-1:
            nxt[-1].append((i+1)*M+j)
        if j<M-1:
            nxt[-1].append(i*M+j+1)
"""
for i in range(N):
    for j in range(M):
        if i % 2 == 0:
            pos=vector( (-W*(N-1-i)+W*(i))/(N-1),0, (-L*(M-1-j)+L*(j))/(M-1) )
        else:
            pos=vector( (-W*(N-1-i)+W*(i))/(N-1),0, (-L*(M-1-j)+L*(j))/(M-1)+L/(M-1) )
        #bal.append(sphere(pos=pos, radius = size, color=(random(), random(), random())))
        #bal.append(sphere(pos=pos, radius = size, color=(1, 0, 1.0*i/(N-1))))
        bal.append(sphere(pos=pos, radius = size, color=(1, 0, 0.0*i/(N-1))))
        nxt.append([])
        if i>0:
            nxt[-1].append((i-1)*M+j)
        if j>0:
            nxt[-1].append(i*M+j-1)
        if i<N-1:
            nxt[-1].append((i+1)*M+j)
        if j<M-1:
            nxt[-1].append(i*M+j+1)
        if i % 2 == 0 :
            if i>0 and j<M-1:
                nxt[-1].append((i-1)*M+j+1)
            if i<N-1 and j>0:
                nxt[-1].append((i+1)*M+j-1)
        else:
            if i>0 and j>0:
                nxt[-1].append((i-1)*M+j-1)
            if i<N-1 and j<M-1:
                nxt[-1].append((i+1)*M+j+1)
for i in range(N*M):
    stri.append([])
    stor.append([])
    for j in nxt[i]:
        stri[-1].append(helix(pos=bal[i].pos, axis=bal[j].pos-bal[i].pos, radius=size/2,thickness =size/5))
        stor[-1].append((bal[j].pos-bal[i].pos).mag)
        if i>j:
            stri[-1][-1].visible=False
v=[]
for i in range(N*M):
    v.append(vector(0,0,0))
v=np.array(v)
pos=[]
for i in range(N*M):
    pos.append(vector(0,0,0))
pos=np.array(pos)
dt=0.001
t=0
while True:
    t+=dt
    rate(1/dt)
    for i in range(N*M):
        pos[i]=bal[i].pos
    will_pos=pos+v*dt

    now_rot_speed=rot_speed*sin(t/chg_T*3.1415926)
    pointer_rad+=pointer_scaler*now_rot_speed
    pointer.pos=(cos(pointer_rad)*RR,sin(pointer_rad)*RR,RR)
    pointer.axis=(cos(pointer_rad)*RR/3,sin(pointer_rad)*RR/3,0)

    for i in range(N*M):
        if (will_pos[i][0]**2+will_pos[i][1]**2)>=R**2:
            v[i]=v[i]
            a=vector(v[i])
            will_pos[i][2]=0
            b=vector(-will_pos[i])
            pos_norm=will_pos[i]/np.sum(will_pos[i]**2)**0.5
            """c=(a-2*(a-a.dot(b)/b.mag*b.norm()))
            c_pos=np.dot(c,pos_norm)*pos_norm
            """

            c=a.dot(b)/b.mag*pos_norm
            #print(a,c,a.dot(b)/b.mag)
            c=a+c*1.01

            #c=c-c_pos+c_pos*0.2
            vec=np.cross(-will_pos[i],[0,0,1])
            vec=vec/np.sum(vec**2)**0.5
            #print(vec)
            rot=vec*max(-np.dot(-pos_norm,g),0)/g.mag * now_rot_speed
            #print(t)
            #if int(t/0.25) %2 == 0:
            #    rot*=-1
            #print(rot)
            c+=rot
            col=bal[i].color
            #print(col)
            if col[1]<1.0:
                col=(col[0],col[1]+np.sum(rot**2)**0.5/10000,col[2])
            bal[i].color=col
            #print(i,-pos_norm,g,-np.dot(-pos_norm,g),max(-np.dot(-pos_norm,g),0))
            v[i][0]=c.x
            v[i][1]=c.y
            v[i][2]=c.z
            #dot=v[i][0]*will_pos[i][0]+v[i][1]*will_pos[i][1]+v[i][2]*will_pos[i][2]
            #dot/
        if will_pos[i][2]>R:
            v[i][2]=-abs(v[i][2])
        if will_pos[i][2]<-R:
            v[i][2]=abs(v[i][2])
    for i in range(N*M):
        if (np.sum(v[i]**2))**0.5>water and pos[i][1]<water_level:
            v[i]=v[i]/((np.sum(v[i]**2))**0.5)*water
        if (np.sum(v[i]**2))**0.5>air:
            v[i]=v[i]/((np.sum(v[i]**2))**0.5)*air
    for i in range(N*M):
        for j in range(len(nxt[i])):
            now=np.sum((pos[i]-pos[nxt[i][j]])**2)**0.5
            #if(now > stor[i][j]):
            F=-K*(now-stor[i][j])
            if F>100:
                F=100
            if F<-100:
                F=-100
            F=F*((pos[i]-pos[nxt[i][j]])/now)
            if(now > stor[i][j]*str_max):
                ve=v[i].dot((pos[i]-pos[nxt[i][j]]))/np.sum((pos[i]-pos[nxt[i][j]])**2)
                if ve > 0:
                    v[i]-=ve*(pos[i]-pos[nxt[i][j]])*0.2
            v[i]+=F*dt
    for i in range(N*M):
        for j in range(N*M):
            if i!=j:
                r=np.sum((pos[i]-pos[j])**2)**0.5
                if r < KQQ_D:
                    v[i]+=KQQ*((pos[i]-pos[j])/r)/r/r
    #print(v)
    for i in range(N*M):
        if (np.sum(v[i]**2))**0.5>water and pos[i][1]<water_level:
            v[i]=v[i]/((np.sum(v[i]**2))**0.5)*water
        if (np.sum(v[i]**2))**0.5>air:
            v[i]=v[i]/((np.sum(v[i]**2))**0.5)*air
    v=v+dt*g
    pos=pos+v*dt
    for i in range(N*M):
        if (pos[i][0]**2+pos[i][1]**2)>=R**2:
            pos_norm=pos[i]/((pos[i][0]**2+pos[i][1]**2)**0.5)*R
            pos[i][0]=pos_norm[0]
            pos[i][1]=pos_norm[1]
            #print(pos[i])
    for i in range(N*M):
        bal[i].pos=pos[i]
    for i in range(N*M):
        for j in range(len(nxt[i])):
            stri[i][j].pos=bal[i].pos
            stri[i][j].axis=bal[nxt[i][j]].pos-bal[i].pos

quit()
