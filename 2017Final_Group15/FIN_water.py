from visual import *
from random import random
from visual.graph import *
import numpy as np

R=0.3

scene.width = 1000
scene.height = 700
cir = shapes.circle(radius=R+R*0.1)
circ = shapes.circle(radius=R)
wash = extrusion(pos=[(0, 0, R), (0, 0, -R)], shape=cir-circ, color=color.white)
bottom = extrusion(pos=[(0,0,-R), (0, 0, -R)], shape=circ, color=color.white)


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
KQQ=0.00001
KQQ_D=0.03
R-=size
rot_speed=100
water=10
air=35
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
        bal.append(sphere(pos=pos, radius = size, color=(1, 0, 1.0*i/(N-1))))
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
g=vector(0,-9.8,0)*100
t=0
while True:
    t+=dt
    rate(1/dt)
    for i in range(N*M):
        pos[i]=bal[i].pos
    will_pos=pos+v*dt
    for i in range(N*M):
        if (will_pos[i][0]**2+will_pos[i][1]**2)>=R**2 and (will_pos[i][0]**2+will_pos[i][1]**2) > (pos[i][0]**2+pos[i][1]**2):
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
            rot=vec*rot_speed*max(-np.dot(-pos_norm,g),0)/g.mag * sin(t/0.25*3.14)
            #print(t)
            #if int(t/0.25) %2 == 0:
            #    rot*=-1
            #print(rot)
            c+=rot
            if(np.sum(rot**2)**0.5>0.002):
                col=bal[i].color
                #print(col)
                if col[1]<1.0:
                    col=(col[0],col[1]+0.002,col[2])
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
        for j in range(len(nxt[i])):
            now=np.sum((pos[i]-pos[nxt[i][j]])**2)**0.5
            #if(now > stor[i][j]):
            F=-K*(now-stor[i][j])
            if F>100:
                F=100
            if F<-100:
                F=-100
            F=F*((pos[i]-pos[nxt[i][j]])/now)
            if(now > stor[i][j]*0.9):
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
        if (np.sum(v[i]**2))**0.5>water and v[i][1]<-R/2:
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
