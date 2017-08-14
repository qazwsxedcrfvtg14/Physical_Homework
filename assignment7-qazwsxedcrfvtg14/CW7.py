from visual import *
from random import random
from visual.graph import *
import numpy as np

#random.seed(7122)

N = 50                              # number of the He atoms
#N=2
L = ((24.4E-3/(6E23))*N)**(1/3.0)/2 # side length of the cubic container box

m, size = 4E-3/6E23, 310E-12        # He atoms mass, radius are made 10 times bigger for easiear collision
radius=np.array([size for v in range(N)])
rad=np.array([[-L+size,-L+size,-L+size] for v in range(N)])
L_size=L-size                       # L - size, used many times in the program
k, T = 1.38E-23, 298.0              # k = Boltzmann constant, and T = temperature in unit K

t, dt = 0, 0.5E-13                  # dt = 0.5E-13 second
#t, dt = 0, 0.1E-15                  # dt = 0.5E-13 second
vrms = (3*k*T/m)**0.5               # root mean square velocity at T
atoms = []                          # list to contain the 50 atoms

# histogram initialization, more on this see http://vpython.org/contents/docs/graph.html
deltav = 100.
vdist = gdisplay(x=800, y=0, ymax = N*deltav/1000.,width=500, height=300, xtitle='v', ytitle='dN')
theory = gcurve(color=color.cyan) # plot the theoretical speed distribution
dv = 10.
for v in arange(0.,3001.+dv,dv):
    theory.plot(pos=(v,(deltav/dv)*N*4.*pi*((m/(2.*pi*k*T))**1.5)*exp((-0.5*m*v**2)/(k*T))*v**2*dv))
observation = ghistogram(bins=arange(0.,3000.,deltav),accumulate=1, average=1, color=color.red)


# initialization of display, setting up for the random position distribution and random velocity direction of atoms
scene = display(width=800, height=800,background=(0.2,0.2,0))
container = box(length = 2*L, height = 2*L, width = 2*L, opacity=0.2, color = color.yellow )



# Initialize atom position and velocity
for i in range(N):
    position = vector(-L_size+2*L_size*random(),-L_size+2*L_size*random(),-L_size+2*L_size*random())
    if i== N-1:
        atom = sphere(pos=position, radius = size, color=color.yellow, make_trail = True, retain = 600)
    else:
        atom = sphere(pos=position, radius = size, color=(random(), random(), random()))
        theta, phi = pi*random(), 2*pi*random()
    atom.m, atom.v = m, vector(vrms*sin(theta)*cos(phi), vrms*sin(theta)*sin(phi), vrms*cos(theta))
    atoms.append(atom)
#atoms[0].v=[-vrms,0,0.0000000001]
#atoms[0].pos=[size*2,size*2,size*2]
#atoms[1].v=[vrms,vrms,vrms]
#atoms[1].pos=[-size*2,-size*2,-size*2]

    
def vcollision(a,b,vel,pos):
    '''
    Function to find the velocities of atoms after each collision
    '''
    #v1prime = vel[a] - 2 * m/(m+m) *(pos[a]-pos[b]) * dot (vel[a]-vel[b], pos[a]-pos[b]) / abs(pos[a]-pos[b])**2
    #v2prime = vel[b] - 2 * m/(m+m) *(pos[b]-pos[a]) * dot (vel[b]-vel[a], pos[b]-pos[a]) / abs(pos[b]-pos[a])**2
    v1=vector(vel[a])
    v2=vector(vel[b])
    p1=vector(pos[a])
    p2=vector(pos[b])
    v1prime = v1 - (p1-p2) * dot (v1-v2, p1-p2) / abs(p1-p2)**2
    v2prime = v2 - (p2-p1) * dot (v2-v1, p2-p1) / abs(p2-p1)**2
    vel[a]=v1prime
    vel[b]=v2prime
    #print(vel[a],vel[b],v1prime,v2prime,vel[a]-vel[b],pos[a]-pos[b],dot (vel[a]-vel[b], pos[a]-pos[b]) , abs(pos[a]-pos[b])**2)
    #return v1prime, v2prime
    
def back_time(a,b,vel,pos):
    l=0.
    r=10.
    for i in range(100):
        mid=(l+r)/2.
        if sum((pos[a]-vel[a]*mid*dt-pos[b]+vel[b]*mid*dt)**2)<4*(size**2):
            l=mid
        else:
            r=mid
    #print(sum((pos[a]-pos[b])**2),sum((pos[a]-vel[a]*r*dt-pos[b]+vel[b]*r*dt)**2)-4*(size**2))
    return r

print("PreRun!")
ddt=dt*50
for run in range(500):
    rate(10000)
    #calculate new positions for all atoms and plot histogram
    pos=[]
    vel=[]
    for i in range(N):
        pos.append(atoms[i].pos)
        vel.append(atoms[i].v)
    pos=np.array(pos)
    vel=np.array(vel)
    pos+=ddt*vel
    outside = less_equal(pos,-L+size) # walls closest to origin 
    r = pos-pos[:,newaxis] # all pairs of atom-to-atom vectors
    rmag = sqrt(sum(square(r),-1)) # atom-to-atom scalar distances
    hit = less_equal(rmag,radius+radius[:,newaxis])-identity(N)
    hitlist = sort(nonzero(hit.flat)[0]).tolist() # i,j encoded as i*N+j
    for i in hitlist:
        a=i/N
        b=i%N
        if a<b:
            if(dot(vel[b]-vel[a],pos[b]-pos[a])>0):
                continue
            vcollision(a,b,vel,pos)
    v1 = vel*outside
    vel = vel-v1+abs(v1) # force p component inward 
    p1=pos*outside
    pos=rad*outside*2-p1+pos*logical_not(outside)
    outside = greater_equal(pos,L-size) # walls farther from origin 
    v1 = vel*outside
    vel = vel-v1-abs(v1) # force p component inward
    p1=pos*outside
    pos=-rad*outside*2-p1+pos*logical_not(outside)
    for i in range(N):
        atoms[i].pos=pos[i]
        atoms[i].v=vel[i]

print("PreRun Finish!")


fp_cnt=0.
fp_tot=0.
fp_now_t=np.array([0 for i in range(N)])
dt_cnt=0.
dap = 0
avp=0
pr=N*m*(vrms**2.)*1./3./((2.*L)**3.)
fr=((2.*L)**3.)/(1.41421*3.14159*((2.*size)**2.)*N)
print("theory pressure:", pr)
print("theory mean free path:", fr)
now_t=0
while True:
    dt_cnt+=1
    if dt_cnt%1000==0:
        if fp_cnt==0:
            print(0)
        else:
            #print(fp_tot/fp_cnt*N,fp_tot/fp_cnt,fp_tot,fp_cnt,fr)
            av_Force = dap / now_t
            avp = av_Force / ( 4.* L * L)
            print("average pressure: ", avp)
            print("average mean free path: ", fp_tot/fp_cnt)
            #print(fp_tot/fp_cnt)
    now_t += dt
    rate(10000)
    #calculate new positions for all atoms and plot histogram
    v=[]
    pos=[]
    vel=[]
    for i in range(N):
        pos.append(atoms[i].pos)
        vel.append(atoms[i].v)
    
    pos=np.array(pos)
    vel=np.array(vel)
    #for i in range(N):
    pos+=dt*vel
    fp_now_t=fp_now_t+dt
    #while True:
    for ti in range(1):
        brk=True
        outside = less_equal(pos,-L+size) # walls closest to origin 
        r = pos-pos[:,newaxis] # all pairs of atom-to-atom vectors
        rmag = sqrt(sum(square(r),-1)) # atom-to-atom scalar distances
        hit = less_equal(rmag,radius+radius[:,newaxis])-identity(N)
        hitlist = sort(nonzero(hit.flat)[0]).tolist() # i,j encoded as i*N+j
        for i in hitlist:
            a=i/N
            b=i%N
            if a<b:
                if(dot(vel[b]-vel[a],pos[b]-pos[a])>0):
                    continue
                t=back_time(a,b,vel,pos)
                #if(b==N-1):
                fp_now_t[a]-=dt*t
                fp_now_t[b]-=dt*t
                fp_tot+=fp_now_t[a]*sqrt(sum(vel[a]**2))
                fp_tot+=fp_now_t[b]*sqrt(sum(vel[b]**2))
                fp_cnt+=2

                pos[a]-=dt*vel[a]*t
                pos[b]-=dt*vel[b]*t
                #print(sum((pos[a]-pos[b])**2),sum((pos[a]-pos[b])**2)-(2*size)**2,(2*size)**2)
                vcollision(a,b,vel,pos)
                
                pos[a]+=dt*vel[a]*t
                pos[b]+=dt*vel[b]*t
                #if(b==N-1):
                fp_now_t[a]=dt*t
                fp_now_t[b]=dt*t
                
        v1 = vel*outside
        vel = vel-v1+abs(v1) # force p component inward 
        

        dap += 2. * m * sum(vel*outside*np.array([[True,False,False] for lop in range(N)]))


        p1=pos*outside
        pos=rad*outside*2-p1+pos*logical_not(outside)
        #pos=rad*outside
        outside = greater_equal(pos,L-size) # walls farther from origin 
        
        v1 = vel*outside
        vel = vel-v1-abs(v1) # force p component inward
        p1=pos*outside
        pos=-rad*outside*2-p1+pos*logical_not(outside)
    for i in range(N):
        atoms[i].pos=pos[i]
        atoms[i].v=vel[i]
    #print(dap,t)
    for i in range(N):
        #### calculate new positions for atoms

        v.append(mag(atoms[i].v))
    observation.plot(data=v)

        #### find collisions between atoms, and then handle their collisions

        #### find collisions between atoms and walls, and then handle their collision
        #### calculate the momentum transferred to the walls and obtain the pressure
        #### print the averaged pressure on the walls every 1000*dt
