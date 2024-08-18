#This program simulates the an ideal gas (or a system with properties similar to one), where all particles start with the same speed
#Initial data: box dimension;
#Initial data about particles: nr of particles (or density), initial total energy (or initial speed), radius of particles

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

"""
L=input('Lenght of the box ')
L=int(L)
l=input('Width of the box')
l=int(l)
N=input('Number of particles')
N=int(N)
v=input('Initial speed')
v=int(v)
r=input('Radius of particles')
r=int(r)
T=input ('Total time of simulation') ##note: make the simulation stop when a certain distribution is reached
T=int(T)
dt=input('Time increment')
dt=float(dt)
"""

L=10
l=10
N=1000
V=1
r=0.1
T=1000
dt=0.1
F=T/dt
bins = 30
##Wall collision- condition: if the distance between the center of a particle and the wall is smaller than the radius and the particle is moving towards certain wall + velocity condition

#1-upper wall
#2-right wall
#3-lower wall
#4-left wall

#Wall collision condition
def cond1(y,vy):
    if (l-y<=r and vy>0):
        return 1
    else:
        return 0

def cond3(y,vy):
    if(y<=r and vy<0):
        return 1
    else:
        return 0
def cond2(x, vx):
    if(L-x<=r and vx>0):
        return 1
    else:
        return 0
def cond4(x, vx):
    if(x<r and vx<0):
        return 1
    else:
        return 0

#Wall collision consequences
def after1or3(vy):
    return(-vy)

def after2or4(vx):
    return(-vx)

def distance(x1,y1,x2,y2):
    return np.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))

#particle collision condition
def condcol(x1, y1, x2, y2, vx1, vy1, vx2, vy2):
    d=distance(x1,y1,x2,y2)
    futured=distance(x1+vx1*dt,y1+vy1*dt,x2+vx2*dt,y2+vy2*dt)
    if(d<=2*r and futured<d):
        #print("Collision detected")
        return 1
    else:
        return 0

def modrelV(v1x,v1y,v2x,v2y):
    v=math.sqrt((v1x-v2x)*(v1x-v2x)+(v1y-v2y)*(v1y-v2y))
    return v

##used for the velocities after collision:
def interm_dot_prod(x1,y1,x2,y2,vx1,vy1,vx2,vy2):
    vrel=np.array([vx1-vx2,vy1-vy2])
    drelvect=np.array([x1-x2,y1-y2])
    dot_prod=np.dot(vrel,drelvect)
    drelscal=distance(x1,y1,x2,y2)
    interm=dot_prod/(drelscal*drelscal)
    return interm

def V1x_after_col(x1,y1,x2,y2,vx1,vy1,vx2,vy2):
    interm=interm_dot_prod(x1,y1,x2,y2,vx1,vy1,vx2,vy2)
    vx1after=vx1-interm*(x1-x2)
    return vx1after

def V1y_after_col(x1,y1,x2,y2,vx1,vy1,vx2,vy2):
    interm=interm_dot_prod(x1,y1,x2,y2,vx1,vy1,vx2,vy2)
    vy1after=vy1-interm*(y1-y2)
    return vy1after

def V2x_after_col(x1,y1,x2,y2,vx1,vy1,vx2,vy2):
    interm=interm_dot_prod(x1,y1,x2,y2,vx1,vy1,vx2,vy2)
    vx2after=vx2+interm*(x1-x2)
    return vx2after

def V2y_after_col(x1,y1,x2,y2,vx1,vy1,vx2,vy2):
    interm=interm_dot_prod(x1,y1,x2,y2,vx1,vy1,vx2,vy2)
    vy2after=vy2+interm*(y1-y2)
    return vy2after

#for intial velocity
deg=360*np.random.random(N)
v=V*np.random.random(N)
vx=v*np.cos(np.radians(deg))
vy=v*np.sin(np.radians(deg))
#print(vx)
#print(vy)
#print(v)

#for initial positions:
x=r+(L-r)*np.random.random(N)
y=r+(l-r)*np.random.random(N)
def position_updates(x, y, vx, vy):
    for i in range(0, N):
        if(cond1(y[i],vy[i])==1):
            vy[i]=-vy[i]
        if(cond3(y[i],vy[i])==1):
            vy[i]=-vy[i]
        if(cond2(x[i], vx[i])==1):
            vx[i]=-vx[i]
        if(cond4(x[i], vx[i])==1):
            vx[i]=-vx[i]

        for j in range(i+1,N):
            if(condcol(x[i], y[i], x[j], y[j], vx[i], vy[i], vx[j], vy[j])==1):
                vx[i]=V1x_after_col(x[i],y[i],x[j],y[j],vx[i],vy[i],vx[j],vy[j])
                vy[i]=V1y_after_col(x[i],y[i],x[j],y[j],vx[i],vy[i],vx[j],vy[j])
                vx[j]=V2x_after_col(x[i],y[i],x[j],y[j],vx[i],vy[i],vx[j],vy[j])
                vy[j]=V2y_after_col(x[i],y[i],x[j],y[j],vx[i],vy[i],vx[j],vy[j])
                v[i]=np.sqrt(vx[i]*vx[i]+vy[i]*vy[i])
                v[j]=np.sqrt(vx[j]*vx[j]+vy[j]*vy[j])


    x=x+vx*dt
    y=y+vy*dt
    #print(x[50])
    return x, y, vx, vy, v


fig, ax = plt.subplots()


##de revenit si fct 2 vers- una in care calc dinainte tot, una cum face gpt
n, bins, patches = ax.hist(v, bins=bins, range=(0, np.max(v)), color='green', alpha=0.7)
ax.set_xlim(0, V)
ax.set_ylim(0, N//2)
ax.set_xlabel('Speed')
ax.set_ylabel('Number of Particles')
ax.set_title('Speed Distribution')

def animate(frame):
    global x, y, vx, vy, v
    # Update the positions of the particles
    x, y, vx, vy, v = position_updates(x, y, vx, vy)
    ax.cla()  # Clear the previous histogram
    ax.hist(v, bins=bins, range=(0, np.max(v)), color='green', alpha=0.7)
    ax.set_xlim(0, V)
    ax.set_ylim(0, N//2)
    #ax.set_xlabel('Speed')
    #ax.set_ylabel('Number of Particles')
    #ax.set_title('Speed Distribution')
    # Return the modified plot elements
    return ()

anim = FuncAnimation(fig, animate, frames=1000, interval=50, blit=True)

plt.show()

