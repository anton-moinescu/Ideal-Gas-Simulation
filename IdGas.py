#This program simulates the an ideal gas (or a system with properties similar to one), where all particles start with random velocity
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
V=2
r=0.1
T=1000
dt=0.1
F=T/dt
bins = 60

class Particle:
    def __init__(self, x, y, vx, vy):
        self.x=x
        self.y=y
        self.vx=vx
        self.vy=vy
##Wall collision- condition: if the distance between the center of a particle and the wall is smaller than the radius and the particle is moving towards wall
    def lateral_wall_collision(self):
        if((L-self.x<=r and self.vx>0) or (self.x<r and self.vx<0 )):
            self.vx=-self.vx

    def upper_lower_wall_collision(self):
        if ((l-self.y<=r and self.vy>0)or (self.y<=r and self.vy<0)):
            self.vy=-self.vy

    def distance(self, other_particle):
        return np.sqrt((self.x-other_particle.x)*(self.x-other_particle.x)+(self.y-other_particle.y)*(self.y-other_particle.y))

    def interm_dot_prod(self, other_particle):
        vrel=np.array([self.vx-other_particle.vx,self.vy-other_particle.vy])
        drelvect=np.array([self.x-other_particle.x,self.y-other_particle.y])
        dot_prod=np.dot(vrel,drelvect)
        drelscal=self.distance(other_particle)
        interm=dot_prod/(drelscal*drelscal)
        return interm

    def V1x_after_col(self, other_particle):
        interm=self.interm_dot_prod(other_particle)
        vx1after=self.vx-interm*(self.x-other_particle.x)
        self.vx= vx1after

    def V1y_after_col(self, other_particle):
        interm=self.interm_dot_prod(other_particle)
        vy1after=self.vy-interm*(self.y-other_particle.y)
        self.vy= vy1after

    def V2x_after_col(self, other_particle):
        interm=self.interm_dot_prod(other_particle)
        vx2after=self.vx+interm*(self.x-other_particle.x)
        other_particle.vx= vx2after

    def V2y_after_col(self, other_particle):
        interm=self.interm_dot_prod(other_particle)
        vy2after=self.vy+interm*(self.y-other_particle.y)
        other_particlevy= vy2after

    def condcol(self, other_particle):
        d=self.distance (other_particle)
        if(d<=2*r):
            if(np.sqrt(((self.x+self.vx*dt)-(other_particle.x+other_particle.vx*dt))**2+
                    (self.y+self.vy*dt-(other_particle.y+other_particle.vy*dt))**2)<d):
            #print("Collision detected")
                return True
        else:
            return False


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

#defining the particles array
particles_list=[]
for j in range (0,N):
    particles_list.append(Particle(x[j], y[j], vx[j], vy[j]))

def position_and_velocity_updates(particles_list, v):
    for i in range(0, N):
        particles_list[i].lateral_wall_collision()
        particles_list[i].upper_lower_wall_collision()

        for j in range(i+1,N):
            if(particles_list[i].condcol(particles_list[j])):
                particles_list[i].V1x_after_col(particles_list[j])
                particles_list[i].V1y_after_col(particles_list[j])
                particles_list[i].V2x_after_col(particles_list[j])
                particles_list[i].V2y_after_col(particles_list[j])
                v[i]=np.sqrt(particles_list[i].vx**2+particles_list[i].vy**2)
                v[j]=np.sqrt(particles_list[j].vx**2+particles_list[j].vy**2)

        particles_list[i].x=particles_list[i].x+particles_list[i].vx*dt
        particles_list[i].y=particles_list[i].y+particles_list[i].vy*dt

    #print(x[50])
    return  v, particles_list


fig, ax = plt.subplots()



n, bins, patches = ax.hist(v, bins=bins, range=(0, np.max(v)), color='green', alpha=0.7)
ax.set_xlim(0, V)
ax.set_ylim(0, N//5)
ax.set_xlabel('Speed')
ax.set_ylabel('Number of Particles')
ax.set_title('Speed Distribution')

def animate(frame):
    global v, particles_list
    # Update the positions of the particles
    v, particles_list = position_and_velocity_updates(particles_list, v)
    ax.cla()  # Clear the previous histogram
    ax.hist(v, bins=bins, range=(0, np.max(v)), color='green', alpha=0.7)
    ax.set_xlim(0, V)
    ax.set_ylim(0, N//5)
    #ax.set_xlabel('Speed')
    #ax.set_ylabel('Number of Particles')
    #ax.set_title('Speed Distribution')
    # Return the modified plot elements
    return ()

anim = FuncAnimation(fig, animate, frames=1000, interval=50, blit=True)

plt.show()

