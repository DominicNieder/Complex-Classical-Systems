import numpy as np
import matplotlib.pyplot as plt

### particle object
class Particle:
    def __init__(self, r, x, y, vx, vy):
        self.r=r
        self.x=x
        self.y=y
        self.vx=vx
        self.vy=vy
        
    def __repr__(self):
        return str("This is a particle at...")
     

### determin distance
def distance(p1:Particle,p2:Particle):
    d=(p1.x-p2.x)+(p1.y-p2.y)

### scattering of two hard balls -> iter_vel
def scatter(pos1, pos2, vel1, vel2):
    d=pos1-pos2
    e_r= d/np.sqrt(d[0]**2+d[1]**2)
    print('distance of scattering balls',e_r)

def periodic_boundry(x,y,box):  # impl. periodic boundry condition
    x=x % box[0]
    y=y % box[1]


def iter_pos(particle,box,dt):  # iterates position
    particle.x+=particle.vx*dt
    particle.y+=particle.vy*dt
    periodic_boundry(particle.x,particle.y,box)

def iter_vel(particle):  # iterates velocity
    pass 


def gauss(obj,box,dt):  # gaussian integration sceem
    iter_pos(obj,box,dt)
    iter_vel(obj,dt)
        
def sv(obj,sv_pos, sv_vel, time_step):  # saving positon and velocities
    sv_pos[0][time_step], sv_pos[1][time_step]= obj.x, obj.y
    sv_vel[0][time_step], sv_vel[1][time_step]= (obj.vx, obj.vy)



### Simulation
if __name__ == "__main__":
     
    ### defining particles with attributes and time
    n=1  # number of particles
    r=0.5  # nm
    m=1  # kg
    
    dt=1.0  # arbitrary unit
    T=2000  # number of time steps

    v0=0.5  #nm/ initial velocity
    init_vel=np.random.rand(n)*2*np.pi  # random vel direction between (0,2pi]
    print(init_vel)

    data_traj=np.zeros((4,T,n))

    box = (20.0,20.0)
    
    ### array of all objects
    objs=np.array(
        [Particle(r,
        np.random.randint(1,box[0]-1,1),
        np.random.randint(1,box[1]-1,1),
        v0*np.cos(init_vel[i]),
        v0*np.sin(init_vel[i])
        ) for i in range(1,n)])

    #### Simulation
    for t in range(T):
        gauss(objs,box)   # time integration
        sv(objs, sv_pos, sv_vel, t)    # saving  position & velocity 
    
    ### Plotting trajectory
    print(data_traj)
    plt.plot(sv_pos)
    plt.savefig("trajectory-1P_4.png")
    