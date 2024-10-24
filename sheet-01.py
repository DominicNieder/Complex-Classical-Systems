import numpy as np
import matplotlib.pyplot as plt

from matplotlib import colors
from matplotlib.ticker import PercentFormatter

class Particle:
    def __init__(self, r:float|int, x:float|int, y:float|int, vx:float|int, vy:float|int):
        self.r:float|int=r
        self.x:float|int=x
        self.y:float|int=y
        self.vx:float|int=vx
        self.vy:float|int=vy
        
    def __repr__(self):
         return(str("This is a particle at %0.2f, %0.2f with v=%0.2f,%0.2f" % (self.x,self.y,self.vx,self.vy)) )

### shift to S->S' of which p1'=(box[0]/2, box[1]/2)
def system_shift(p1:Particle,box:tuple[float|int]) -> tuple[float|int]:
    (dx,dy) = (p1.x-box[0]/2 , p1.y-box[1]/2)
    return(dx,dy)

### vec(p2->p1)
def rel_koordiantes(p1:Particle,p2:Particle,box:tuple[float|int]) -> tuple[float|int]:  
    dx, dy = system_shift(p1,box)  # S->S': S' is a box with p1' at its center
    #  print("10,10=",p1.x+dx,p1.y+dy)
    #  print("p2'=",(p2.x+dx)%box[0], (p2.y+dy)%box[1])
    #  print("rel koordiantes=", (p1.x+dx - (p2.x+dx)%box[0] , p1.y+dy - (p2.y+dy)%box[1]))
    return(p1.x +dx - (p2.x-dx)%box[0] , p1.y +dy - (p2.y-dy)%box[1])


### calculating distributions
### velocity-squared distribution
def collect_vabs(vx:list[list[float|int]],vy:list[list[float|int]],interval:tuple[int]) -> list[float]:
    v_abs=np.zeors(len(vx)+(interval[1]-interval[0]))
    for i in range(len(vx)):
        for j in range(interval[0],interval[1]):
            v_abs[i+j-interval[0]]=(np.sqrt(vx[i][j]**2+vy[i][j]**2))
    return v_abs

### velocity y distribution
def collect_v_1D(vx:list[float|int],interval:tuple[int]) -> list[float]:
    l=np.zeros(len(vx)+(interval[1]-interval[0]))
    for i in range(len(vx)):
        for j in range(interval[0],interval[1]):
            l[i+j-interval[0]]=(vx[i][j])
    return l
   
### Maxwell Boltzmann distribution 
def mw_boltzmann_distri(v_intervall:tuple[float|int], T:float|int)->tuple[list[float]]: 
    v=np.linspace(v_intervall[0],v_intervall[1]) 
    p=4*np.pi/(2*np.pi*T)**(3/2)*v**2 * np.exp(v**2/2/T)  # k_B=1, m=1
    return (p,v)

def run():
    log="Initiating\n"
    box:tuple[float|int]=(20,20)
    n:int=50
    r:float|int=1.0
    v0:float|int=0.5
    angles=np.random.rand(n)*2*np.pi
    par:list[Particle]=[Particle(
            r,
            np.random.randint(1,box[0]-1,1)[0],
            np.random.randint(1,box[1]-1,1)[0],
            v0*np.cos(angles[i]),
            v0*np.sin(angles[i])
            )
        for i in range(n)]   

    T:int=2000
    dt:float=1.0
    log = log + "n_particles=" +str(n)+"\ntimesteps="+str(T)+"\ntimestep="+str(dt)+"\nbox="+str(box)+"\n" 

    x_data=np.zeros((n,T))
    y_data=np.zeros((n,T))
    vx_data=np.zeros((n,T))
    vy_data=np.zeros((n,T))
    ### Simulation
    for t in range(T):
        ### iterate all particles
        for i in range(n):
            #  print("iter i=",i)
            ### collision
            for j in range(i+1,n):  # iteration over all others
                #  print("iter j=", j)
                d_vec:tuple[float|int] = rel_koordiantes(par[i],par[j], box)  # rel. distance vector j -> i
                d_abs = d_vec[0]**2+d_vec[1]**2  # distance squared
                #  print("distance vector",d_vec[0],d_vec[1])
                #  print("distance=",d_abs)
                #### collision condition -> kinematics
                if np.sqrt(d_abs) <= par[i].r+par[j].r:
                    d_unit:tuple[float|int] = d_vec/np.sqrt(d_abs)  # collision axis
                    #  print("!collision!: Axis",d_unit)
                    dvix , dviy= -((par[i].vx-par[j].vx)*d_unit[0]+ (par[i].vy-par[j].vy)*d_unit[1]) *d_unit
                    #  print("dv=",(dvix, dviy))
                    par[i].vx , par[i].vy = (par[i].vx+dvix) , ( par[i].vy+dviy)
                    par[j].vx , par[j].vy = (par[j].vy-dviy) , (par[j].vy-dviy)
            #  print("particle",i ," velocities:",par[i].vx,par[i].vy)

            par[i].x = (par[i].x+par[i].vx*dt) % box[0] 
            par[i].y = (par[i].y+par[i].vy*dt) % box[1]
            x_data[i][t]=par[i].x
            y_data[i][t]=par[i].y 
            vx_data[i][t]=par[i].vx
            vy_data[i][t]=par[i].vy
    log=log+"Simulation complete!\n"
    ### plotting kinetic energy
    # e_kin = 1/2* sum([par[k].vx**2+par[k].vy**2 for k in ])

    ### plotting historgam
    n_bins = 20

    # Generate velocity distributionsdistributions 
    dist1 = collect_v_1D(vx_data,(0,200))
    dist2 = collect_v_1D(vx_data,(1800,2000))
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)

    # We can set the number of bins with the *bins* keyword argument.
    axs[0].hist(dist1, bins=n_bins)
    axs[1].hist(dist2, bins=n_bins)
    # Saving figure
    plt.savefig("histograms-vx.png")
    return log

if __name__ == "__main__":
    print(run())
