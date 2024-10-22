import numpy as np
import matplotlib.pyplot as plt

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
    (dx,dy) = (box[0]/2-p1.x , box[1]/2-p1.y)
    return(dx,dy)

### vec(p2->p1)
def rel_koordiantes(p1:Particle,p2:Particle,box:tuple[float|int]) -> tuple[float|int]:  
    dx, dy = system_shift(p1,box)  # S->S': S' is a box with p1' at its center
    print("10,10=",p1.x+dx,p1.y+dy)
    print("p2'=",(p2.x+dx)%box[0], (p2.y+dy)%box[1])
    return(p1.x - (p2.x+dx)%box[0] , p1.y - (p2.y+dy)%box[1])

def run():
    box:tuple[float|int]=(20,20)
    n:int=2
    r:float|int=1.0
    v0:float|int=0.5
    angles=np.random.rand(n)*2*np.pi
    #  par:list[Particle]=[Particle(
    #          r,
    #          np.random.randint(1,box[0]-1,1)[0],
    #          np.random.randint(1,box[1]-1,1)[0],
    #          v0*np.cos(angles[i]),
    #          v0*np.sin(angles[i])
    #          )
    #      for i in range(n)]   
    par:list[Particle]=[Particle(
        1.0,
        0.0,
        0.0,
        2.0,
        0.0
    ), Particle(
        1.0,
        0.0,
        18.0,
        0.0,
        0.0
    )] 

    T:int=50
    x_data=np.zeros((n,T))
    y_data=np.zeros((n,T))
    vx_data=np.zeros((n,T))
    vy_data=np.zeros((n,T))
    dt:float=1.0
    ### Simulation
    for t in range(T):
        ### iterate all particles
        for i in range(n):
            print("iter i=",i)
            ### collision
            for j in range(i+1,n):  # iteration over all others
                print("iter j=", j)
                d_vec:tuple[float|int] = rel_koordiantes(par[i],par[j], box)  # rel. distance vector j -> i
                d_abs = d_vec[0]**2+d_vec[1]**2  # distance squared
                print("distance=",d_abs)
                #### collision condition -> kinematics
                if np.sqrt(d_abs) <= par[i].r+par[j].r:
                    d_unit:tuple[float|int] = d_vec/np.sqrt(d_abs)  # collision axis
                    print("!collision!: Axis",d_unit)
                    dvix , dviy= -((par[i].vx-par[j].vx)*d_unit[0]+ (par[i].vy-par[j].vy)*d_unit[1]) *d_unit
                    print("dv=",(dvix, dviy))
                    par[i].vx , par[i].vy = (par[i].vx+dvix) , ( par[i].vy+dviy)
                    par[j].vx , par[j].vy = (par[j].vy-dviy) , (par[j].vy-dviy)
            print("particle",i ," velocities:",par[i].vx,par[i].vy)
            par[i].x = (par[i].x+par[i].vx*dt) % box[0] 
            par[i].y = (par[i].y+par[i].vy*dt) % box[1]
            x_data[i][t]=par[i].x
            y_data[i][t]=par[i].y 
            vx_data[i][t]=par[i].vx
            vy_data[i][t]=par[i].vy

    #plt.plot(x_data[0]-x_data[1],y_data[0]-y_data[1],label="relative")
    plt.plot(x_data[0],y_data[0],label="p1")
    plt.plot(x_data[1],y_data[1],label="p2")
    plt.xlabel("x [nm]")
    plt.ylabel("y [nm]")
    plt.legend()
    plt.savefig("relative-position-2p-system_5.png")
    

if __name__ == "__main__":
    run()
