import numpy as np
import matplotlib.pyplot as plt

class Particle:
    def __init__(self, r:float, x:float, y:float, vx:float, vy:float):
        self.r:float=r
        self.x:float=x
        self.y:float=y
        self.vx:float=vx
        self.vy:float=vy
        
    def __repr__(self):
         return(str("This is a particle at %0.2f, %0.2f with v=%0.2f,%0.2f" % (self.x,self.y,self.vx,self.vy)) )

# p2 ->p1
def distance(p1:Particle,p2:Particle) -> tuple[float]:  
    return((p1.x-p2.x) , (p1.y-p2.y))


def run():
    box=(20,20)
    n:int=2
    r:float=1.0
    v0:float=0.5
    # angles=np.random.rand(n)*2*np.pi
    angles:list[float]=np.array([0, np.pi/2])
    par:list[Particle]=[Particle(
            r,
            10.0+i,
            10.0,
            # np.random.randint(1,box[0]-1,1)[0],
            # np.random.randint(1,box[1]-1,1)[0],
            v0*np.cos(angles[i]),
            v0*np.sin(angles[i])
            )
        for i in range(n)]    

    T:int=50
    x_data=np.zeros((n,T))
    y_data=np.zeros((n,T))
    vx_data=np.zeros((n,T))
    vy_data=np.zeros((n,T))
    dt:float=1.0

    for t in range(T):
        for i in range(n):  # iteration over all particles
            for j in range(i+1,n):  # iteration over all others
                d_vec = distance(par[i],par[j])  # distance vec j -> i
                d_abs = d_vec[0]**2+d_vec[1]**2  # distance
                d_unit = d_vec/np.sqrt(d_abs)  # dist direct
                if d_abs <= 4*par[0].r**2:  # condition for collision
                    dvix,dviy= -((par[i].vx-par[j].vx)*d_unit[0]+ (par[i].vy-par[j].vy)*d_unit[1]) *d_unit
                    par[i].vx, par[i].vy= par[i].vx+dvix , par[i].vy+dviy
                    par[j].vx, par[j].vy= par[j].vy-dviy , par[j].vy-dviy

            par[i].x = (par[i].x+par[i].vx*dt) % box[0] 
            par[i].y = (par[i].y+par[i].vy*dt) % box[1]
            x_data[i][t]=par[i].x
            y_data[i][t]=par[i].y 
            vx_data[i][t]=par[i].vx
            vy_data[i][t]=par[i].vy

    plt.plot(x_data[0]-x_data[1],y_data[0]-y_data[1],label="relative")
    plt.plot(x_data[0],y_data[0],label="p1")
    plt.plot(x_data[1],y_data[1],label="p2")
    plt.xlabel("x [nm]")
    plt.ylabel("y [nm]")
    plt.legend()
    plt.savefig("relative-position-2p-system_5.png")
    

def diffGeo(c):
    x=np.linspace(-10,10,200)
    y=np.linspace(-10,10,200)
    x_plot=[]
    y_plot=[]
    for i in x:
        for j in y:
            if i**2-j**2==c:
                x_plot.append(i)
                y_plot.append(j)



    plt.plot(x_plot,y_plot)
    plt.savefig("diffgeo2.png")


if __name__ == "__main__":
    run()
