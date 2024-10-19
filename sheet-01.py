import numpy as np
import matplotlib.pyplot as plt

class Particle:
    def __init__(self, r, x, y, vx, vy):
        self.r=r
        self.x=x
        self.y=y
        self.vx=vx
        self.vy=vy
        
    def __repr__(self):
         return(str("This is a particle at %0.2f, %0.2f with v=%0.2f,%0.2f" % (self.x,self.y,self.vx,self.vy)) )

# p2 -> p1
def distance(p1,p2):  
    return(p1.x-p2.x , p1.y-p2.y)

def collision(p1, p2):
    d21=distance(p1,p2)  # distance pointing to p1
    e21=d21/np.sqrt(d21[0]**2 + d21[1]**2)
    theta=np.arctan(d21[1]/d21[0])
    v_rel=(p2.vx-p1.vx, p2.vy-p1.vx)
    vx1=np.cos(theta)*p2.vx
    vy1=np.cos(theta)*p2.vy
    return(vx1, vy1)


def run():
    box=(20,20)
    n=2
    r=1.0
    v0=0.5
    # angles=np.random.rand(n)*2*np.pi
    angles=np.array([0, np.pi/2])
    par=[Particle(
            r,
            10.0+i,
            10.0,
            # np.random.randint(1,box[0]-1,1)[0],
            # np.random.randint(1,box[1]-1,1)[0],
            v0*np.cos(angles[i]),
            v0*np.sin(angles[i])
            )
        for i in range(n)]    

    T=50
    x_data=np.zeros((n,T))
    y_data=np.zeros((n,T))
    vx_data=np.zeros((n,T))
    vy_data=np.zeros((n,T))
    print(type(x_data))
    print(type(x_data[0]))
    x_data[0][0]=1.0
    print(x_data[0][0])
    dt=1.0

    for t in range(T):
        for i in range(n):  # iteration
            for j in range(i,n):
                dvx_ij = par[j].vx - par[i].vx
                dvy_ij = par[j].vy - par[i].vy
                par[i].vx=par[i].vx
                par[i].vy=par[i].vy

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
    

if __name__ == "__main__":
    run()