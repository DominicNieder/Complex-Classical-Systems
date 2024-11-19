import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray  # trying to be typesafe

class Particle:
    def __init__(
        self,
        number: int,
        mass: float | int,
        dimensions: int
    ):
        """
        Issue: No need to seperate x, y axis -> update all functions to take array([[x],[y]])
        """
        self.n: int = number
        self.m: float = float(mass)
        self.d: int = dimensions
        self.r: NDArray[np.float64]= np.zeros((dimensions,number))
        self.vx: NDArray[np.float64]= np.zeros((dimensions,number))
        self.a1: NDArray[np.float64]= np.zeros((dimensions,number))
        self.a2: NDArray[np.float64]= np.zeros((dimensions,number))
        
    def __repr__(self):
        return f"This is a particle at {self.x[0]}, {self.y[0]} with v={self.vx[0]},{self.vy[0]} interacting by Lennard-Jones potential."
    
    ### initiating particles on a grid
    def initiate_positions_on_grid(
            self,
            box: tuple[float|int, float|int]
    ) -> None:
        """
        Initializes 2D positions on a Grid with even spacing.

        Parameters:
                
            n_particles:(int)
                    -> number of particles

            box:(tuple[number,number]) 
                    -> box size (in [nm])
        """
        grid_sections = int(np.ceil(np.sqrt(self.n)))  # find the number of colums & rows

        # even spacing
        x_spacing = box[0]/grid_sections 
        y_spacing = box[1]/grid_sections
        # makes grid coordinates
        x, y= np.meshgrid(
            np.arange(grid_sections) * x_spacing, 
            np.arange(grid_sections) * y_spacing
        )
        x,y= x.flatten()[:self.n], y.flatten()[:self.n]
        self.r= np.linalg.matrix_transpose(np.array([x,y]))   
        pass
    
    ### add a random jitter to the initial positions
    def jitter_posiiotns(
            self,
            jitter_strength: float
    ) -> None:
        """
        Add a jitter with jitter strength (uniform sampling) to the 2D Particle positions.

        Parameters:

            jitter_strength: (float)
                    ->  scaling factor to the random sampling. Initial sampling intervall [-0.5, 0.5).
        """
        self.r+= (np.random.rand(2,self.n)-0.5)* jitter_strength
        pass

    ### initiate velocities accord to the Boltzmann distribution
    def initiate_velocities_temperature(
            self,
            temperature: float|int,
            mass: float|int
    )-> None:
        """
        Initiation method of Particle class: 2D velocity vectors to a Boltzmann distribution to a given temperature.

        Parameters:

            n_particles: (int)
                    -> number of particles in the system

            temperature: (number)
                    -> temperature of the system (in [K])

            mass: (number)
                    -> mass of the particles  (in [kg])

        """
        k_B = 1.380649*10**(-23)  # Boltzmann constant (J/K)

        # Initialize velocities for each particle in two dimensions
        std_dev = np.sqrt(k_B * temperature / mass)
        velocities = np.random.normal(0, std_dev, (self.n, 2))

        # Remove net momentum to ensure the center of mass is stationary
        mean_velocity = np.mean(velocities, axis=0)
        velocities -= mean_velocity

        self.v= velocities    
        pass


class Data:
    def __init__(
            self,
            particle: Particle,
            simulation_time: int,
            positions: NDArray,
            velocities: NDArray,
            potential: function
            ):
        positions = np.zeros((particle.n,))
        pass

    def run_time_analysis(self):
        pass