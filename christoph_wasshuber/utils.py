import numpy as np
from particle import Particle

def neighbors(particle_i, particles):
    neighbor_particles = []
    for j, particle_j in enumerate(particles):
        if np.linalg.norm(particle_i.coordinates - particle_j.coordinates) <= particle_i.radius + particle_j.radius and particle_i != particle_j:
            neighbor_particles.append(particle_j)
    return neighbor_particles

def calculate_coefficients(particles):
    for particle in particles:
        particle.influence_coefficients.append(particle.influence_coefficient())

def get_final_coefficients(particles, radius, coords, neighbors):
    for i, particle in enumerate(particles):
        particle.calculate_prime_charges(Charge(1, radius, coords, i), neighbors)
        calculate_coefficients(particles)