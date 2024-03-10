from particle import Particle
from utils import neighbors, calculate_coefficients, get_final_coefficients

if __name__== "__main__":
    particles = [Particle([0, 0], 1), Particle([1, 1], 1)]
    for i, particle in enumerate(particles):
        neighbors_list = neighbors(particle, particles)
        get_final_coefficients(particles, 1, [0, 0], neighbors_list)
    for particle in particles:
        print(particle.final_coefficient())