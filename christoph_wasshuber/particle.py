import numpy as np
from charge import Charge

class Particle(object):
    def __init__(self, coordinates, radius):
        self.coordinates = np.array(coordinates)
        self.radius = radius
        self.charges = []
        self.influence_coefficients = []

    def assign_charge(self, charge):
        self.charges.append(charge)

    def influence_coefficient(self):
        return sum(charge.value for charge in self.charges)

    def final_coefficient(self):
        return self.influence_coefficients[-1] if self.influence_coefficients else None

    def calculate_prime_charges(self, charge0, particles):
        prime_charges = []
        for j, particle in enumerate(particles):
            if self == particle:
                continue
            else:
                d = np.linalg.norm(self.coordinates - particle.coordinates)
                Q_prime = -(charge0.value * particle.radius) / (d - charge0.radial_distance)
                radial_distance_charge = (particle.radius ** 2) / (d - charge0.radial_distance)
                coordinatesQ_prime = particle.coordinates - ((particle.radius / d) ** 2) * (particle.coordinates - self.coordinates)
                prime_charge = Charge(Q_prime, radial_distance_charge, coordinatesQ_prime, j)
                particle.assign_charge(prime_charge)
                prime_charges.append(prime_charge)
        return prime_charges