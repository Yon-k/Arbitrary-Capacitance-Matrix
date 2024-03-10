import numpy as np

class Charge(object):
    def __init__(self, value, radial_distance, coordinates, particle_index):
        self.value = value
        self.radial_distance = radial_distance
        self.coordinates = coordinates
        self.particle_index = particle_index