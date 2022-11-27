import numpy as np

class Ball:
    def __init__(self, radius, mass, name):
        self.m = mass
        self.r = radius
        self.gamma = np.arcsin(np.sqrt(np.power((self.r + 1.5e-3), 2) - pow(7e-3, 2)) / (self.r + 1.5e-3))
        self.d = self.r * np.sin(self.gamma)
        self.mu_s = np.tan(np.radians(12)) * np.sin(self.gamma)
        self.mu_k = self.mu_s / 1.1     # approx
        self.name = name
