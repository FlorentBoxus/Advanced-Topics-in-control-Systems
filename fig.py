# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 12:09:47 2023

@author: Florent Boxus
"""
import numpy as np
import matplotlib.pyplot as plt

# Load data
t = np.loadtxt(r'C:\Users\Florent Boxus\Documents\master1\Advanced topics\Projet\t.txt')
x = np.loadtxt(r'C:\Users\Florent Boxus\Documents\master1\Advanced topics\Projet\x.txt')
V = np.loadtxt(r'C:\Users\Florent Boxus\Documents\master1\Advanced topics\Projet\V.txt')

# Ensure x and t are 1D arrays
x = np.squeeze(x)
t = np.squeeze(t)

# Create meshgrid from x and t
X, T = np.meshgrid(x, t)

# Plot pcolormesh
plt.pcolormesh(np.transpose(V), cmap='viridis', shading='auto')

# Add colorbar and labels
plt.colorbar(label='V')
plt.xlabel('x')
plt.ylabel('t')
plt.title('Color Plot of V as a Function of x and t')

# Show the plot
plt.show()
print(V[:,2000])
