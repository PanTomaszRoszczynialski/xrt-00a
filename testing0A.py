# -*- coding: utf-8 -*-
"""
Created on Fri Feb 06 17:58:11 2015

@author: Vladimir Putin
"""
import numpy as np
import matplotlib.pyplot as plt


layers = 12
r0 = 0.1
wall = 0.02
rSample = 100

xList = []
alphaList= []
for n in range(layers):
    ms = range(n)
    i6 = range(6)
    for i in i6:
        for m in ms:
            x = 2*(r0+wall) * (n**2 + m**2 - n*m)**0.5
            xList.append(x)
            alpha = np.arcsin( x/rSample )
            alphaList.append(alpha)
            