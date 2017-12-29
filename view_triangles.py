#!/usr/bin/env python

import sys

triangles = []

for line in sys.stdin:
    t = line.strip().split()
    v0 = [0, 0, 0]
    v1 = [0, 0, 0]
    v2 = [0, 0, 0]
    (v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]) = [float(f) for f in t]
    triangles.append([v0, v1, v2])

minv = [1000000, 1000000, 1000000]
maxv = [-1000000, -1000000, -1000000]

for t in triangles:
    for v in t:
        for c in (0, 1, 2):
            minv[c] = min(minv[c], v[c])
            maxv[c] = max(maxv[c], v[c])

# Get the center of the model
center = [0, 0, 0]
for c in (0, 1, 2):
    center[c] = (minv[c] + maxv[c]) / 2

max_dimension = 0
for c in (0, 1, 2):
    max_dimension = max(max_dimension, (maxv[c] - minv[c]));

# The camera looks down -Z towards the model
# FOV is always 45 degrees, camera is moved so that the model fits.
camera = [center[0], center[1], maxv[2] + max_dimension / .9]

for t in triangles:
    for v in t:
        for c in (0, 1, 2):
            v[c] -= camera[c]
    print t[0][0], t[0][1], t[0][2], t[1][0], t[1][1], t[1][2], t[2][0], t[2][1], t[2][2]


