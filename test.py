import sys
from descartes import PolygonPatch
import matplotlib.pyplot as plt
import alphashape

points = [(0., 0.), (0., 1.), (1., 1.), (1., 0.),
          (0.5, 0.25), (0.5, 0.75), (0.25, 0.5), (0.75, 0.5)]

fig, ax = plt.subplots()
ax.scatter(*zip(*points))
plt.show()

alpha_shape = alphashape.alphashape(points, 0.)

fig, ax = plt.subplots()
ax.scatter(*zip(*points))
ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2))
plt.show()

alpha_shape = alphashape.alphashape(points, 2.0)

fig, ax = plt.subplots()
ax.scatter(*zip(*points))
ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2))
plt.show()

alpha_shape = alphashape.alphashape(points, 3.5)

fig, ax = plt.subplots()
ax.scatter(*zip(*points))
ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2))
plt.show()

alpha_shape = alphashape.alphashape(points, 5.0)
print(alpha_shape)