import numpy as np
import matplotlib.pyplot as plt

import UBKLibpp
import UBKLib


def getMagneticAmplitude(x, y):
    return UBKLibpp.calculateB(x, y, 1)


potential_levels = np.linspace(-100, 100, 100)

equipotentials = UBKLib.generate_equipotentials(
    UBKLib.cross_tail_potential, 
    potential_levels,
    x_bounds=(-15.0, 15.0),
    y_bounds=(-15.0, 15.0),
    resolution=1000
)

w0_points = UBKLib.contour_w0_points(equipotentials, getMagneticAmplitude)

upper_points = UBKLib.parse_upper_contour_w0_points(w0_points)
lower_points = UBKLib.parse_lower_contour_w0_points(w0_points)

plt.figure(figsize=(10, 10))

for level_idx, level_contours in enumerate(equipotentials):
    for contour in level_contours:
        plt.plot(contour[:, 0], contour[:, 1], 'b-', linewidth=0.5, alpha=0.5)

for point in upper_points:
    plt.plot(point.x, point.y, 'x', color='red')

for point in lower_points:
    plt.plot(point.x, point.y, 'x', color='blue')

plt.xlabel('X')
plt.ylabel('Y')
plt.title('Equipotential Contours')
plt.grid(True, alpha=0.3)
plt.axis('equal')  # Keep aspect ratio square
plt.show()
