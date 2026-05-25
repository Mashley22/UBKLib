import UBKLibpp
import UBKLib
import matplotlib.pyplot as plt

K_VALUES = [1, 10, 100, 1000]


def POTENTIAL(x, y):
    return UBKLib.volland_stern_potential(x, y, kp=0)


grid = UBKLib.Grid(
    K_VALUES,
    POTENTIAL,
    UBKLibpp.calculateB,
    resolution=100
)

if __name__ == "__main__":
    grid.calc_potential_grid()
    grid.calc_magnetic_amp_grid_parallel()
    plt.plot(grid.magnetic_amp_grids[0].flatten(), grid.potential_grid.flatten(), '^', alpha=0.1)
    plt.show()
