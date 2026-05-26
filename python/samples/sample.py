import UBKLibpp
import UBKLib
import matplotlib.pyplot as plt
import argparse
import multiprocessing as mp
import numpy as np

K_VALUES = [1, 10, 100, 1000]
RESOLUTION = 1000
THREAD_COUNT = mp.cpu_count()

ELECTRON_REST_MASS = 510.99895069
PROTON_REST_MASS = 939272.0894329

MU = 0.224
CHARGE = -1
K_IDX = 0

grid = None


def HAMILTONIAN(b, u):
    return UBKLib.relatavistic_hamiltonian(b, u, MU, CHARGE, ELECTRON_REST_MASS)


def POTENTIAL(x, y):
    return UBKLib.volland_stern_potential(x, y, kp=0)


def UB_TRAJECTORY(b, total_energy):
    return UBKLib.relativistic_ub_trajectory(
        b,
        MU,
        CHARGE,
        ELECTRON_REST_MASS,
        CHARGE * total_energy
    )


def gen_data():
    global grid
    grid = UBKLib.FullGrid(
        K_VALUES,
        POTENTIAL,
        UBKLibpp.calculateB,
        resolution=RESOLUTION
    )
    grid.calc_potential_grid()
    grid.calc_magnetic_amp_grid_parallel()


def save_data(filepath):
    grid.save(filepath)


def load_data(filepath):
    global grid
    grid = UBKLib.FullGrid.load(filepath)


def plot():
    grid.calc_cross_products_grids()
    grid.calc_hamiltonian_grids(HAMILTONIAN)
    lcds = grid.find_lcds_contours(1, 50)
    w0_contours = grid.find_cross_product_zeros()[K_IDX]
    plt.plot(
        grid.magnetic_amp_grids[K_IDX].discrete.flatten(),
        grid.potential_grid.discrete.flatten(),
        '^', alpha=0.01, color='grey'
    )
    for contour in w0_contours:
        plt.plot([x.B for x in contour], [x.U for x in contour])
    
    b_space = np.array([
        grid.magnetic_amp_grids[K_IDX](lcds[K_IDX][1][i], lcds[K_IDX][2][i])
        for i in range(len(lcds[K_IDX][1]))
    ])
    u_space = np.array([
        grid.potential_grid(lcds[K_IDX][1][i], lcds[K_IDX][2][i])
        for i in range(len(lcds[K_IDX][1]))
    ])
    max_idx = np.argmax(u_space)
    min_idx = np.argmin(u_space)
    plt.plot(
        [b_space[min_idx], b_space[max_idx]],
        [u_space[min_idx], u_space[max_idx]],
        'x', color='purple'
    )

    b_space = np.linspace(b_space[min_idx], b_space[max_idx], 1000)

    plt.plot(
        b_space,
        UB_TRAJECTORY(b_space, lcds[K_IDX][0]),
        color='purple', linestyle='--'
    )
    plt.show()

    for contour in w0_contours:
        plt.plot([x.x for x in contour], [x.y for x in contour])

    plt.plot(
        lcds[K_IDX][1], lcds[K_IDX][2],
        color='purple', linestyle='--'
    )
    mask = ~np.isnan(grid.magnetic_amp_grids[K_IDX].discrete)
    plt.plot(
        grid.x_grid[mask],
        grid.y_grid[mask],
        '^', alpha=0.01, color='grey'
    )
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", nargs="?", choices=["data", "plot"])
    parser.add_argument("-f", "--filepath", default="grid_data.npz", 
                        help="Path to save/load grid data (default: grid_data.npz)")
    args = parser.parse_args()

    if args.mode is None:
        gen_data()
        save_data(args.filepath)
        plot()
    elif args.mode == "data":
        gen_data()
        save_data(args.filepath)
    elif args.mode == "plot":
        load_data(args.filepath)
        plot()

