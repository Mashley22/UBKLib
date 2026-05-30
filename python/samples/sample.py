import matplotlib.pyplot as plt
import argparse
import multiprocessing as mp
import numpy as np
from matplotlib.patches import Wedge

import UBKLibpp
import UBKLib
import UBKLib.types

K_VALUES = [1, 10, 100, 1000]
RESOLUTION = 1000
THREAD_COUNT = mp.cpu_count()

ELECTRON_REST_MASS = 510.99895069
PROTON_REST_MASS = 939272.0894329

MU = 0.224
CHARGE = -1
K_IDX = 0

grid = None

LOWER_CONTOUR_COLOUR = "blue"
UPPER_CONTOUR_COLOUR = "red"

COLOURS = [
    "#4daf4a",  # green
    "#984ea3",  # purple
    "#ff7f00",  # orange
    "#ffd700",  # gold / yellow
    "#a65628",  # brown
    "#999999",  # gray
    "#6b8e23",  # olive drab
    "#daa520",  # goldenrod
    "#8b4513",  # saddle brown
    "#9acd32",  # yellow green
]


def HAMILTONIAN(b, u):
    return UBKLib.relatavistic_hamiltonian(b, u, MU, CHARGE, ELECTRON_REST_MASS)


def POTENTIAL(x, y):
    return UBKLib.volland_stern_potential(x, y)


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
    lcds = grid.find_lcds_contours(50)
    w0_contours = grid.find_cross_product_zeros()[K_IDX]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    mask = grid.x_grid.flatten() > 0
    ax2.plot(
        grid.magnetic_amp_grids[K_IDX].discrete.flatten()[mask],
        grid.potential_grid.discrete.flatten()[mask],
        'x', alpha=0.01, color='grey'
    )
    ax2.set_xlabel('B (nT)')
    ax2.set_ylabel('U (kV)')
    idx = 0
    for i in range(len(w0_contours)):
        contour = w0_contours[i]
        tpe = UBKLib.w0_contour_type(contour)
        if tpe == UBKLib.types.W0ContourType.NONE:
            color = COLOURS[idx]
            idx += 1
        if tpe == UBKLib.types.W0ContourType.LOWER:
            color = LOWER_CONTOUR_COLOUR
        if tpe == UBKLib.types.W0ContourType.UPPER:
            color = UPPER_CONTOUR_COLOUR
        ax2.plot([x.B for x in contour], [x.U for x in contour], color=color)

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
    ax2.plot(
        [b_space[min_idx], b_space[max_idx]],
        [u_space[min_idx], u_space[max_idx]],
        'x', color='purple'
    )

    b_space = np.linspace(b_space[min_idx], b_space[max_idx], 1000)
    ax2.plot(
        b_space,
        UB_TRAJECTORY(b_space, lcds[K_IDX][0]),
        color='purple', linestyle='--'
    )
    
    idx = 0
    for i in range(len(w0_contours)):
        contour = w0_contours[i]
        tpe = UBKLib.w0_contour_type(contour)
        if tpe == UBKLib.types.W0ContourType.NONE:
            color = COLOURS[idx]
            idx += 1
        if tpe == UBKLib.types.W0ContourType.LOWER:
            color = LOWER_CONTOUR_COLOUR
        if tpe == UBKLib.types.W0ContourType.UPPER:
            color = UPPER_CONTOUR_COLOUR
        contour = w0_contours[i]
        ax1.plot([x.x for x in contour], [x.y for x in contour], color=color)
    ax1.plot(
        [lcds[K_IDX][1][min_idx], lcds[K_IDX][1][max_idx]],
        [lcds[K_IDX][2][min_idx], lcds[K_IDX][2][max_idx]],
        'x', color='purple'
    )

    ax1.plot(
        lcds[K_IDX][1], lcds[K_IDX][2],
        color='purple', linestyle='--'
    )
    mask = ~np.isnan(grid.magnetic_amp_grids[K_IDX].discrete)
    ax1.plot(
        grid.x_grid[mask],
        grid.y_grid[mask],
        'x', alpha=0.01, color='grey'
    )
    dayside = Wedge((0, 0), 1, 270, 90, facecolor='white', edgecolor='black', linewidth=2)
    ax1.add_patch(dayside)

    nightside = Wedge((0, 0), 1, 90, 270, facecolor='black', edgecolor='black', linewidth=2)
    ax1.add_patch(nightside)

    ax1.set_aspect('equal')
    ax1.set_xlabel(r'x ($R_E$)')
    ax1.set_ylabel(r'y ($R_E$)')
    ax1.margins(0)
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
