import UBKLibpp
import UBKLib
import matplotlib.pyplot as plt
import argparse
import multiprocessing as mp

K_VALUES = [1, 10, 100, 1000]
RESOLUTION = 1000
THREAD_COUNT = mp.cpu_count()

grid = None


def POTENTIAL(x, y):
    return UBKLib.volland_stern_potential(x, y, kp=0)


def gen_data():
    global grid
    grid = UBKLib.Grid(
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
    grid = UBKLib.Grid.load(filepath)


def plot():
    plt.plot(grid.magnetic_amp_grids[0].flatten(), grid.potential_grid.flatten(), '^', alpha=0.1)
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

