import UBKLib
import UBKLib.utils
import UBKLib.analytical_case
import numpy as np
import matplotlib.pyplot as plt
import argparse

POINT_COUNT = 100_000

ELECTRON_REST_MASS = 510.99895069
PROTON_REST_MASS = 939272.0894329

PERCENTILE = 100

RUN_COUNT = 1000

MU_VALS = np.linspace(0, 3, 100)

data_loc = 'run_to_run.npy'

upper_col = 'blue'
lower_col = 'red'

analytical = UBKLib.analytical_case.VollandSternAndDipole()

analytical_trajectories = [
    UBKLib.continuous_lcds_ub_trajectory(
        analytical.lower_contour,
        analytical.upper_contour,
        lambda x: UBKLib.relativistic_ub_trajectory(x, mu, -1, ELECTRON_REST_MASS),
    ) for mu in MU_VALS
]

electron_analytical_lower = np.array([analytical.radius(traj.lower_intercept.B) for traj in analytical_trajectories])
electron_analytical_upper = np.array([analytical.radius(traj.upper_intercept.B) for traj in analytical_trajectories])

analytical_trajectories = [
    UBKLib.continuous_lcds_ub_trajectory(
        analytical.lower_contour,
        analytical.upper_contour,
        lambda x: UBKLib.relativistic_ub_trajectory(x, mu, 1, PROTON_REST_MASS),
    ) for mu in MU_VALS
]

proton_analytical_lower = np.array([analytical.radius(traj.lower_intercept.B) for traj in analytical_trajectories])
proton_analytical_upper = np.array([analytical.radius(traj.upper_intercept.B) for traj in analytical_trajectories])

proton_upper = [[] for _ in range(len(MU_VALS))]
proton_lower = [[] for _ in range(len(MU_VALS))]
electron_upper = [[] for _ in range(len(MU_VALS))]
electron_lower = [[] for _ in range(len(MU_VALS))]

proton_upper_stddev = []
proton_lower_stddev = []
electron_upper_stddev = []
electron_lower_stddev = []

proton_upper_avg = []
proton_lower_avg = []
electron_upper_avg = []
electron_lower_avg = []


def save_data():
    global proton_lower_stddev
    global proton_upper_stddev
    global proton_lower_avg
    global proton_upper_avg

    global electron_lower_stddev
    global electron_upper_stddev
    global electron_lower_avg
    global electron_upper_avg

    data = {
        "proton_lower_stddev": proton_lower_stddev,
        "proton_upper_stddev": proton_upper_stddev,
        "proton_lower_avg": proton_lower_avg,
        "proton_upper_avg": proton_upper_avg,
        "electron_lower_stddev": electron_lower_stddev,
        "electron_upper_stddev": electron_upper_stddev,
        "electron_lower_avg": electron_lower_avg,
        "electron_upper_avg": electron_upper_avg,
    }

    np.save(data_loc, data, allow_pickle=True)


def load_data():
    global proton_lower_stddev
    global proton_upper_stddev
    global proton_lower_avg
    global proton_upper_avg

    global electron_lower_stddev
    global electron_upper_stddev
    global electron_lower_avg
    global electron_upper_avg

    data = np.load(data_loc, allow_pickle=True).item()

    proton_lower_stddev = data["proton_lower_stddev"]
    proton_upper_stddev = data["proton_upper_stddev"]
    proton_lower_avg = data["proton_lower_avg"]
    proton_upper_avg = data["proton_upper_avg"]

    electron_lower_stddev = data["electron_lower_stddev"]
    electron_upper_stddev = data["electron_upper_stddev"]
    electron_lower_avg = data["electron_lower_avg"]
    electron_upper_avg = data["electron_upper_avg"]


def gen_data():
    
    global proton_lower_stddev
    global proton_upper_stddev
    global proton_lower_avg
    global proton_upper_avg

    global electron_lower_stddev
    global electron_upper_stddev
    global electron_lower_avg
    global electron_upper_avg

    for i in range(RUN_COUNT):
        points = {}
        points['x'], points['y'] = UBKLib.utils.random_xy(r_max=16, n=POINT_COUNT)
        points['B'] = UBKLib.equatorial_dipole_amplitude(points['x'], points['y'])
        points['U'] = UBKLib.volland_stern_potential(points['x'], points['y'])
        
        bins = UBKLib.utils.generate_inv_cuberoot_bins(
            np.min(points['B']),
            np.max(points['B']),
            int(np.cbrt(POINT_COUNT))
        )
        bins2 = np.linspace(
            np.min(points['B']),
            np.max(points['B']),
            num=int(np.cbrt(POINT_COUNT)),
            endpoint=True
        )
        bins = np.concatenate([bins, bins2])
        bins = np.sort(bins)[::2]

        upper_points, lower_points = UBKLib.w0_points_from_cloud(
            points['x'],
            points['y'],
            points['B'],
            points['U'],
            bins,
            PERCENTILE
        )
           
        upper = UBKLib.generate_ub_spline([x.pos for x in upper_points])
        upper_x, upper_y = UBKLib.generate_realSpace_splines([x.pos for x in upper_points])
        
        lower = UBKLib.generate_ub_spline([x.pos for x in lower_points])
        lower_x, lower_y = UBKLib.generate_realSpace_splines([x.pos for x in lower_points])

        trajectories = [
            UBKLib.continuous_lcds_ub_trajectory(
                lower,
                upper,
                lambda x: UBKLib.relativistic_ub_trajectory(x, mu, -1, ELECTRON_REST_MASS),
            ) for mu in MU_VALS
        ]
        
        for i, traj in enumerate(trajectories):
            electron_upper[i].append(np.hypot(upper_x(traj.upper_intercept.B), upper_y(traj.upper_intercept.B)))
            electron_lower[i].append(np.hypot(lower_x(traj.lower_intercept.B), lower_y(traj.lower_intercept.B)))

        trajectories = [
            UBKLib.continuous_lcds_ub_trajectory(
                lower,
                upper,
                lambda x: UBKLib.relativistic_ub_trajectory(x, mu, 1, PROTON_REST_MASS),
            ) for mu in MU_VALS
        ]

        for i, traj in enumerate(trajectories):
            proton_upper[i].append(np.hypot(upper_x(traj.upper_intercept.B), upper_y(traj.upper_intercept.B)))
            proton_lower[i].append(np.hypot(lower_x(traj.lower_intercept.B), lower_y(traj.lower_intercept.B)))

    electron_lower_stddev = np.array([np.std(ele) for ele in electron_lower])
    electron_upper_stddev = np.array([np.std(ele) for ele in electron_upper])

    proton_lower_stddev = np.array([np.std(ele) for ele in proton_lower])
    proton_upper_stddev = np.array([np.std(ele) for ele in proton_upper])

    electron_lower_avg = np.array([np.average(ele) for ele in electron_lower])
    electron_upper_avg = np.array([np.average(ele) for ele in electron_upper])

    proton_lower_avg = np.array([np.average(ele) for ele in proton_lower])
    proton_upper_avg = np.array([np.average(ele) for ele in proton_upper])


def plot():
    fig, (axLower, axUpper) = plt.subplots(1, 2, figsize=(12, 4), sharey=True, gridspec_kw={'wspace': 0.1})
    axLower.set_xlabel(r"$\mu (kV/nT)$")
    axLower.set_ylabel(r"$\Delta R (Re)$")
    axLower.set_title("Electron error at the W = 0 contours")
    val = electron_upper_avg - electron_analytical_upper
    axLower.plot(MU_VALS, val + electron_upper_stddev, color=upper_col, label='upper contour')
    axLower.plot(MU_VALS, val - electron_upper_stddev, color=upper_col)
    axLower.fill_between(MU_VALS, val + electron_upper_stddev, val - electron_upper_stddev, alpha=0.3, color=upper_col)

    val = electron_lower_avg - electron_analytical_lower
    axLower.plot(MU_VALS, val + electron_lower_stddev, color=lower_col, label='lower contour')
    axLower.plot(MU_VALS, val - electron_lower_stddev, color=lower_col)
    axLower.fill_between(MU_VALS, val + electron_lower_stddev, val - electron_lower_stddev, alpha=0.3, color=lower_col)

    axUpper.set_xlabel(r"$\mu (kV/nT)$")
    axUpper.set_title("Proton error at the W = 0 contours")
    axUpper.legend()
    val = proton_upper_avg - proton_analytical_upper
    axUpper.plot(MU_VALS, val + proton_upper_stddev, color=upper_col)
    axUpper.plot(MU_VALS, val - proton_upper_stddev, color=upper_col)
    axUpper.fill_between(MU_VALS, val + proton_upper_stddev, val - proton_upper_stddev, alpha=0.3, color=upper_col)

    val = proton_lower_avg - proton_analytical_lower
    axUpper.plot(MU_VALS, val + proton_lower_stddev, color=lower_col)
    axUpper.plot(MU_VALS, val - proton_lower_stddev, color=lower_col)
    axUpper.fill_between(MU_VALS, val + proton_lower_stddev, val - proton_lower_stddev, alpha=0.3, color=lower_col)
    fig.legend()
    plt.tight_layout(rect=[0, 0, 1, 0.88])
    axLower.margins(x=0, y=0)
    axUpper.margins(x=0, y=0)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", nargs="?", choices=["data", "plot"])
    args = parser.parse_args()
    if args.mode is None:
        gen_data()
        save_data()
        plot()
    elif args.mode == "data":
        gen_data()
        save_data()
    elif args.mode == "plot":
        load_data()
        plot()
