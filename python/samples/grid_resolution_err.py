import UBKLib
import UBKLib.utils
import UBKLib.analytical_case
import numpy as np
import matplotlib.pyplot as plt

RESOLUTIONS = [200, 500, 1000, 1500, 2000]
COLORS = [
    "#e41a1c",  # red
    "#377eb8",  # blue
    "#4daf4a",  # green
    "#984ea3",  # purple
    "#ff7f00",  # orange
]

ELECTRON_REST_MASS = 510.99895069
PROTON_REST_MASS = 939272.0894329

MU_VALS = np.linspace(0, 2000, 200)

R_MAX = 15


def proton_hamiltonian(b, u, mu):
    return UBKLib.relatavistic_hamiltonian(b, u, mu, 1, PROTON_REST_MASS)


def proton_ub_trajectory(b, mu, h):
    return UBKLib.relativistic_ub_trajectory(b, mu, 1, PROTON_REST_MASS, h)


def electron_hamiltonian(b, u, mu):
    return UBKLib.relatavistic_hamiltonian(b, u, mu, -1, ELECTRON_REST_MASS)


def electron_ub_trajectory(b, mu, h):
    return UBKLib.relativistic_ub_trajectory(b, mu, -1, ELECTRON_REST_MASS, -h)


def equatorial_dipole(x, y):
    r = np.hypot(x, y)
    if r > R_MAX or r < 1.5:
        return [np.nan]
    return [UBKLib.equatorial_dipole_amplitude(x, y)]


analytical = UBKLib.analytical_case.VollandSternAndDipole()

if __name__ == "__main__":
    grids = [UBKLib.FullGrid(
        [0],
        UBKLib.volland_stern_potential,
        equatorial_dipole,
        resolution=res) for res in RESOLUTIONS
    ]
    analytical_trajectory = [
        UBKLib.continuous_lcds_ub_trajectory(
            analytical.lower_contour,
            analytical.upper_contour,
            lambda b: proton_ub_trajectory(b, mu, 0)
        ) for mu in MU_VALS
    ]
    analytical_lower_radii = np.array([analytical.radius(traj.lower_intercept.B) for traj in analytical_trajectory])
    analytical_upper_radii = np.array([analytical.radius(traj.upper_intercept.B) for traj in analytical_trajectory])
    
    figUB, axs = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    axs[0, 1].sharey(axs[0, 0])
    axs[1, 1].sharey(axs[1, 0])
    axs[0, 1].tick_params(labelleft=False)
    axs[1, 1].tick_params(labelleft=False)

    for ele in axs:
        for x in ele:
            x.margins(x=0.1, y=0.1)
            x.axhline(y=0, color='black', linestyle='--')

    axLowerDelX = axs[0][0]
    axUpperDelX = axs[0][1]

    axLowerDelX.set_title("Lower contour")
    axUpperDelX.set_title("Upper contour")

    axLowerDelTheta = axLowerDelX.twinx()
    axUpperDelTheta = axUpperDelX.twinx()
    axLowerDelTheta.sharey(axUpperDelTheta)
    axLowerDelTheta.tick_params(labelright=False)   

    axLowerDelX.set_ylabel(r"$\Delta$x $(R_E)$")
    axUpperDelTheta.set_ylabel(r"$\Delta\theta$ (radians)")

    axLowerDelR = axs[1][0]
    axUpperDelR = axs[1][1]
    axLowerDelR.set_xlabel(r"$\mu$ (kV/nT)")
    axUpperDelR.set_xlabel(r"$\mu$ (kV/nT)")
    axLowerDelR.set_ylabel(r"$\Delta$R absolute $(R_E)$")
    axLowerDelR_rel = axLowerDelR.twinx()
    axUpperDelR_rel = axUpperDelR.twinx()
    axLowerDelR_rel.tick_params(labelright=False)
    axUpperDelR_rel.set_ylabel(r"$\Delta$R relative (%)")

    for j, grid in enumerate(grids):
        grid.calc_potential_grid()
        grid.calc_magnetic_amp_grids()
        grid.calc_cross_products_grids()
        w0_contours = grid.find_cross_product_zeros()[0]
        lower_contour_pos = []
        upper_contour_pos = []
        for mu in MU_VALS:

            grid.calc_hamiltonian_grids(lambda b, u: proton_hamiltonian(b, u, mu))

            lcds = grid.find_lcds_contours(50)[0]
            
            min_idx = np.argmin(lcds[2])
            max_idx = np.argmax(lcds[2])
            lcds_max = [lcds[1][min_idx], lcds[2][min_idx]]
            lcds_min = [lcds[1][max_idx], lcds[2][max_idx]]

            lower_contour_pos.append(lcds_min)
            upper_contour_pos.append(lcds_max)
        
        axLowerDelX.plot(MU_VALS, [x[0] for x in lower_contour_pos], color=COLORS[j], label=f"{RESOLUTIONS[j]}")
        axUpperDelX.plot(MU_VALS, [x[0] for x in upper_contour_pos], color=COLORS[j])

        axLowerDelTheta.plot(
            MU_VALS,
            [np.arcsin(x[0] / x[1]) for x in lower_contour_pos],
            color=COLORS[j], alpha=0.75, linestyle="--"
        )
        axUpperDelTheta.plot(
            MU_VALS,
            [np.arcsin(x[0] / x[1]) for x in upper_contour_pos],
            color=COLORS[j], alpha=0.75, linestyle='--'
        )
        
        r_lower = np.array([np.hypot(x[0], x[1]) for x in lower_contour_pos])
        axLowerDelR.plot(MU_VALS, r_lower - analytical_lower_radii, color=COLORS[j])
        axLowerDelR_rel.plot(
            MU_VALS,
            (r_lower - analytical_lower_radii) / analytical_lower_radii * 100,
            color=COLORS[j], alpha=0.75, linestyle="--"
        )

        r_upper = np.array([np.hypot(x[0], x[1]) for x in upper_contour_pos])
        axUpperDelR.plot(MU_VALS, r_upper - analytical_upper_radii, color=COLORS[j])
        axUpperDelR_rel.plot(
            MU_VALS,
            (r_upper - analytical_upper_radii) / analytical_upper_radii * 100,
            color=COLORS[j], alpha=0.75, linestyle="--"
        )

    figUB.legend(title="grid resolution")
    figUB.suptitle("Effect of grid resolution on proton LCDS trajectory errors")
    plt.tight_layout(rect=[0, 0, 0.92, 1])  
    plt.show()
