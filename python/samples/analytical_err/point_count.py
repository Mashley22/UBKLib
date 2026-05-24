import UBKLib
import UBKLib.utils
import UBKLib.analytical_case
import numpy as np
import matplotlib.pyplot as plt

point_counts = [50_000, 100_000, 500_000, 1_000_000]
colors = [
    "#e41a1c",  # red
    "#377eb8",  # blue
    "#4daf4a",  # green
    "#984ea3",  # purple
    "#ff7f00",  # orange
]

ELECTRON_REST_MASS = 510.99895069
PROTON_REST_MASS = 939272.0894329

PERCENTILE = 100

INNER_LIM = 1.25
OUTER_LIM = 15

point_clouds = []

lower_contours = []
upper_contours = []

analytical = UBKLib.analytical_case.VollandSternAndDipole()

for count in point_counts:
    temp = {}
    temp['x'], temp['y'] = UBKLib.utils.random_xy(r_max=16, n=count)
    point_clouds.append(temp)


for points in point_clouds:
    points['B'] = UBKLib.equatorial_dipole_amplitude(points['x'], points['y'])
    points['U'] = UBKLib.volland_stern_potential(points['x'], points['y'])

    bins = UBKLib.utils.generate_inv_cuberoot_bins(
        np.min(points['B']),
        np.max(points['B']),
        int(np.cbrt(len(points['B'])))
    )
    bins2 = np.linspace(
        np.min(points['B']),
        np.max(points['B']),
        num=int(np.cbrt(len(points['B']))),
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
    
    temp = {}

    temp['UB'] = UBKLib.generate_ub_spline([x.pos for x in upper_points])
    temp['x'], temp['y'] = UBKLib.generate_realSpace_splines([x.pos for x in upper_points])

    upper_contours.append(temp)

    temp = {}

    temp['UB'] = UBKLib.generate_ub_spline([x.pos for x in lower_points])
    temp['x'], temp['y'] = UBKLib.generate_realSpace_splines([x.pos for x in lower_points])

    lower_contours.append(temp)


lin = np.linspace(UBKLib.field_models.dipole.SURFACE_STRENGTH / (15 ** 3), UBKLib.field_models.dipole.SURFACE_STRENGTH / (1.05 ** 3), 1000)

fig, axs = plt.subplots(2, 2, figsize=(12, 8), sharex=True)

axLower = axs[0, 0]
axUpper = axs[0, 1]

for i in range(len(lower_contours)):
    axLower.plot(lin, lower_contours[i]['UB'](lin) - analytical.lower_contour(lin), alpha=0.7, color=colors[i], label=point_counts[i])
axLower.axhline(y=0, color='black', linestyle='--')
axLower.set_title(r'Lower contour $\Delta$U(B)')
axLower.set_ylabel(r'$\Delta$ U (kV)')

for i in range(len(lower_contours)):
    axUpper.plot(lin, upper_contours[i]['UB'](lin) - analytical.upper_contour(lin), alpha=0.7, color=colors[i])
axUpper.axhline(y=0, color='black', linestyle='--')
axUpper.set_title(r'Upper contour $\Delta$U(B)')

axLower = axs[1, 0]
axUpper = axs[1, 1]
axs[1, 0].sharey(axs[1, 1])


for i in range(len(lower_contours)):
    axUpper.plot(lin, np.hypot(upper_contours[i]['x'](lin), upper_contours[i]['y'](lin)) - analytical.radius(lin), alpha=0.7, color=colors[i])
axUpper.axhline(y=0, color='black', linestyle='--')
axUpper.set_xlabel('B (nT)')
axUpper.set_title(r'Upper contour $\Delta$R(B)')

for i in range(len(lower_contours)):
    axLower.plot(lin, np.hypot(lower_contours[i]['x'](lin), lower_contours[i]['y'](lin)) - analytical.radius(lin), alpha=0.7, color=colors[i])
axLower.axhline(y=0, color='black', linestyle='--')
axLower.set_xlabel('B (nT)')
axLower.set_title(r'Lower contour $\Delta$R(B)')
axLower.set_ylabel(r'$\Delta$R (Re)')
fig.legend(title="point count")
plt.show()

mu_vals = np.linspace(0, 30, 100)

analytical_trajectories = []

for mu in mu_vals:
    analytical_trajectories.append(
        UBKLib.continuous_lcds_ub_trajectory(
            analytical.lower_contour,
            analytical.upper_contour,
            lambda x: UBKLib.relativistic_ub_trajectory(x, mu, -1, ELECTRON_REST_MASS),
        )
    )

numerical_trajectories = []

for i in range(len(point_counts)):
    temp_trajectories = []
    for mu in mu_vals:
        temp_trajectories.append(
            UBKLib.continuous_lcds_ub_trajectory(
                analytical.lower_contour,
                analytical.upper_contour,
                lambda x: UBKLib.relativistic_ub_trajectory(x, mu, -1, ELECTRON_REST_MASS),
            )
        )

    numerical_trajectories.append(temp_trajectories)


fig, axs = plt.subplots(2, 2, figsize=(12, 8), sharex=True, gridspec_kw={"wspace": 0.6})
axLower = axs[0, 0]
axUpper = axs[0, 1]

analytical_radii = np.array([analytical.radius(traj.lower_intercept.B) for traj in analytical_trajectories])
numerical_radii = []

for i in range(len(point_counts)):
    numerical_radii.append(np.array(
        [
            np.hypot(
                lower_contours[i]['x'](traj.lower_intercept.B),
                lower_contours[i]['y'](traj.lower_intercept.B)
            ) for traj in numerical_trajectories[i]
        ]
    ))

for i in range(len(point_counts)):
    absolute_error = numerical_radii[i] - analytical_radii
    axLower.plot(mu_vals, absolute_error, color=colors[i], 
                 alpha=0.8, linewidth=2, label=f'{point_counts[i]}')

axLower.axhline(y=0, color='black', linestyle='--', linewidth=1)
axLower.set_xlabel(r'$\mu$ (KeV/nT)', fontsize=12)
axLower.set_ylabel(r'Absolute Error $\Delta R$ (Re)', fontsize=12)
axLower.tick_params(axis='y')
axLower.grid(True, alpha=0.3)

ax2 = axLower.twinx()

legend_handles = []

for i in range(len(point_counts)):
    relative_error = 100 * (numerical_radii[i] - analytical_radii) / analytical_radii
    ax2.plot(mu_vals, relative_error, color=colors[i], 
             alpha=0.4, linewidth=1.5, linestyle='--')

ax2.axhline(y=0, color='black', linestyle=':', linewidth=0.8)
ax2.set_ylabel(r'Relative Error (%)', fontsize=12)
ax2.tick_params(axis='y')

axLower.set_title("Electron error at lower W = 0 contour")

analytical_radii = np.array([analytical.radius(traj.upper_intercept.B) for traj in analytical_trajectories])
numerical_radii = []

for i in range(len(point_counts)):
    numerical_radii.append(np.array(
        [
            np.hypot(
                upper_contours[i]['x'](traj.upper_intercept.B),
                upper_contours[i]['y'](traj.upper_intercept.B)
            ) for traj in numerical_trajectories[i]
        ]
    ))


for i in range(len(point_counts)):
    absolute_error = numerical_radii[i] - analytical_radii
    axUpper.plot(mu_vals, absolute_error, color=colors[i], 
                 alpha=0.8, linewidth=2)

axUpper.axhline(y=0, color='black', linestyle='--', linewidth=1)
axUpper.set_xlabel(r'$\mu$ (KeV/nT)', fontsize=12)
axUpper.set_ylabel(r'Absolute Error $\Delta R$ (Re)', fontsize=12)
axUpper.tick_params(axis='y')
axUpper.grid(True, alpha=0.3)

ax3 = axUpper.twinx()

legend_handles = []

for i in range(len(point_counts)):
    relative_error = 100 * (numerical_radii[i] - analytical_radii) / analytical_radii
    ax3.plot(mu_vals, relative_error, color=colors[i], 
             alpha=0.4, linewidth=1.5, linestyle='--')

ax3.axhline(y=0, color='black', linestyle=':', linewidth=0.8)
ax3.set_ylabel(r'Relative Error (%)', fontsize=12)
ax3.tick_params(axis='y')


axUpper.set_title("Electron error at upper W = 0 contour")

analytical_trajectories = []

axLower = axs[1, 0]
axUpper = axs[1, 1]

for mu in mu_vals:
    analytical_trajectories.append(
        UBKLib.continuous_lcds_ub_trajectory(
            analytical.lower_contour,
            analytical.upper_contour,
            lambda x: UBKLib.relativistic_ub_trajectory(x, mu, 1, PROTON_REST_MASS),
        )
    )

numerical_trajectories = []

for i in range(len(point_counts)):
    temp_trajectories = []
    for mu in mu_vals:
        temp_trajectories.append(
            UBKLib.continuous_lcds_ub_trajectory(
                analytical.lower_contour,
                analytical.upper_contour,
                lambda x: UBKLib.relativistic_ub_trajectory(x, mu, 1, PROTON_REST_MASS),
            )
        )

    numerical_trajectories.append(temp_trajectories)

analytical_radii = np.array([analytical.radius(traj.lower_intercept.B) for traj in analytical_trajectories])
numerical_radii = []

for i in range(len(point_counts)):
    numerical_radii.append(np.array(
        [
            np.hypot(
                lower_contours[i]['x'](traj.lower_intercept.B),
                lower_contours[i]['y'](traj.lower_intercept.B)
            ) for traj in numerical_trajectories[i]
        ]
    ))

for i in range(len(point_counts)):
    absolute_error = numerical_radii[i] - analytical_radii
    axLower.plot(mu_vals, absolute_error, color=colors[i], 
                 alpha=0.8, linewidth=2)

axLower.axhline(y=0, color='black', linestyle='--', linewidth=1)
axLower.set_xlabel(r'$\mu$ (KeV/nT)', fontsize=12)
axLower.set_ylabel(r'Absolute Error $\Delta R$ (Re)', fontsize=12)
axLower.tick_params(axis='y')
axLower.grid(True, alpha=0.3)

ax2 = axLower.twinx()

legend_handles = []

for i in range(len(point_counts)):
    relative_error = 100 * (numerical_radii[i] - analytical_radii) / analytical_radii
    ax2.plot(mu_vals, relative_error, color=colors[i], 
             alpha=0.4, linewidth=1.5, linestyle='--')

ax2.axhline(y=0, color='black', linestyle=':', linewidth=0.8)
ax2.set_ylabel(r'Relative Error (%)', fontsize=12)
ax2.tick_params(axis='y')

axLower.set_title("Lower contour trajectory error")

analytical_radii = np.array([analytical.radius(traj.upper_intercept.B) for traj in analytical_trajectories])
numerical_radii = []

for i in range(len(point_counts)):
    numerical_radii.append(np.array(
        [
            np.hypot(
                upper_contours[i]['x'](traj.upper_intercept.B),
                upper_contours[i]['y'](traj.upper_intercept.B)
            ) for traj in numerical_trajectories[i]
        ]
    ))


for i in range(len(point_counts)):
    absolute_error = numerical_radii[i] - analytical_radii
    axUpper.plot(mu_vals, absolute_error, color=colors[i], 
                 alpha=0.8, linewidth=2)

axUpper.axhline(y=0, color='black', linestyle='--', linewidth=1)
axUpper.set_xlabel(r'$\mu$ (KeV/nT)', fontsize=12)
axUpper.set_ylabel(r'Absolute Error $\Delta R$ (Re)', fontsize=12)
axUpper.tick_params(axis='y')
axUpper.grid(True, alpha=0.3)

ax3 = axUpper.twinx()


for i in range(len(point_counts)):
    relative_error = 100 * (numerical_radii[i] - analytical_radii) / analytical_radii
    ax3.plot(mu_vals, relative_error, color=colors[i], 
             alpha=0.4, linewidth=1.5, linestyle='--')

ax3.axhline(y=0, color='black', linestyle=':', linewidth=0.8)
ax3.set_ylabel(r'Relative Error (%)', fontsize=12)
ax3.tick_params(axis='y')


axUpper.set_title("Upper contour trajectory error")
fig.legend(title="point count")
plt.show()
