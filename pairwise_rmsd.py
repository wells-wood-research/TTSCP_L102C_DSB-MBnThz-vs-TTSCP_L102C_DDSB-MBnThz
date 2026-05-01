import MDAnalysis as mda
from MDAnalysis.analysis import align, diffusionmap
import glob
import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# SETTINGS
# -------------------------
systems = {
    "DDSB": {
        "traj": "TTSCP_L102C_DDSB_simoutput/*.dcd",
        "top": "TTSCP_L102C_DDSB_simoutput/TTSCP_L102C_DDSB_energyminimized_solvated.pdb"
    },
    "DSB": {
        "traj": "TTSCP_L102C_DSB_simoutput/*.dcd",
        "top": "TTSCP_L102C_DSB_simoutput/TTSCP_L102C_DSB_energyminimized_solvated.pdb"
    }
}

SIM_TIME = 200  # ns

# -------------------------
# STYLE
# -------------------------
LABEL = 24
TICKS = 18

fig, axes = plt.subplots(5, 2, figsize=(14, 24))

# -------------------------
# LOOP SYSTEMS
# -------------------------
for col, (system, info) in enumerate(systems.items()):
    files = sorted(glob.glob(info["traj"]))

    print(f"Processing {system}...")

    for row in range(5):

        if row >= len(files):
            axes[row, col].axis("off")
            continue

        u = mda.Universe(info["top"], files[row])

        align.AlignTraj(u, u, select="name CA", in_memory=True).run()

        dm = diffusionmap.DistanceMatrix(u, select="name CA").run()
        mat = dm.results.dist_matrix

        n = mat.shape[0]
        time = np.linspace(0, SIM_TIME, n)

        # -------------------------
        # TIME-RESCALED HEATMAP
        # -------------------------
        im = axes[row, col].imshow(
            mat,
            cmap="viridis",
            vmin=0,
            vmax=5,
            extent=[0, SIM_TIME, SIM_TIME, 0],  # 🔥 TIME-BASED AXES
            aspect="auto"
        )

        axes[row, col].tick_params(axis='both', labelsize=TICKS)

# -------------------------
# GLOBAL LABELS
# -------------------------
fig.supxlabel("Time (ns)", fontsize=LABEL)
fig.supylabel("Time (ns)", fontsize=LABEL)

# -------------------------
# COLORBAR
# -------------------------
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label("RMSD (Å)", fontsize=LABEL)
cbar.ax.tick_params(labelsize=TICKS)

# -------------------------
# LAYOUT
# -------------------------
plt.subplots_adjust(
    wspace=0.25,
    hspace=0.25,
    right=0.9
)

plt.savefig("pairwise_rmsd_time_based.png", dpi=300, bbox_inches="tight")
plt.show()
