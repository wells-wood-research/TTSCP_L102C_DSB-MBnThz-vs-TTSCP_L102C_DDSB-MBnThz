import mdtraj as md
import numpy as np
import glob
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

simulation_ns = 200

# -------------------------
# STYLE SETTINGS
# -------------------------
TITLE = 28
LABEL = 22
TICKS = 20

# -------------------------
# FIGURE SETUP
# -------------------------
fig, axes = plt.subplots(5, 2, figsize=(18, 28))

# -------------------------
# MAIN LOOP
# -------------------------
for col, (system, info) in enumerate(systems.items()):
    files = sorted(glob.glob(info["traj"]))

    if len(files) == 0:
        print(f"No trajectories found for {system}")
        continue

    print(f"Processing {system}...")

    for i, f in enumerate(files):
        traj = md.load(f, top=info["top"])
        n_frames = traj.n_frames
        time = np.linspace(0, simulation_ns, n_frames)

        # -------------------------
        # SELECTIONS
        # -------------------------
        ca_idx = traj.topology.select("name CA")
        backbone_idx = traj.topology.select("backbone")
        mbnthz_idx = np.arange(1538, 1588)

        ca_traj = traj.atom_slice(ca_idx)
        ca_traj.superpose(ca_traj, 0)

        # -------------------------
        # 1. RMSD (CA)
        # -------------------------
        rmsd = md.rmsd(ca_traj, ca_traj, 0) * 10
        axes[0, col].plot(time, rmsd, linewidth=1.5)

        # -------------------------
        # 2. Rg (CA)
        # -------------------------
        rg = md.compute_rg(ca_traj) * 10
        axes[1, col].plot(time, rg, linewidth=1.5)

        # -------------------------
        # 3. RMSF (CA)
        # -------------------------
        rmsf = md.rmsf(ca_traj, ca_traj, 0) * 10
        axes[2, col].plot(rmsf, linewidth=1.5)

        # -------------------------
        # 4. SASA (Backbone ONLY)
        # -------------------------
        backbone_traj = traj.atom_slice(backbone_idx)
        sasa_backbone = md.shrake_rupley(backbone_traj, probe_radius=0.15)
        sasa_backbone = sasa_backbone.sum(axis=1) * 100  # nm² → Å²
        axes[3, col].plot(time, sasa_backbone, linewidth=1.5)

        # -------------------------
        # 5. SASA (MBnThz ONLY)
        # -------------------------
        mbnthz_traj = traj.atom_slice(mbnthz_idx)
        sasa_mbnthz = md.shrake_rupley(mbnthz_traj, probe_radius=0.15)
        sasa_mbnthz = sasa_mbnthz.sum(axis=1) * 100  # nm² → Å²
        axes[4, col].plot(time, sasa_mbnthz, linewidth=1.5)

    # -------------------------
    # TITLES + AXIS FORMATTING
    # -------------------------
    titles = [
        "RMSD",
        "Rg",
        "RMSF",
        "SASA (Backbone)",
        "SASA (MBnThz)"
    ]

    for row in range(5):
        axes[row, col].set_title(
            f"{titles[row]}",
            fontsize=TITLE,
            fontweight='bold'
        )

        if row != 2:
            axes[row, col].set_xlabel("Time (ns)", fontsize=LABEL, labelpad=8)
        else:
            axes[row, col].set_xlabel("Residue", fontsize=LABEL, labelpad=8)

        axes[row, col].tick_params(axis='both', labelsize=TICKS, length=5, width=1)

# -------------------------
# Y LABELS (BOTH COLUMNS)
# -------------------------
ylabel_texts = [
    "Å",
    "Å",
    "Å",
    "Å²",
    "Å²"
]

for row in range(5):
    for col in range(2):
        axes[row, col].set_ylabel(ylabel_texts[row], fontsize=LABEL)

# -------------------------
# Y-LIMITS (YOUR SPECIFIED VALUES)
# -------------------------
for col in range(2):
    axes[0, col].set_ylim(0, 3.5)      # RMSD
    axes[1, col].set_ylim(13.4, 14.4)  # Rg
    axes[2, col].set_ylim(0.5, 3.6)    # RMSF
    axes[3, col].set_ylim(6800, 7800)  # SASA protein
    axes[4, col].set_ylim(580, 760)    # SASA MBnThz

# -------------------------
# LAYOUT FIX
# -------------------------
plt.subplots_adjust(hspace=0.6, wspace=0.4)

plt.savefig("all_plots.png", dpi=300, bbox_inches='tight')
plt.show()
