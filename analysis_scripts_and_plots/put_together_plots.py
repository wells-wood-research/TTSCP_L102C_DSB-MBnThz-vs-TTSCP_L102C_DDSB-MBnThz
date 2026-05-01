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
# FIGURE SETUP (5 rows × 2 cols)
# -------------------------
fig, axes = plt.subplots(5, 2, figsize=(16, 24), sharex=False)

# -------------------------
# LOOP SYSTEMS
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
        mbnthz_idx = np.arange(1543, 1589)

        ca_traj = traj.atom_slice(ca_idx)

        # Align CA atoms
        ca_traj.superpose(ca_traj, 0)

        # -------------------------
        # 1. RMSD (CA)
        # -------------------------
        rmsd = md.rmsd(ca_traj, ca_traj, 0) * 10
        axes[0, col].plot(time, rmsd, label=f"rep{i+1}")

        # -------------------------
        # 2. Rg (CA)
        # -------------------------
        rg = md.compute_rg(ca_traj) * 10
        axes[1, col].plot(time, rg, label=f"rep{i+1}")

        # -------------------------
        # 3. RMSF (CA)
        # -------------------------
        rmsf = md.rmsf(ca_traj, ca_traj, 0) * 10
        axes[2, col].plot(rmsf, label=f"rep{i+1}")

        # -------------------------
        # 4. SASA (protein backbone ONLY — FIXED)
        # -------------------------
        backbone_traj = traj.atom_slice(backbone_idx)
        sasa_backbone = md.shrake_rupley(backbone_traj, probe_radius=0.15)
        sasa_backbone = sasa_backbone.sum(axis=1) * 100  # convert nm² → Å²
        axes[3, col].plot(time, sasa_backbone, label=f"rep{i+1}")

        # -------------------------
        # 5. SASA (MBnThz ONLY — FIXED)
        # -------------------------
        mbnthz_traj = traj.atom_slice(mbnthz_idx)
        sasa_mbnthz = md.shrake_rupley(mbnthz_traj, probe_radius=0.15)
        sasa_mbnthz = sasa_mbnthz.sum(axis=1) * 100  # convert nm² → Å²
        axes[4, col].plot(time, sasa_mbnthz, label=f"rep{i+1}")

    # -------------------------
    # LABELS PER COLUMN
    # -------------------------
    titles = ["RMSD (Å)", "Rg (Å)", "RMSF (Å)", "SASA Backbone (Å²)", "SASA MBnThz (Å²)"]

    for row in range(5):
        axes[row, col].set_title(f"{titles[row]} - {system}", fontsize=16)
        axes[row, col].legend()

        if row != 2:
            axes[row, col].set_xlabel("Time (ns)")
        else:
            axes[row, col].set_xlabel("Residue Index")

# Shared Y labels
axes[0, 0].set_ylabel("RMSD (Å)")
axes[1, 0].set_ylabel("Rg (Å)")
axes[2, 0].set_ylabel("RMSF (Å)")
axes[3, 0].set_ylabel("SASA (Å²)")
axes[4, 0].set_ylabel("SASA (Å²)")

plt.tight_layout()
plt.savefig("all_plots.png", dpi=300)
plt.show()
