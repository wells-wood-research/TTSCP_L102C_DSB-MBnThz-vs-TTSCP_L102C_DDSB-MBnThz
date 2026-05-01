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
# STYLE
# -------------------------
TITLE = 28
LABEL = 22
TICKS = 20

colors = {
    "DDSB": "blue",
    "DSB": "orange"
}

# -------------------------
# STORAGE
# -------------------------
data = {
    "DDSB": {"rmsd": [], "rg": [], "rmsf": [], "sasa_bb": [], "sasa_lig": []},
    "DSB": {"rmsd": [], "rg": [], "rmsf": [], "sasa_bb": [], "sasa_lig": []},
}

time_ref = None

# -------------------------
# LOAD DATA
# -------------------------
for system, info in systems.items():
    files = sorted(glob.glob(info["traj"]))

    print(f"Processing {system}...")

    for f in files:
        traj = md.load(f, top=info["top"])
        n_frames = traj.n_frames
        time = np.linspace(0, simulation_ns, n_frames)

        if time_ref is None:
            time_ref = time

        ca_idx = traj.topology.select("name CA")
        backbone_idx = traj.topology.select("backbone")
        mbnthz_idx = np.arange(1538, 1588)

        ca_traj = traj.atom_slice(ca_idx)
        ca_traj.superpose(ca_traj, 0)

        # -------------------------
        # METRICS
        # -------------------------
        data[system]["rmsd"].append(md.rmsd(ca_traj, ca_traj, 0) * 10)
        data[system]["rg"].append(md.compute_rg(ca_traj) * 10)
        data[system]["rmsf"].append(md.rmsf(ca_traj, ca_traj, 0) * 10)

        backbone_traj = traj.atom_slice(backbone_idx)
        sasa_bb = md.shrake_rupley(backbone_traj, probe_radius=0.15).sum(axis=1) * 100
        data[system]["sasa_bb"].append(sasa_bb)

        lig_traj = traj.atom_slice(mbnthz_idx)
        sasa_lig = md.shrake_rupley(lig_traj, probe_radius=0.15).sum(axis=1) * 100
        data[system]["sasa_lig"].append(sasa_lig)

# -------------------------
# HELPERS
# -------------------------
def stats(arr_list):
    arr = np.array(arr_list)
    return arr.mean(axis=0), arr.min(axis=0), arr.max(axis=0)

# -------------------------
# METRIC CONFIG (TITLE UNCHANGED)
# -------------------------
metrics = {
    "rmsd": ("RMSD — mean / min / max", "Å"),
    "rg": ("Rg — mean / min / max", "Å"),
    "rmsf": ("RMSF — mean / min / max", "Å"),
    "sasa_bb": ("SASA (Protein) — mean / min / max", "Å²"),
    "sasa_lig": ("SASA (MBnThz) — mean / min / max", "Å²")
}

# -------------------------
# PLOT FUNCTION
# -------------------------
def plot_metric(metric):

    title, unit = metrics[metric]

    plt.figure(figsize=(8, 6))

    for system in systems:
        mean, minv, maxv = stats(data[system][metric])

        if metric == "rmsf":
            x = np.arange(len(mean))
            plt.plot(x, mean, color=colors[system], linewidth=2)
            plt.fill_between(x, minv, maxv, color=colors[system], alpha=0.25)
        else:
            plt.plot(time_ref, mean, color=colors[system], linewidth=2)
            plt.fill_between(time_ref, minv, maxv, color=colors[system], alpha=0.25)

    # -------------------------
    # TITLE (UNCHANGED)
    # -------------------------
    plt.title(
        title,
        fontsize=TITLE,
        fontweight='bold',
        pad=22
    )

    # -------------------------
    # Y-AXIS LABELS (UNITS ONLY)
    # -------------------------
    plt.xlabel(
        "Residue Index" if metric == "rmsf" else "Time (ns)",
        fontsize=LABEL,
        labelpad=14
    )
    plt.ylabel(unit, fontsize=LABEL, labelpad=14)

    # -------------------------
    # Y-LIMITS
    # -------------------------
    if metric == "rmsd":
        plt.ylim(0, 3.5)
    elif metric == "rg":
        plt.ylim(13.4, 14.4)
    elif metric == "rmsf":
        plt.ylim(0.5, 3.6)
    elif metric == "sasa_bb":
        plt.ylim(6800, 7800)
    elif metric == "sasa_lig":
        plt.ylim(580, 760)

    # -------------------------
    # TICKS
    # -------------------------
    plt.tick_params(axis='both', labelsize=TICKS, length=6, width=1.2)

    # -------------------------
    # LAYOUT FIX
    # -------------------------
    plt.subplots_adjust(top=0.86, left=0.16)

    plt.tight_layout()
    plt.savefig(f"{metric}_ddsb_dsb.png", dpi=300, bbox_inches='tight')
    plt.show()

# -------------------------
# RUN ALL PLOTS
# -------------------------
for m in metrics:
    plot_metric(m)
