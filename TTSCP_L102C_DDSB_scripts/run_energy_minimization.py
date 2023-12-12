#!/usr/bin/env pythn
# coding: utf-8

# In[ ]:


from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout


platform=Platform.getPlatformByName("CUDA")

TOTAL_SIM_LENGTH = 20*nanoseconds
TIME_STEP = 2*femtosecond
TOTAL_STEPS = int(TOTAL_SIM_LENGTH/TIME_STEP)
LOG_FREQUENCY_STEPS = int(TOTAL_STEPS / 500)

# Load in Amber files (prepared using AmberTools)
prmtop = AmberPrmtopFile('TTSCP_L102C_DDSB_preenergyminimized.prmtop')
inpcrd = AmberInpcrdFile('TTSCP_L102C_DDSB_preenergyminimized.inpcrd')
#pdb_file = 'TTSCP_L102C_DSB_energyminimized.pdb'
#pdb = PDBFile(pdb_file)
#k_rmsd = 100.0*kilocalories_per_mole/angstroms**2
#rmsd_receptor_0 = 0
#backbone_receptor=list(range(1538, 1588))
#print(backbone_receptor)
# Set up simulation "environment"
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        constraints=HBonds)
system.addForce(MonteCarloBarostat(1*bar, 298*kelvin))
#rmsd_cv =RMSDForce(pdb.positions, backbone_receptor) # RMSD biasing/harmonic force to act on the backbone of the 'receptor' chains
#rmsd_receptor = CustomCVForce('0.5*k_rmsd*(rmsd-rmsd_target)^2')
#rmsd_receptor.addCollectiveVariable('rmsd', rmsd_cv)
#rmsd_receptor.addGlobalParameter('rmsd_target', rmsd_receptor_0)
#rmsd_receptor.addGlobalParameter('k_rmsd', k_rmsd)
#system.addForce(rmsd_receptor)
# Set up integrator to integrate forces
integrator = LangevinIntegrator(298*kelvin, 1/picosecond, TIME_STEP)

# Pull together all the elements for the simulation
simulation = Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

# Energy minimisation
simulation.minimizeEnergy()

# Setup reportersi

simulation.reporters.append(DCDReporter('TTSCP_L102C_DDSB_energyminimization.dcd', LOG_FREQUENCY_STEPS))
simulation.reporters.append(
    StateDataReporter(
        stdout,
        LOG_FREQUENCY_STEPS,
        step=True,
        totalSteps=TOTAL_STEPS,
        potentialEnergy=True,
        temperature=True,
        kineticEnergy=True,
        speed=True,
        volume=True,
        progress=True
    )
)
simulation.reporters.append(
    StateDataReporter(
        'TTSCP_L102C_DDSB_energyminimization.txt',
        LOG_FREQUENCY_STEPS,
        step=True,
        totalSteps=TOTAL_STEPS,
        potentialEnergy=True,
        temperature=True,
        kineticEnergy=True,
        speed=True,
        volume=True,
        progress=True
    )
)
simulation.reporters.append(CheckpointReporter('TTSCP_L102C_DDSB_energyminimization.chk', 5000))

# Run the simulation!
simulation.step(TOTAL_STEPS)


# In[ ]:




