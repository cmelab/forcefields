import mbuild as mb
import foyer
import hoomd
from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield


esp_ff = foyer.Forcefield(forcefield_files="/bsuhome/mpaul/scratch/forcefields/xml_files/PIDTBT_nC16.xml")

# We have to add the underscore to the names manually if we are using foyer XML files without SMARTS definitions
PIDTBT_nC16 = mb.load("/bsuhome/mpaul/scratch/forcefields/mol2files/PIDTBT_nC16_typed.mol2")

for p in PIDTBT_nC16.particles():
    p.name = f"_{p.name}"

box = mb.fill_box(compound=PIDTBT_nC16, n_compounds=1, box=[10,10,10])
box_pmd = esp_ff.apply(box)

snapshot, forcefield, refs = create_hoomd_forcefield(box_pmd, auto_scale=True, r_cut=2.5)




#setting cpu and simulation 
cpu = hoomd.device.CPU()
sim = hoomd.Simulation(device=cpu,seed=0)
sim.create_state_from_snapshot(snapshot)

#setting the integrator
kt = 1.2
free_particle = hoomd.filter.Tags(tags=[0,1])  #letting hoomd know which particles to update the positions of
integrator = hoomd.md.Integrator(dt = 0.0001)
nvt = hoomd.md.methods.NVT(filter=hoomd.filter.All(), kT=1.2, tau=0.5) #what is tau?
integrator.forces = forcefield
integrator.methods.append(nvt)
sim.operations.integrator = integrator
sim.state.thermalize_particle_momenta(filter=free_particle, kT=kt)
# Set up GSD writer
gsd_writer = hoomd.write.GSD(
    trigger=hoomd.trigger.Periodic(int(2e2)),
    filename="traj_test.gsd",  #name the output file
    mode="wb"
)
sim.operations.writers.append(gsd_writer)


sim.run(1e4)
