# YASARA MACRO
# TOPIC:       3. Molecular Dynamics
# TITLE:       Analyzing the ligand binding energy during a molecular dynamics simulation
# REQUIRES:    Dynamics
# AUTHOR:      Elmar Krieger
# LICENSE:     GPL
# DESCRIPTION: This macro calculates the binding energy of a ligand during a molecular dynamics simulation, including the time average. By default, the receptor must be object 1, the ligand must be object 2 (or change the settings). More positive energies indicate better binding, negative energies DO NOT indicate no binding, see the 'BindEnergy' command in the user manual for details, also for words of caution, a much more sophisticated and reliable approach is in preparation. 

# The structure to analyze must be present with a .sce extension.
# You can either set the target structure by clicking on Options > Macro > Set target,
# by providing it as command line argument (see docs at Essentials > The command line),
# or by uncommenting the line below and specifying it directly.
# MacroTarget = '/home/myprojects/1h8a/1h8a'

# Forcefield to use (these are all YASARA commands, so no '=' used)
# Use YASARA2 in YASARA Structure to include a quality Z-score
ForceField AMBER14,SetPar=Yes

# Set the object numbers of receptor and ligand for binding energy calculation.
# If the ligand is not in a separate object, set ligobj also to 1 and specify
# the ligand residue with 'ligres' further below:
recobj=1
ligobj=1
# If the ligand is part of the receptor object, specify its residue ID(s),
# e.g. ligres='ATP 301' for a single residue ATP, ligres='ATP 301 Mol A' if the ligand ATP 301 is
# present multiple times in an oligomer to choose only the one in molecule A, or
# ligres='Thr 32-Ala 35' for a peptide, or ligres='Mol C' if the ligand has a unique molecule name C.
ligres='Mol D'

# First snapshot to be analyzed, increase number to ignore an equilibration period.
firstsnapshot=0

# Choose method, either boundary elements ('BoundaryFast') or Poisson-Boltzmann ('PBS').
# In the latter case you get almost 'MM/PBSA', just without the entropy term from normal mode analysis.
# When comparing with other MM/PBSA calculations, make sure to use the same 'surfcost' below.
# Note that MM/PBSA requires a cuboid cell and does not work with dodecahedral/triclinic cells.
method='PBS'

# Temperature to use if method='PBS' above, 'BoundaryFast' is not temperature-dependent.
Temp 298K

# This is a guesstimate of the entropic cost of exposing one A^2 to the solvent in kJ/mol
# (The water molecules in the first solvent shell have to reorganize and lose
#  conformational freedom, which decreases their entropy and thus costs energy.
#  The value below does not include VdW interaction energies and is thus larger
#  than e.g. the 0.03 kJ/(mol*A^2) sometimes mentioned in textbooks). 
# We divide by 6.02214199e20 to obtain J and multiply with JToUnit
# to get the currently selected energy unit
surfcost=(0.65e0/6.02214199e20)*JToUnit
  
# No change required below this point
# ===================================

RequireVersion 18.11.10

# Do we have a target?
if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target structure"

Clear
Console Off
# Ensure that the surface resolution is high enough for solvation calculations
SurfPar Resolution=3,Molecular=Numeric
  
# Load scene with water or other solvent
waterscene = FileSize (MacroTarget)_water.sce
solventscene = FileSize (MacroTarget)_solvent.sce
if waterscene
  LoadSce (MacroTarget)_water
elif solventscene 
  LoadSce (MacroTarget)_solvent
else
  RaiseError 'Could not find initial scene file (MacroTarget)_water.sce. You must run a simulation with the macro md_run first'

# Check that recobj and ligobj are OK
objlist() = ListObj all
for type in 'rec','lig'
  if (type)obj not in objlist
    RaiseError '(type)obj is currently set to ((type)obj), but this object is not defined, please adapt the (type)obj setting in this macro'
  objatoms = CountAtom Obj ((type)obj)
  if !objatoms
    RaiseError '(type)obj is currently set to ((type)obj), but this object does not contain any atoms, please adapt the (type)obj setting in this macro'

# Check if cell needs to be changed from triclinic to cuboid to run a Poisson-Boltzmann calculation
cellchanged=0
cellinfo() = Cell
if method=='PBS' and (cellinfo(3)!=90 or cellinfo(4)!=90 or cellinfo(5)!=90)
  cellchanged=1
cellpos() = PosObj SimCell

# Look for membrane
memobjects = CountObj Membrane
if memobjects>1
  RaiseError 'Found (memobjects) membrane objects while at most one was expected'
  
if ligres!=''
  # Ligand is part of receptor object and needs to be extracted, perform consistency checks
  ligresidues = CountRes (ligres)
  if ligresidues==0
    RaiseError 'No ligand residues matching (ligres) were found'
  # Make sure that the object is chosen correctly
  if ligobj!=recobj or recobj!=1
    RaiseError 'If a ligand residue is selected, ligand and receptor must both be in object 1, please set recobj and ligobj to 1'
  if Objects!=3+memobjects
    RaiseError 'The scene must contain 3-4 objects: Receptor/Ligand, SimCell, Membrane [optionally] and Water'
  # Check that it's a single continuous selection
  ligatoms = CountAtom Res (ligres)
  first,last = SpanAtom Res (ligres)   
  if last-first+1!=ligatoms
    ListRes (ligres)
    RaiseError 'The selection (ligres) matches (ligatoms) atoms, but these are not continuous in the soup, ranging from atom (first) to (last) instead' 
  # Make sure the ligand is not covalently bound
  boundatoms = CountAtom not Res (ligres) with bond to Res (ligres)
  if boundatoms
    RaiseError 'The ligand residue (ligres) is bound to (boundatoms) other atoms, YASARA cannot calculate binding energies of covalently bound ligands'
  # Extract ligand residues
  splitpoints=''+first
  if last+1<=Atoms
    splitpoints=splitpoints+' (last+1)'
  SplitAtom (splitpoints)
  objlist() = SplitObj (recobj),Center=Yes,Atom (splitpoints),Keep=AtomNum
  ligobj=2
  if count objlist==3
    recobj='1 3'
elif ligobj==recobj
  RaiseError 'Both recobj and ligobj are set to (ligobj), but the ligand residue selection ligres is empty. Please correct it'
  
ShowMessage "Preparing analysis, please wait..."
Wait 1

# Create a table with 8 columns
MakeTab Final,Dimensions=2,Columns=8

# Determine trajectory format
for format in 'xtc','mdcrd','sim'
  found = FileSize (MacroTarget).(format)
  if found
    break
  
i=00000+firstsnapshot
emin=1e99
last=0
while !last
  # Load next snapshot from SIM, XTC or MDCRD trajectory
  if format=='sim'
    sim = FileSize (MacroTarget)(i).sim
    if not sim
      break
    LoadSim (MacroTarget)(i)
  else
    # Set last (end of file) to 1 if last snapshot loaded
    last = Load(format) (MacroTarget),(i+1)
  Sim Pause
  # Add time in picoseconds to table
  simtime = Time
  ShowMessage 'Analyzing snapshot (0+i) at (0+(simtime/1000)) ps'
  Wait 1
  for type in 'cmp','rec','lig'
    # Keep only the relevant objects
    if type=='cmp'
      RemoveObj not (recobj) (ligobj)
      if cellchanged
        # Choose a cuboid cell for Poisson-Boltzmann 
        Cell Auto,Extension=15
    else
      RemoveObj not ((type)obj)
    # Energies are calculated with net-charge 0 to avoid overestimation of net-charge
    # effects when comparing binding energies of molecules with different charges.
    ChargeObj all,0
    # Start a simulation so that receptor and ligand are transferred into the coordinate
    # system of the simulation cell and the surface area is not simply the sum of both.
    Sim on
    epot(type) = Energy
    esolcoulomb,esolvdw = SolvEnergy (method)
    molsurf = Surf molecular
    esol(type)=esolcoulomb+esolvdw+molsurf*surfcost
    print 'epot(type)=(epot(type)), (esolcoulomb),(esolvdw),surfcost=(molsurf*surfcost), esol(type)=(esol(type)), total=(epot(type)+esol(type))'
    AddObj all
    # Restore original cell position, shape will be set by LoadSim
    PosObj SimCell,(cellpos)
  # Calculate result: energy of separated compounds - energy of complex
  Tabulate (simtime/1000),(epotrec+esolrec+epotlig+esollig-epotcmp-esolcmp),(epotrec),(esolrec),(epotlig),(esollig),(epotcmp),(esolcmp)
  # Next snapshot
  i=i+1

if i==firstsnapshot
  RaiseError "This macro is meant to analyze a molecular dynamics trajectory created with md_run, but none was found in this directory"

# Tabulate average binding energy
Tabulate 'Average'
vallist() = Tab Final,Column=2
Tabulate (mean vallist)
# Create a table header
header='____Time[ps]__BindEnergy=  EpotRecept+ EsolvRecept+  EpotLigand+ EsolvLigand- EpotComplex- EsolvComplex [(EnergyUnit)]'
SaveTab Final,(MacroTarget)_bindenergy,Format=Text,Columns=8,NumFormat=12.3f,'(header)'
HideMessage

# Exit YASARA if this macro was provided as command line argument in console mode and not included from another macro
if runWithMacro and ConsoleMode and !IndentationLevel
  Exit
