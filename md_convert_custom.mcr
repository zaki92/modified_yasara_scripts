# YASARA MACRO
# TOPIC:       3. Molecular Dynamics
# TITLE:       Convert between Sim, XTC, MDCrd and PDB simulation trajectories
# REQUIRES:    Dynamics 9.5.10
# AUTHOR:      Elmar Krieger
# LICENSE:     GPL
# DESCRIPTION: This macro converts an existing MD trajectory between various formats. Supported are conversions between YASARA Sim, GROMACS XTC and AMBER MDCrd trajectories, as well as conversion to PDB files

# Parameter section - adjust as needed
# ====================================

# The trajectory to convert must be present with a .sim, .xtc or .mdcrd extension.
# The starting scene *_water.sce is also required.
# You can either set the target by clicking on Options > Macro > Set target,
# by providing it as command line argument (see docs at Essentials > The command line),
# or by uncommenting the line below and specifying it directly.
#MacroTarget = 'c:\MyProject\1crn'

# Source format (srcformat) can be 'sim' (see SaveSim/LoadSim), 'xtc' (see SaveXTC/LoadXTC)
# or 'mdcrd' (see SaveMDCrd/LoadMDCrd).
# Destination format (dstformat) can be 'sim', 'xtc', 'mdcrd', 'pdb' (a series of PDB files)
# or 'pdbw' (a series of wrapped PDB files, where all atoms are inside the cell
# and potentially wrapped around periodic boundaries (i.e. broken molecules)).
# If one is left empty, YASARA will ask for the formats interactively.
srcformat='sim'
dstformat='pdb'

# Flag if water object should be included (1) or not (0)
waterincluded=0

# Forcefield to use 
ForceField AMBER14

# In case the simulation was run without 'CorrectDrift on' and the solute diffused
# through a periodic boundary, you can keep it centered here by specifying the number
# of an atom close to the core of the solute, which will be kept at the cell center.
# If both srcformat and dstformat are 'sim', you can correct an existing trajectory.
central=0

# If you want to convert only every Nth snapshot, set 'skippedsnapshots' to N-1.
# E.g. skippedsnapshots=1 will skip one snapshot for each converted one, and thus
# convert every 2nd snapshot.
skippedsnapshots=9

# No changes required below this point!

# Do we have a target?
if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target structure"
# Do we have source and destination format?
if srcformat=='' or dstformat==''
  # No, ask user interactively
  srcformat,dstformat,waterincluded,skippedsnapshots =
    ShowWin Type=Custom,Title="Select formats to convert trajectory",Width=600,Height=400,
            Widget=Text,        X= 20,Y= 50,Text="Choose the source format (the existing trajectory)",
            Widget=RadioButtons,Options=3,Default=1,
                                X= 20,Y= 70,Text="YASARA _S_im format",
                                X= 20,Y=106,Text="GROMACS _X_TC format",
                                X= 20,Y=142,Text="AMBER _M_DCrd format",
            Widget=Text,        X= 20,Y=188,Text="Choose the destination format (the trajectory to create)",
            Widget=RadioButtons,Options=5,Default=2,
                                X= 20,Y=208,Text="YASARA S_i_m format",
                                X= 20,Y=244,Text="GROMACS X_T_C format",
                                X= 20,Y=280,Text="AMBER M_D_Crd format",
                                X= 20,Y=316,Text="A series of _P_DB files",
                                X= 20,Y=352,Text="_W_rapped PDB files (bonds may cross periodic boundaries)",
            Widget=CheckBox,    X=310,Y=106,Text="_I_nclude water object",Default=Yes,
            Widget=Text,        X=310,Y=228,Text="Skipped snapshots per converted",
            Widget=NumberInput, X=310,Y=248,Text="_s_napshot",Default=0,Min=0,Max=1000,
            Widget=Button,      X=548,Y=348,Text="_O_ K"
  # Convert format from integer to string
  formatlist='sim','xtc','mdcrd','pdb','pdbw'
  srcformat=formatlist(srcformat)
  dstformat=formatlist(dstformat)
if srcformat!='sim' and srcformat!='xtc' and srcformat!='mdcrd'
  RaiseError "The source format must be either 'sim', 'xtc' or 'mdcrd'"
if srcformat==dstformat
  if srcformat=='xtc'
    RaiseError "A conversion from 'xtc' to 'xtc' is not possible"
  elif srcformat=='mdcrd'
    RaiseError "A conversion from 'mdcrd' to 'mdcrd' is not possible"
  elif not central
    RaiseError "A conversion from 'sim' to 'sim' only makes sense if an atom should be kept centered in the cell ('central')"
# Speed up conversion using short dummy cutoff and no longrange forces
Cutoff 2.62
Longrange None
# Load starting structure
LoadSce (MacroTarget)_water
bnd = Boundary
if waterincluded
  selection='Atom all'
else
  selection='Atom all and Obj !Water'
if dstformat=='pdb' or dstformat=='pdbw'
  # Since we can only save one object per PDB file, join all objects with atoms to the first with atoms
  JoinObj (selection),Atom 1
# Number of first snapshot (counting starts at 0 also for XTC)
first=00000
# Fix all atoms: after saving a snapshot, we proceed by one simulation step.
# If the simulation has a problem, this may blow it up. Fixing atoms avoids this.
FixAll
# Backwards compatibility: Starting with YASARA version 12.8.1, XTC trajectories no longer contain a number in the filename
old = FileSize (MacroTarget)00000.xtc
if old
  RenameFile (MacroTarget)00000.xtc,(MacroTarget).xtc
# Load the first two snapshots to calculate the saving interval (steps)
for i=0 to 1
  if srcformat=='sim'
    LoadSim (MacroTarget)(first+i)
  else
    Load(srcformat) (MacroTarget),(first+i+1)
  t(i) = Time
steps=0+(t1-t0)
TimeStep 1,1
if dstformat=='xtc' or dstformat=='mdcrd'
  # SaveXTC/SaveMDCrd does currently only extend but not overwrite an existing trajectory
  DelFile (MacroTarget).(dstformat)
  Save(dstformat) (MacroTarget),Steps=(steps*(1+skippedsnapshots)),(selection)
elif dstformat=='sim'
  SaveSim (MacroTarget)00000,Steps=(steps*(1+skippedsnapshots))
else
  # When saving PDB files, we can save time by pausing
  Sim Pause
Console Off
# As saving is done automatically for sim, xtc and mdcrd, all we need to do is load
# the snapshots from the source trajectory. This will adjust the current
# simulation time, which in turn triggers the save events. The 'Wait 1' is
# essential to suspend the macro for one cycle, so that the simulation
# continues and can be saved. (Snapshots are only saved while the simulation is running!).
i=first
do
  if srcformat=='sim'
    LoadSim (MacroTarget)(i)
    found = FileSize (MacroTarget)(i+skippedsnapshots+1).sim
  else
    last = Load(srcformat) (MacroTarget),(i+1)
    found=!last
  t = Time
  ShowMessage 'Converting (srcformat) snapshot (i) from (srcformat) to (dstformat) format, time is (0+t/1000) ps, (steps) fs/snapshot...'
  # Just in case the trajectory is not continuous, adjust the time
  Time ((0.+i-first)*steps)
  if central
    # Keep a chosen atom at the center of the cell (Cell returns center as values 7-9) 
    _,_,_,_,_,_,cen() = Cell
    pos() = PosAtom (central)
    MoveAtom all,(cen-pos)
  if dstformat=='pdb' or dstformat=='pdbw'
    # To unwrap the soup, we need to stop the simulation
    Sim Off
    # By transfering the soup into the cell, we make sure that the current cell
    # is copied into the saved PDB file as CRYST1 (see SavePDB docs). This also
    # sets a new transformation history, so we don't have to turn off transformation.
    TransferObj Atom 1,SimCell,Local=Fix
    if dstformat=='pdb'
      # And transform back into a non-periodic cell to avoid shifts
      Boundary Wall
      # Wall boundaries require a cuboid cell. Get the cell center to define a cell
      # that keeps the original atom coordinates but is cuboid.
      _,_,_,_,_,_,cen() = Cell
      Cell (cen*2),90,90,90
    Sim On
    # SavePDB expects an object selection, so 'Atom 1' selects the object with the first atom
    SavePDB Atom 1,(i)(i/(1+skippedsnapshots))
    Sim Off
    Boundary (bnd)
    Time (t)
  # We need to proceed to save
  Wait 1
  if srcformat!='sim'
    # Check if XTC/MDCrd trajectory ends before the next used snapshot in case we have 'skippedsnapshots'
    j=1
    while found and j<=skippedsnapshots
      last = Load(srcformat) (MacroTarget),(i+j+1)
      found=!last
      j=j+1
  i=i+skippedsnapshots+1
while found    
HideMessage
Sim Off
FreeAll

# Exit YASARA if this macro was provided as command line argument in console mode and not included from another macro
if runWithMacro and ConsoleMode and !IndentationLevel
  Exit
