 # YASARA MACRO
# TOPIC:       3. Molecular Dynamics
# TITLE:       Analyzing a molecular dynamics trajectory
# REQUIRES:    Dynamics
# AUTHOR:      Elmar Krieger
# LICENSE:     GPL
# DESCRIPTION: This macro analyzes a simulation and creates a table with energies and RMSDs from the starting structure, as well as the minimum energy structure and the time average structure with B-factors calculated from the root mean square fluctuations. If you want to do your own custom analysis, search for 'Example:'.

# The structure to analyze must be present with a .sce extension.
# You can either set the target structure by clicking on Options > Macro > Set target,
# or by uncommenting the line below and specifying it directly.
#MacroTarget = 'c:\MyProject\1crn'

# Forcefield to use (these are all YASARA commands, so no '=' used)
# Use YASARA2 in YASARA Structure to include a quality Z-score
ForceField AMBER14,SetPar=Yes

# Number of the main object whose RMSDs from the starting conformation will be calculated
# If the protein is an oligomer, check the documentation of the 'Sup' command at 'analyzing a simulation' to avoid pitfalls.
currobj=1

# The B-factors calculated from the root-mean-square fluctuations can be too large to fit them
# into the PDB file's B-factor column. Replace e.g. 1.0 with 0.1 to scale them down to 10%
bfactorscale=1.0

# Flag to save a PDB file of the solute snapshots for further analysis
pdbsaved=0 

# Selection of atoms to include for 'Calpha' RMSD calculation (includes C1* of nucleic acids to consider DNA/RNA)
casel='CA Protein or C1* NucAcid'

# Selection of atoms for which the dynamic cross-correlation matrix (DCCM) should be visualized.
# Here are some typical examples:
# dccmsel=''                - Don't calculate the DCCM, the default
# dccmsel='Atom CA Protein' - Calculate the DCCM for protein Calpha atoms
# dccmsel='Res Protein'     - Calculate the DCCM for protein residue centers
dccmsel=''

# Selection of atoms for which the radial distribution function (RDF) should be visualized
# Here are some examples:   
# rdfsel='' - Don't calculate the RDF, the default
# rdfsel='O Res HOH,O Res HOH,Bins=40,BinWidth=0.25'
#        - Calculate the RDF of water in 40 bins, each 0.25 A wide (thus up to 10 A).
# rdfsel='CG Res Asp 120,ND1 Res His 200,Bins=20,BinWidth=0.5'
#        - Calculate the RDF between two specific atoms in 20 bins, each 0.25 A wide
#          (thus again up to 10 A). Note that you may have to save more snapshots than
#          usually in md_run.mcr to avoid problems with sparse data and noisy RDF results.  
rdfsel=''

# First snapshot to be analyzed, increase number to ignore an equilibration period.
firstsnapshot=0

# All snapshots will be superposed on this reference snapshot to calculate RMSDs etc.
# The starting structure is snapshot 0. Having run the macro once, you can also change
# refsnapshot=X to refsnapshot='average' to superpose on the time average structure.
refsnapshot=0

# No change required below this point
# ===================================

# Do we have a target?
if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target structure"

Clear
Console Off
SurfPar Molecular=Numeric

# Load scene with water or other solvent
waterscene = FileSize (MacroTarget)_water.sce
solventscene = FileSize (MacroTarget)_solvent.sce
if waterscene
  LoadSce (MacroTarget)_water
elif solventscene 
  LoadSce (MacroTarget)_solvent
else
  RaiseError 'Could not find initial scene file (MacroTarget)_water.sce. You must run a simulation with the macro md_run first'

calphas = CountAtom (casel)
if calphas>0 and calphas<3
  # We cannot superpose 1 or 2 Calpha atoms
  casel='None'

ShowMessage "Preparing analysis, please wait..."
Wait 1

# See if structure validation checks should be done (require YASARA2 force field)
checked=0
fof = ForceField
if fof=='YASARA2'
  checked=1

# Backwards compatibility: Starting with YASARA version 12.8.1, XTC trajectories no longer contain a number in the filename
old = FileSize (MacroTarget)00000.xtc
if old
  RenameFile (MacroTarget)00000.xtc,(MacroTarget).xtc

# Determine trajectory format
for format in 'xtc','mdcrd','sim'
  found = FileSize (MacroTarget).(format)
  if found
    break
  
if refsnapshot=='average'
  # We superpose on the time average structure
  filename='(MacroTarget)_average.pdb'
  exists = FileSize (filename)
  if !exists
    RaiseError "No time-average structure has been calculated yet, cannot superpose onto it. Run the macro once with refsnapshot=0 (or any other number), then run again with refsnapshot='average'"
  startobj = LoadPDB (filename)
  SupObj (startobj),(currobj)
else
  if refsnapshot
    # We superpose on a certain snapshot
    if format=='sim'
      LoadSim (MacroTarget)(00000+refsnapshot)
    else
      Load(format) (MacroTarget),(refsnapshot+1)
  # Duplicate the intial object for RMSD calculation
  startobj = DuplicateObj (currobj)
RemoveObj (startobj)

i=00000+firstsnapshot
emin=1e99
last=0
while !last
  # Load next snapshot from SIM or XTC trajectory
  if format=='sim'
    sim = FileSize (MacroTarget)(i+1).sim
    if not sim
      last=1
    LoadSim (MacroTarget)(i)
  else
    # Set last (end of file) to 1 if last snapshot loaded
    last = Load(format) (MacroTarget),(i+1)
  Sim Pause
  # Add time in picoseconds to table
  simtime = Time
  ShowMessage 'Analyzing snapshot (0+i) at (0+(simtime/1000)) ps'
  Wait 1
  Tabulate (simtime/1000)
  # Get energy components, tabulate some and individual components 1-6 (exlucing packing energies)
  elist() = EnergyAll All
  Tabulate (sum elist)
  for j=1 to 6
    Tabulate (elist(j))
  # Prepare to save the minimum energy structure (ignoring solvent-solvent interactions)
  e = EnergyObj (currobj)
  #
  # The following examples provide a few hints for other things to analyze while
  # the simulation is active (i.e. considering periodic boundaries).
  # Just remove the two characters '# ' to uncomment the command(s) of each example.
  # Each value you tabulate becomes a new column in the result table, between energies and RMSDs.
  #
  # Example: Measure the distance between the carboxyl group of Glu 123 (Cdelta) in molecule/chain A
  #          and the guanidinium group of Arg 345 (Czeta) in molecule/chain B:
  # Tabulate Distance CD Res Glu 123 Mol A, CZ Res Arg 345 Mol B
  #
  # Example: Measure the angle between three Calpha atoms:
  # Tabulate Angle CA Res Glu 10, CA Res Asp 30, CA Res Lys 40
  #
  # Calculate the molecule surface area of the solute object 'currobj'
  Tabulate SurfMol E and Protein (currobj), accessible #calculates SASA for Mol B (change to desired MOL)
  Tabulate SurfMol B and Protein (currobj), accessible #calculates SASA for Mol C (change to desired MOL)
  Tabulate SurfMol B E and Protein (currobj), accessible #calculates SASA for atoms 1-7158 (change last number to include all protein atoms)

  #Tabulate SurfMol B and Protein (currobj), vdw #uncomment lines to get VDW Surface area
  #Tabulate SurfMol C and Protein (currobj), vdw
  #Tabulate SurfAtom 1-7158 (currobj), vdw

  #Tabulate SurfMol B and Protein (currobj), molecular #uncomment lines to get molecular surface area
  #Tabulate SurfMol C and Protein (currobj), molecular
  #Tabulate SurfAtom 1-7158 (currobj), molecular
  
  #WriteTable "SurfAtom 1-3000, accessible",'SASA Mol A','SASA Mol A'
  #  WriteTable "SurfAtom 3000-12399, accessible",'SASA Mol D','SASA Mol D'
  #  WriteTable "SurfAtom 1-12399, accessible", 'SASA Total', 'SASA Total'

  # Example: Count the number of atoms in molecule E that are hydrogen-bonded to H,L 
  hbolist() = ListHBoAtom Mol B and protein, Mol E and protein
  Tabulate (count hbolist)
  #
  # Example: Check if the carboxyl group of Glu 123 and the guanidinium group of Arg
  #          345 are currently hydrogen-bonded (1) or not (0). The last row in the
  #          table then contains the average, i.e. the probability of being H-bonded. 
  # hbolist() = ListHBoAtom Sidechain Res Glu 123, Sidechain Res Arg 345
  # Tabulate (count hbolist>0)
  #
  # Example: Count the number of salt-bridges involving Lys,Arg and Asp,Glu
  # lysbridges = CountRes Lys Atom NZ with distance<4 from Asp Glu Atom OD? OE?
  # argbridges = CountRes Arg Atom NE NH? with distance<3.5 from Asp Glu Atom OD? OE?
  # Tabulate (lysbridges+argbridges)
  #
  # Example: Count the number of residues involved in hydrophobic interactions with Phe 13
  # intlist() = ListIntRes Phe 13,all,Type=Hydrophobic,Exclude=5
  # Tabulate (count intlist)
  #
  # Example: Calculate the current potential energy of residue Glu 123:
  # Tabulate EnergyRes Glu 123 
  #
  # Example: Calculate the secondary structure content. This adds five values to the
  #          Table: the percentage of helix, sheet, turn, coil and 3-10 helix. 
  # Tabulate SecStr
  # 
  if checked
    # Calculate structure validation Z-scores, using the formula from homology modeling
    for type in 'dihedrals','packing1d','packing3d'
      zscore(type) = CheckObj (currobj),(type)
    zscore=zscoredihedrals*0.145+zscorepacking1d*0.390+zscorepacking3d*0.465 
  if rdfsel!=''
    # Collect data for radial distribution function
    BinDistance (rdfsel)
  Sim Off
  if e<emin
    # Save minimum energy structure
    emin=e
    SaveSce (MacroTarget)_energymin
    SavePDB (currobj),(MacroTarget)_energymin
  if last
    # Save last structure
    SaveSce (MacroTarget)_last
    SavePDB (currobj),(MacroTarget)_last  
  if pdbsaved
    # Save a PDB file of the solute
    SavePDB (currobj),(MacroTarget)_(i)
  # The following examples provide a few hints for other things to analyze while the
  # simulation is *not* active (periodic boundaries removed) but the starting structure
  # has not yet been added to the soup to perform RMSD calculations
  #
  # Example: Measure the radius of gyration of the main object
  # Tabulate RadiusObj (currobj),Center=Mass,Type=Gyration
  #
  # Example: Calculate the angle between the transmembrane helix formed by residues 35-65
  # in Mol A and the membrane normal (which is parallel to the Y-axis, since the membrane is
  # parallel to the XZ-plane) using formula angle=acos(dotpro(dir,Y)/(len(dir)*len(Y))) 
  # _,_,_,dir() = GroupLine CA Protein Res 35-65 Mol A
  # Tabulate (acos dir2)
  # 
  # Example: Calculate the angle between the secondary structure elements formed by
  # residues 106-140 and 149-169 in Mol B, in the range -180 to 180 (see GroupAngle docs)
  # Tabulate GroupAngle CA Protein Res 106-140 Mol B,CA Protein Res 149-169 Mol B,Range=360
  #
  # Example: Calculate the Coulomb and VdW binding energies between the solute 'currobj' and the solvent
  #          EBind = ESolute+ESolvent-ETotal. For ligand binding look at md_analyzebindenergy macro.
  # etotvdw,etotcoulomb = Energy VdW,Coulomb
  # RemoveObj Water
  # esoluvdw,esolucoulomb = Energy VdW,Coulomb
  # AddObj Water
  # RemoveObj (currobj)
  # esolvvdw,esolvcoulomb = Energy VdW,Coulomb
  # AddObj (currobj)
  # Tabulate (esoluvdw+esolvvdw-etotvdw),(esolucoulomb+esolvcoulomb-etotcoulomb)
  #
  # Add starting structure to the soup to perform RMSD calculations etc.
  AddObj (startobj)
  #
  # The following examples provide a few hints for other things to analyze while the
  # simulation is *not* active (i.e. things that can't be analyzed with periodic boundaries)
  # Just remove the two characters '# ' to uncomment the command(s) of each example.
  # Each value you tabulate becomes a new column in the result table, between energies and RMSDs.
  #
  # Example: Calculate the backbone RMSD for residues 1-120:
  #Tabulate SupAtom Backbone Mol A Obj (currobj),Backbone Mol A Obj (startobj)
  #Tabulate SupAtom Backbone Mol E Obj (currobj),Backbone Mol E Obj (startobj)
    
#Tabulate SupAtom Backbone Res 3569-3632 Obj (currobj),Backbone Res 3569-3632 Obj (startobj)
  #
  # Example: Calculate the sidechain heavy-atom RMSD:
  # Tabulate SupAtom Sidechain Element !H Obj (currobj),Sidechain Element !H Obj (startobj)
  #
  # Example: Measure how far residue HP6 86 moved since the start (after superposing on Calphas) 
  # SupAtom (casel) and Obj (currobj),(casel) and Obj (startobj)
  # Tabulate GroupDistance Res HP6 86 Obj (currobj),Res HP6 86 Obj (startobj) 
  #
  # Example: Measure the distance between two centers of mass, e.g. the loop from
  #          residue Ala 205 to Glu 210, and the ligand NAD: 
  # cenA() = GroupCenter Res Ala 205 - Glu 210 Obj (currobj), Type=Mass
  # cenB() = GroupCenter Res NAD Obj (currobj), Type=Mass
  # Tabulate (norm (cenA-cenB))
  #
  # Add CA, backbone and heavy atom RMSDs to table
  # Trick: If the solute object contains neither CA nor backbone atoms,
  # we simply assign the SupAtom error code (=0)
  carmsd = SupAtom (casel) and Obj (currobj),(casel) and Obj (startobj)
  bbrmsd = SupAtom Backbone Obj (currobj),Backbone Obj (startobj)
  aarmsd = SupAtom Element !H Obj (currobj),Element !H Obj (startobj)
  # Add the Calpha, backbone and all-atom RMSDs to the table
  Tabulate (carmsd),(bbrmsd),(aarmsd)    
  if checked
    Tabulate (zscore)
  # Add the current atom positions to internal table to obtain RMSF and average positions
  AddPosAtom Obj (currobj)
  RemoveObj (startobj)
  # Next snapshot
  i=i+1

if i==firstsnapshot
  RaiseError "This macro is meant to analyze a molecular dynamics trajectory created with md_run, but none was found in this directory"

# Calculate how many columns the table has, maybe the user tabulated private analysis data
vallist() = Tab Default
columns = count vallist/(i-firstsnapshot)
# Create a new final table with the now known number of columns
MakeTab Final,Dimensions=2,Columns=(columns)
# Put the original table 'Default' into the new one
Tabulate Tab Default 
# Loop over the columns (except 1, the time) to calculate the average values
Tabulate 'Average'
for i=2 to columns
  vallist() = Tab Final,Column=(i)
  Tabulate (mean vallist)
# Create a table header
header='____Time[ps] Energy[(EnergyUnit)]_____Bond _______Angle ____Dihedral ___Planarity _____Coulomb _________VdW '
# You can add your own headers here, 12 characters and a space per column 
for i=1 to columns-11-checked
  header=header+'_Your result '
header=header+'_RMSDs[A]:CA ____Backbone __HeavyAtoms'
if checked
  header=header+' QualityScore'
SaveTab Final,(MacroTarget)_analysis,Format=Text,Columns=(columns),NumFormat=12.3f,(header)

# Calculate time-average structure and set B-factors according to RMSF
AveragePosAtom Obj (currobj)
# The dummy assignment '_ =' ensures that B-factors are not printed to the console
_ = RMSFAtom Obj (currobj),Unit=BFactor
if bfactorscale!=1.0
  # Scale B-factors so that they fit into the PDB format
  first,last=SpanAtom Obj (currobj)
  for i=first to last
    bf = BFactorAtom (i)
    BFactorAtom (i),(bf*bfactorscale)    
if refsnapshot!='average'
  # The time average structure has incorrect covalent geometry and should be energy minimized
  SavePDB (currobj),(MacroTarget)_average

"""
# Example: Create an additional RMSF file, in case B-factors are too large for the PDB format
DelTab default
first,last = SpanAtom Obj (currobj)
rmsflist() = RMSFAtom Obj (currobj)
for i=first to last
  res = ListAtom (i),Format='RESNAME,RESNUM'
  Tabulate (i),(res),(rmsflist(i-first+1))
SaveTab default,(MacroTarget)_rmsf,Format=Text,Columns=4,NumFormat=8.2f,"RMSF Table"  
"""  

if rdfsel!=''
  # Additionally calculate and show the radial distribution function (RDF).
  MakeTab RDF,Dimensions=1
  Tabulate RDF
  SaveTab RDF,(MacroTarget)_rdf,Format=Text,Columns=1,NumFormat=6.3f,
          'Radial distribution function with parameters (rdfsel)'
  obj = ShowTab RDF,Width=1.0,Range=10,MinCol=Blue,MaxCol=Yellow
  PosObj (obj),X=0,Y=-12,Z=31
  RotateObj (obj),X=-90
  SaveSce (MacroTarget)_rdf
  DelObj (obj)

if dccmsel!=''
  # Additionally calculate and show the dynamic cross-correlation matrix (DCCM).
  # This matrix correlates the displacements from the time average structure,
  # see the documentation of the 'DCCM' command for details.
  # First get the number of selected units, i.e. the rows/columns in the matrix 
  units = Count(dccmsel)
  if !units
    RaiseError 'The DCCM selection (dccmsel) did not match any atoms'
  # Take the time average structure as the start object to superpose onto
  DelObj (startobj)
  startobj = DuplicateObj (currobj)
  RemoveObj (startobj)
  # Loop over the snapshots a second time to calculate the displacements from the time average
  i=00000+firstsnapshot
  while 1
    # See if next snapshot is present
    sim = FileSize (MacroTarget)(i).sim
    if not sim
      break
    # Yes, load it
    LoadSim (MacroTarget)(i)
    Sim Pause
    ShowMessage 'Calculating dynamic cross-correlation matrix, analyzing snapshot (0+i)...'
    Wait 1
    Sim Off
    # Superpose snapshot on the time average structure
    AddObj (startobj)
    SupAtom (casel) and Obj (currobj),(casel) and Obj (startobj)
    # Add the current displacements to an internal table to obtain the DCCM
    AddDisp(dccmsel) Obj (currobj),(dccmsel) Obj (startobj)
    RemoveObj (startobj)
    # Next snapshot
    i=i+1
  # Store the DCCM in a table
  MakeTab DCCM,Dimensions=2,Columns=(units)
  Tabulate DCCM
  # Visualize the DCCM
  pointwidth=1.
  height=5.
  dccmobj1 = ShowTab DCCM,Width=(pointwidth),Range=(height),Min=-1,Max=1.0
  # By default, ShowTab shows the minimum at Z=0, move so that correlation 0 is at Z=0
  MoveObj (dccmobj1),Z=(height*0.5)
  # Visualize the zero level with a flat DCCM wireframe
  dccmobj2 = ShowTab DCCM,Width=(pointwidth),Range=0,Min=-1,Max=1.0
  dccmobj3 = ShowWireObj (dccmobj2),Static,Mesh=Solid
  DelObj (dccmobj2)
  PointPar Radius=0.5,Plastic=No
  NameObj (dccmobj3),ZeroLevel
  RotateObj (dccmobj1) (dccmobj3),X=180
  # Create a text object with the residue names and the table header
  textwidth=pointwidth*units*2
  textobj1 = MakeTextObj Units,Width=(textwidth),Height=(textwidth)
  Font Arial,Height=(pointwidth*0.6),Color=Yellow,Depth=0.5,DepthCol=Red
  idlist() = List(dccmsel) Obj (currobj),Format='MOLNAME RESName RESNUM' 
  for i=1 to units
    PosText X=(textwidth*0.5+pointwidth*-0.5*(units+6)),
            Y=(textwidth*0.5+pointwidth*(0.5*units-i)),justify=left
    Print (idlist(i))
  # Duplicate the labels at the bottom
  textobj2 = DuplicateObj (textobj1)
  RotateObj (textobj2),Z=90
  DelObj not (dccmobj1) (dccmobj3) (join textobj)
  RenumberObj all
  # Save the matrix
  SaveTab DCCM,(MacroTarget)_dccm,Format=Text,Columns=(units),NumFormat=6.3f,
          'Dynamic Cross-Correlation Matrix for (units) selected units'
  # Save the visualized matrix
  SaveSce (MacroTarget)_dccm
HideMessage

# Exit YASARA if this macro was provided as command line argument in console mode
if runWithMacro and ConsoleMode
  Exit
