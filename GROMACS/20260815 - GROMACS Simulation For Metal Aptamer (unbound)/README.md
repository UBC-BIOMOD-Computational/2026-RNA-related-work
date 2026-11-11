# GROMACS SIMULATION FOR METAL APTAMER (unbound to metal ion ligand)

Purpose: Simulation of lead (pb) RNA aptamer in GROMACS

## Table of Contents

## Simulation

### Simulation step 1: obtain RNA 3d structural file

- You could technically get this by taking a sequence, running it in a secondary structure prediction model, then running the secondary structure in a tertiary structure model but we already have an experimental result for our sequence https://www.rcsb.org/structure/1LDZ 
    - So why reinvent the wheel and risk being inaccurate?? We should just use official docs

- Whenever you download items from pdb- always check where the experimental structure was derived from. Like was it NMR, was it xray crystallography?  
    - For us, it was NMR: see this paper, which has the exact pdb id tied to it. https://www.sciencedirect.com/science/article/abs/pii/S0022283698921829?via%3Dihub 

- And then check clashscores - it should be blue not red

Result: get the .pdb or .cif file from the website.
- If you don’t see the file extension you will have to turn on “Show Hidden Extensions” in your file system. If this is confusing ask Google. 

Inspect: inspect pdb file 
- Commands: `show all models`
- Commands: `hide all models`



-----------------------------------------------------------------

### Simulation step 2: cleanup pdb file

#### Step 2.1 Keep one NMR model (typically Model 1)
- save model1.pdb models #1.1 or close #1.2-1.25
    -	save via:
        - UI
        - the command line
    
    
    - Note: You may have to remove `ENDMDL` => this a pdb file note when one pdb file has multiple models, so if you only have one model, it should be fine, but just double check by opening pdb file in text editor and ctrl-f for ENDMDL


#### Step 2.2 Remove waters and lead ion. 

- Why do we have to remove Water and metal ions?
    - Ans: Force fields!!! 
    - For Molecular Dynamics Simulations, you generally re-solvate the entire system using your chosen water model, rather than blindly retaining every PDB water.
    - However, structurally important waters may be deliberately retained if there's good evidence they are important.

Commands
- `show atoms`
- `show solvent`
- `select :C` selects Carbon ions
- `select :Pb` selects Lead ions

Do we have any water or metal ions?

Note: water will usually show in X-Ray or Cryo EM structures:
- NMR structure will not have crystallographic waters in the same sense as an X-ray crystal structure.

#### Step 2.3 


- Check all atoms are nucleic acid residues and names
    - Are there any missing backbones?
- `select nucleic`
- `info residues`
    - visually inspect for break
- `clashes`
    - detects two atoms are closer together than expected based on their van der Waals radii.
    - Can be an issue if you have a lot of clashes
- `hide all` -> `clashes reveal true`
- `show atoms`

You can also open the pdb file

### Step 3 Prep System

#### Step 3.1 - Check you have GROMACS, correct force fields

0. open up your ubuntu shell specifically (windows). or terminal shell (mac)
1. `gmx --version` : do you have GROMACS installed
2. `gmx --version | grep "Data prefix"`: Find which folder is GROMACS is installed on your computer 
3. `ls /usr/local/gromacs/share/gromacs/top`: This is the topology folder
    - Note: replace `/usr/local/gromacs` with the result of step 2 above
    - macOS: the GROMACS directory is often somewhere under Homebrew, such as `/opt/homebrew/...` on Apple Silicon or `/usr/local/...` on Intel Macs, so don't assume the same path.
    
4. Search specifically that we have RNA-capable force field
```bash
find /usr/local/gromacs/share/gromacs/top \
  -name "rna.rtp" \
  -o -name "nucleic.rtp" \
  -o -name "*ol3*" \
  -o -name "*bsc1*"
```
The standard amber99sb-ildn is primarily a protein force field. The RNA improvements (parmbsc0, parmbsc1, OL3, etc.) are distributed as modified AMBER force fields. To use them, they must already exist in your GROMACS share/gromacs/top directory or be installed separately.


#### Step 3.2 - CONVERT PDB TO GMX SYSTEM: 0 Get force field files - AKA getting topology
OUTPUT: `pdb2gmx` means GROMACS successfully converted your RNA into a force-field-compatible structure/topology

1. Go into the folder with your model in the command line
    - if you are on windows WSL, you will have to get into your mounted folde
    - `cd /mnt/c` => bring you to your windows folder
    - continue navigating the file system until you get to the folder

2. Run the following: -> replace the line `-f 1LDZ_model_1.pdb` with your file
```bash 
gmx pdb2gmx \
    -f 1LDZ_model_1.pdb \
    -o processed.gro \
    -p topol.top \
    -i posre.itp
```
- pdb2gmx is a GROMACS program that converts your molecular structure into a form GROMACS can simulate.

3. Then choose the line that says
```bash
AMBER OL3 or ff99bsc0χOL3.
```

This is called getting topology and will produce .top file. expected output files:
•	processed.gro = Your processed 3D structure after pdb2gmx 
•	topol.top  = The system topology
    - This is telling simulations what kinds of forces it have to know:
        - Which force-field parameters are being used
            - ie: Bonds/interactions, Molecule types, Number of molecules, RNA, water, ions, etc.
•	posre.itp = position-restraint file
    - use `head posre.itp` or `head <any file name>` to see the first few lines of each file
    - in our case the fx, fy and fz is the restraint strength in the x, y, and z directions.
    - 1000 is pretty strong
        - It means the atom is strongly encouraged to stay near its starting position.

Troubleshooting:
- `gmx pdb2gmx -f 1LDZ_model_1.pdb -o processed.gro -p topol.top -i posre.itp -ignh`

#### Step 3.3: prep system: Prepare simulation box 
Output: => boxed.gro

```
gmx editconf \
    -f processed.gro \
    -o boxed.gro \
    -c \
    -d 1.0 \
    -bt dodecahedron
```
- `-c` put the RNA at the center of the box
- `-d`	1.0 nm padding (leave 1.0 nm between the RNA and the box boundary)
- `bt` dodecahedron box




### Simulation: prep system = solvate and neutralize 
Add back waters and neutralize for rna bases
SOLVENT

Find tip3
find /usr/local/gromacs/share/gromacs -name "tip3p.gro" -o -name "spc216.gro"
/usr/local/gromacs/share/gromacs/top/spc216.gro

gmx solvate -cp boxed.gro -cs tip3p.gro -o solv.gro -p topol.top

gmx solvate -cp boxed.gro -cs /usr/local/gromacs/share/gromacs/top/spc216.gro -o solv.gro -p topol.top

got solv.gro and topol.top


IONS
First make a tpr
gmx grompp \
    -f ions.mdp \
    -c solv.gro \
    -p topol.top \
    -o ions.tpr

You need empty ions.mdp

Then
gmx genion \
    -s ions.tpr \
    -o solv_ions.gro \
    -p topol.top \
    -neutral \
    -conc 0.15
Select
SOL
when prompted.
For RNA,
0.15 M NaCl is a reasonable physiological concentration unless you're reproducing a specific experiment with different salt conditions




Why add ions at all?
Because RNA is highly negatively charged. Every phosphate on the backbone carries negative charge. If you only:
1.	solvate 
2.	minimize 
then you will have a box of water around a strongly charged RNA with no counterions. That is usually not a realistic physical environment for RNA MD.
Adding ions does two things:
•	Neutralizes the system so the long-range electrostatics are physically sensible. 
•	Optionally adds salt concentration like 0.15 M, which better mimics solution conditions. 
Without ions, you may get:
•	overly strong repulsion between phosphates, 
•	unrealistic RNA collapse or expansion, 
•	unstable equilibration, 
•	poor comparison to experiment. 
Can you skip ions? sure


### Simulation: prep system = energy minimization

gmx grompp \
    -f minim.mdp \
    -c solv_ions.gro \
    -p topol.top \
    -o em.tpr

gmx mdrun -deffnm em
Check
Potential energy

Maximum force

No LINCS warnings


### Simulation: prep system = heat up to stable equilibrium phase to prepare for simulation
NVT
Heat from approximately
0 K

↓

300 K
using position restraints.
Example
100–500 ps

NPT

Maintain
300 K
1 bar
Typical
500 ps–2 ns
while gradually relaxing restraint


Simulation: run simulation finally
Typical first run
100 ns
Longer if studying RNA folding or conformational changes.


Okay results:
### step Analyze results 



## Extra Information

### Resolution

Resolution is a measure of how much structural detail an experimental method can reliably distinguish or determine.

HOWEVER, this is very different between different experimental methods.


- For X-ray & cryo-EM:
    - Resolution = how finely spatial features can be distinguished in the experimental map.
    - Important: 1.7 Å resolution does NOT mean atoms are 1.7 Å apart.
- NMR:
    - No equivalent single map-resolution number.


-	Resolution tells you the smallest spacing between two features that can still be distinguished
-	1.7 =/= saying atoms are 1.7A apart
o	C–C bond is about 1.54 A, but 1.7 can determine carbon positions because 1.7 is about minimum distinguishable spacing in the reconstructed electron density

-	X ray
    - X-RAY RANGE: angstroms
        - 3.5A = alpha helix beta sheets
        - 2.5A = side chains mostly, you can get back bone
        - 2A = water appears
        - 1.5 = you can get atomic detail

    - Method:
        - X-rays shot a 3d lattice of proteins => wow you get diffraction pattern! 
        - Downstream software analysis sees the diffraction pattern
        - Diffraction pattern needs to get converted into electron density
            - Using fourier transforms
            - So you get lots of peaks, where peaks = regions containing many electrons WOAHHHH  
        - Crystallographer then goes in and maps blobs to carbon, nitrogen, phenyl rings etc etc
        - Software iteratively adjusts atom positions checking for chemistry and physics to align in the structure

-	Cryo EM
    - Resolution RANGE: 2-4
        - 2.5A = X ray quality
        - 2 side chains are clear. H-bond networks are apparent
            
    - Note: Resolution isn't always uniform
        - Because maybe the core part of the crystal is 2.4A, but then the binding interface is 2.7A, and then on proteins, maybe a flexible tail = 5A, and flexible tail is not visible 

    - Method:
        - Electrons passing through individual frozen proteins.
        - Reconstructs a 3D electrostatic potential map (often informally called a density map).
        - Each protein: is frozen in vitreous (non-crystalline) ice,  has a random orientation, is an individual molecule.
        - Each protein produces a tiny, very noisy 2D projection.
        - Software is given millions of 2d images, needs to project into 3D:
            - Software groups similar views/ie: orientations, aligns them, and reconstructs the 3D structure.
    
-	NMR: Nuclear magnetic resonance
    - What it is
        - NMR doesn't “see” the molecule as an image. It measures how atomic nuclei behave in a magnetic field
        - Think of it as solving a huge 3D puzzle from many distance and angle measurements. 
        - Outputs:  set of distance, angle, and orientation constraints

    - How NMR works
        - Sample is dissolved in solution. 
        - Molecules tumble freely instead of being crystallized. 
        - Sample is placed inside a very strong magnetic field (typically 600–1200 MHz spectrometers). 
        - Atomic nuclei (commonly ¹H, ¹³C, ¹⁵N) behave like tiny magnets. 
        - Radiofrequency pulses perturb these nuclear spins. 
        - Nuclei emit tiny radio signals as they relax. 
        - Software Fourier transforms these time-domain signals into spectra. 
            - Just like X-ray Crystallography uses Fourier transforms on diffraction patterns, NMR uses Fourier transforms on FIDs (Free Induction Decays).
    - NMR structures usually contain multiple models, so if they all overlap well: good!!


### Force Fields:

- How to think about force fields?
    - What are we actually trying to solve?
        - Imagine you have an RNA molecule with 100 atoms.
        - Every atom is pulling or pushing on every other atom.
        
Newton's Law is F=ma

So if we know the force on every atom, we can compute

acceleration
velocity
new position

repeat 500 billion times

That's literally molecular dynamics.

So the real question becomes:

Where do the forces come from?
Quantum mechanics!!

Electrons determine chemistry.
Electrons determine
	bonds 
	charges 
	hydrogen bonds 
	steric repulsion 
	resonance 
	polarization 
The true equations are the Schrödinger equation.
Unfortunately...
For a protein of 50,000 atoms
the Schrödinger equation is essentially impossible to solve every femtosecond.
You'd need absurd amounts of computation.

Force fields are a giant approximation
Instead of solving quantum mechanics every timestep, people asked:
Can we just build equations that behave similarly?
That's exactly what a force field is.
Instead of
"simulate electrons"
it says
"pretend atoms are connected by springs and electrical charges."

SAME THING FOR VIDEO GAMES PHYSICS ENGINES!!

Suppose you're making Minecraft.
You could simulate
	atomic nuclei 
	electron clouds 
	chemical bonds 
or...
you could simply say
blocks can't overlap
and
gravity pulls downward

So every atom interaction becomes one equation:
	bonds stretching = spring = AKA Hooke's Law
how do we understand it- so in real life technically the atom on a carbon-carbon bond for a base pair COULD come apart. But its so rare, that in the simulation we just make it so that it could never come apart, just stretches
	bond angles: H2O likes 104.5 degrees, so the force field pushes it back if the angle of the atoms change. 
	Dihedral angles – RNA Imagine A-B-C-D
	Rotate D around the B-C bond.
	Some rotations are good and others are impossible, is there a force field that can capture likes and dislikes?
	Electrostatics
	Van der waals

Equations are mostly the same across all force fields. But like: 
Carbon-carbon bond
AMBER might use
1.526 Å
CHARMM might use
1.529 Å

Spring constant for 
WHY DO WE HAVE DIFFERENT FORCE FIELDS:
Suppose you optimize a force field so RNA helices match experiments beautifully.
You might accidentally make protein side chains too rigid.
Someone else tunes the parameters for proteins.
Another group focuses on membrane lipids.
Another on DNA.
Another on intrinsically disordered proteins.

Isn’t this so FUNNY

So all you computational biologist – what did you learn???

ALL MODELS ARE WRONG BUT SOME ARE USEFUL. SO WHENEVER WE DO ANYTHING, ASK YOURSELF WHAT THE SIMULATION IS ACCOMPLISHING. And experimental is king, always.

