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

You can also open the pdb file

### Step 3 Prep System

#### Step 3.1 - Get force field files

The standard amber99sb-ildn is primarily a protein force field. The RNA improvements (parmbsc0, parmbsc1, OL3, etc.) are distributed as modified AMBER force fields that must already exist in your GROMACS share/gromacs/top directory or be installed separately.

First we check if you have Amber14SB + OL3/parmbsc1 installed

```bash 
gmx pdb2gmx \
    -f 1LDZ_model1.pdb \
    -o processed.gro \
    -p topol.top \
    -i posre.itp
```
Choose
Force field:
Amber14SB + OL3/parmbsc1

Water:
TIP3P
If your installation only shows
amber99sb-ildn
I would not use it for RNA.
Instead install an RNA-capable force field first.

Verify it is installed:

Check your force field files:
gmx pdb2gmx -h

Step 2. Look inside the force field directory
Find where GROMACS stores force fields.
On Linux/WSL this is commonly:
gmx --version
Then list the force fields:
ls /usr/local/gromacs/share/gromacs/top


find /usr/local/gromacs/share/gromacs/top -name "rna.rtp" -o -name "nucleic.rtp" -o -name "*ol3*" -o -name "*bsc1*"
 


RG	Internal guanosine	Middle of the chain
RG3	3′-terminal guanosine	Last nucleotide (free 3′-OH)
RGN	5′-terminal guanosine	First nucleotide (free 5′ end)

Install it if you are missing it:



Download the official Amber force field conversion from the MacKerell/Lemkul GROMACS force field repository or the AMBER/GROMACS community distribution.
These provide directories like
amber14sb.ff
containing
rna.rtp
dna.rtp
aminoacids.rtp
ffbonded.itp
ffnonbonded.itp
with the OL3/parmbsc1 corrections already incorporated.


For RNA in GROMACS today, I would recommend one of these, in order of preference:
Force field	Recommendation
AMBER OL3 + parmbsc1	⭐ Best-established choice for RNA aptamers and folded RNA
AMBER14SB + OL3	⭐ Excellent if available in your installation
CHARMM36m nucleic acid	Very good alternative
AMBER99SB-ILDN	Suitable for proteins, but not recommended by itself for RNA because it lacks the later RNA backbone corrections (OL3/parmbsc1).


This is called getting topology and will produce .top file. expected output files:
•	processed.gro 
•	topol.top 
•	posre.itp


gmx pdb2gmx -f clean.pdb -o processed.gro -p topol.top -i posre.itp


Select force fields Select your modified AMBER force field and choose the TIP3P explicit water
•	Use the AMBER99SB-ILDN force field or AMBER14sb with parmbsc1
•	These are better for nucleic acids than default 
### Simulation: prep system: => boxed.gro
gmx editconf \
    -f processed.gro \
    -o boxed.gro \
    -c \
    -d 1.0 \
    -bt dodecahedron
Typical:
•	1.0 nm padding 
•	dodecahedron box


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


