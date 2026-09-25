import MDAnalysis as mda, numpy as np

def shape(pos, m):
    x = pos - np.average(pos, axis=0, weights=m)
    T = (x * m[:, None]).T @ x / m.sum()          # gyration tensor
    l = np.linalg.eigvalsh(T)
    rg = np.sqrt(l.sum())
    kappa2 = 1 - 3*(l[0]*l[1] + l[1]*l[2] + l[2]*l[0]) / l.sum()**2  # 0 = sphere, 1 = rod
    return rg, kappa2

u = mda.Universe("design_A_ensemble.pdb")       # multi-model PDB
atoms = u.select_atoms("all")
out = np.array([shape(atoms.positions, atoms.masses) for ts in u.trajectory])
print("Rg %.1f ± %.1f Å | kappa2 %.2f ± %.2f" % (*out.mean(0)[:1], out[:, 0].std(),
      out.mean(0)[1], out[:, 1].std()))