"""Run inside ChimeraX with redocked complex #1 and crystal complex #2.
Uses the verified heavy-atom mapping for the two supplied PDB files.
Fits receptors only (all corresponding CA atoms), then measures ligand RMSD.
Does not rename atoms, change internal coordinates, or overwrite input files.
"""
import numpy as np
from chimerax.core.commands import run
from chimerax.core.errors import UserError

dock_names = 'C1 C2 C3 C4 C5 C6 C7 C8 N1 C9 C10 O1 C11 C12 N2 C13 C14 C15 N3 N4 O2 C16 C17 N5 C18 N6 C19 C20 C21 C22 F1 C23 F2 C24 C25 O3 C26 C27 O4'.split()
ref_names = 'C4 C14 C6 C11 C7 C10 C12 C13 N1 C3 C37 O36 C31 C32 N33 C34 C35 C30 N29 N5 O38 C39 C2 N8 C15 N9 C16 C17 C18 C19 F28 C22 F23 C21 C20 O24 C25 C26 O27'.split()
models = {m.id: m for m in session.models.list() if hasattr(m, 'atoms')}
if (1,) not in models or (2,) not in models:
    raise UserError('Expected redocking model #1 and crystal model #2.')
dock, ref = models[(1,)], models[(2,)]
def ligand_atoms(model, number, names, residue_name):
    residues = [r for r in model.residues if r.number == number and r.name == residue_name]
    if len(residues) != 1:
        raise UserError('Expected one %s %s in model #%s' % (residue_name, number, model.id_string))
    heavy = [a for a in residues[0].atoms if a.element.number > 1]
    lookup = {a.name: a for a in heavy}
    if len(heavy) != 39 or len(lookup) != 39 or set(lookup) != set(names):
        raise UserError('Ligand atoms differ from the verified input PDB. Stopping.')
    return [lookup[n] for n in names]
a = ligand_atoms(dock, 900, dock_names, 'UNK')
b = ligand_atoms(ref, 1101, ref_names, 'CVQ')
if any(x.element.number != y.element.number for x,y in zip(a,b)):
    raise UserError('Element mismatch in ligand atom mapping.')
run(session, 'align #1/A@CA toAtoms #2/A@CA matchNumbering true')
x = np.array([dock.scene_position * atom.coord for atom in a], dtype=float)
y = np.array([ref.scene_position * atom.coord for atom in b], dtype=float)
rmsd = float(np.sqrt(np.mean(np.sum((x-y)**2, axis=1))))
session.logger.info('CVQ LIGAND RMSD: %.6f angstroms; 39 heavy-atom pairs; receptor CA alignment; NO separate ligand fit.' % rmsd)
