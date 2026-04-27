from openmm.app import PDBFile
from pdbfixer import PDBFixer
import mdtraj as md

def preparar_estructura(input_pdb, output_cleaned, output_ca):
    fixer = PDBFixer(input_pdb)
    fixer.removeHeterogens(keepWater=False)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_cleaned, 'w'))
    pdb = md.load(output_cleaned)
    indices_ca = [atom.index for atom in pdb.topology.atoms if atom.name == "CA"]
    solo_ca = pdb.atom_slice(indices_ca)
    solo_ca.save(output_ca)
    print(f"{input_pdb}: {len(indices_ca)} carbonos alfa → {output_ca}")

preparar_estructura("1z98_assembly.pdb", "protein_abierto.pdb", "ca_only_abierto.pdb")
preparar_estructura("2b5f.pdb",          "protein_cerrado.pdb", "ca_only_cerrado.pdb")
