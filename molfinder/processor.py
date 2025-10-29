import re
from collections import defaultdict
import ovito
from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier
from ovito.data import ParticleType
from rdkit import Chem
from .molecule import Molecule

class LammpsProcessor:
    """
    Processes a LAMMPS file to identify and classify molecules.
    """
    def __init__(self, filepath, atom_mapping_str):
        self.filepath = filepath
        self.user_atom_mapping_str = atom_mapping_str

    def _ensure_element_names(self, data):
        """
        Assigns element names to particle types if they are missing.
        """
        type_property = data.particles_.particle_types_

        # Check if names are already present (e.g., from a dump file with 'element' column)
        if all(pt.name for pt in type_property.types):
            return

        # Strategy 1: Parse from LAMMPS data file 'Masses' section
        try:
            with open(self.filepath, 'r') as f:
                content = f.read()
            masses_section = re.search(r'^\s*Masses\s*\n\n(.*?)\n\n', content, re.DOTALL | re.MULTILINE)
            if masses_section:
                lines = masses_section.group(1).strip().split('\n')
                for line in lines:
                    match = re.search(r'^\s*(\d+)\s+[\d.]+\s*#\s*([a-zA-Z]+)', line.strip())
                    if match:
                        type_id = int(match.group(1))
                        symbol = match.group(2).capitalize()
                        pt = type_property.type_by_id_(type_id)
                        if pt and not pt.name:
                            pt.name_ = symbol
                if all(pt.name for pt in type_property.types): return
        except Exception:
            pass

        # Strategy 2: Use user-provided comma-separated string
        if self.user_atom_mapping_str:
            symbols = [s.strip().capitalize() for s in self.user_atom_mapping_str.split(',')]
            for i, symbol in enumerate(symbols):
                type_id = i + 1
                pt = type_property.type_by_id_(type_id)
                if pt and not pt.name:
                    pt.name_ = symbol
            return

        raise ValueError("Could not determine element names for all particle types.")

    def _create_rdkit_mol(self, atom_indices, data):
        """Creates an RDKit molecule from a set of atom indices."""
        mol = Chem.RWMol()
        atom_map = {}  # Maps original particle index to RDKit atom index

        particle_types = data.particles.particle_types
        positions = data.particles.positions
        bonds = data.particles.bonds

        for idx in atom_indices:
            symbol = particle_types.type_by_id(data.particles['Particle Type'][idx]).name
            rdkit_idx = mol.AddAtom(Chem.Atom(symbol))
            atom_map[idx] = rdkit_idx
        
        for a, b in bonds.topology:
            if a in atom_map and b in atom_map:
                mol.AddBond(atom_map[a], atom_map[b], Chem.BondType.SINGLE)
        
        try:
            Chem.SanitizeMol(mol)
            return mol
        except Exception:
            return mol.GetMol()

    def process(self):
        """
        Executes the full analysis pipeline.
        """
        pipeline = import_file(self.filepath)
        data = pipeline.compute()

        # Step 1: Ensure particle types have element names.
        self._ensure_element_names(data)

        # Step 2: Create bonds and identify molecules.
        bonds_modifier = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise, intra_molecule_only=True)
        defaults = ParticleType.load_defaults()
        type_list = list(data.particles.particle_types.types)
        for i in range(len(type_list)):
            for j in range(i, len(type_list)):
                type_a = type_list[i]
                type_b = type_list[j]
                radius_a = defaults.get(type_a.name, {}).get('covalent_radius', 0.8)
                radius_b = defaults.get(type_b.name, {}).get('covalent_radius', 0.8)
                bonds_modifier.set_pairwise_cutoff(type_a.name, type_b.name, radius_a + radius_b)

        pipeline.modifiers.append(bonds_modifier)
        data = pipeline.compute() # Re-compute to get molecule identifiers.

        # Step 3: Group atoms by the new 'Molecule' identifier.
        molecule_atoms = defaultdict(list)
        if 'Molecule' in data.particles:
            for i, mol_id in enumerate(data.particles['Molecule']):
                if mol_id > 0: # 0 means not part of a molecule
                    molecule_atoms[mol_id].append(i)

        # Step 4: Convert groups to RDKit molecules and generate report.
        molecule_groups = defaultdict(lambda: {'count': 0, 'molecules': []})
        for mol_id, atom_indices in molecule_atoms.items():
            if len(atom_indices) > 1:
                rdkit_mol = self._create_rdkit_mol(atom_indices, data)
                if rdkit_mol:
                    mol_obj = Molecule(rdkit_mol)
                    molecule_groups[mol_obj.formula]['count'] += 1
                    if not any(m.smiles == mol_obj.smiles for m in molecule_groups[mol_obj.formula]['molecules']):
                        molecule_groups[mol_obj.formula]['molecules'].append(mol_obj)

        return dict(sorted(molecule_groups.items(), key=lambda item: item[1]['count'], reverse=True))
