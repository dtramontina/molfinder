from collections import defaultdict
from ovito.io import import_file
from ovito.modifiers import ClusterAnalysisModifier, CreateBondsModifier
from ovito.data import Bonds
from rdkit import Chem
from .molecule import Molecule

class LammpsProcessor:
    """
    Processes a LAMMPS file to identify and classify molecules.
    """
    def __init__(self, filepath, atom_mapping_str):
        """
        Initializes the processor.
        
        Args:
            filepath (str): Path to the LAMMPS data or dump file.
            atom_mapping_str (str): Comma-separated string like "1:C, 2:H".
        """
        self.filepath = filepath
        self.atom_type_map = self._parse_atom_mapping(atom_mapping_str)
        self.symbol_map = {v: k for k, v in self.atom_type_map.items()}

    def _parse_atom_mapping(self, mapping_str):
        """Parses the atom mapping string into a dictionary."""
        try:
            atom_map = dict(item.split(':') for item in mapping_str.split(','))
            return {int(k.strip()): v.strip().capitalize() for k, v in atom_map.items()}
        except (ValueError, IndexError):
            raise ValueError("Invalid atom type mapping format. Use '1:C, 2:H, ...'")

    def _create_rdkit_mol(self, cluster_data):
        """Creates an RDKit molecule from a cluster of atoms and bonds."""
        mol = Chem.RWMol()
        atom_map = {}  # Maps OVITO atom index to RDKit atom index

        # Add atoms
        for ovito_idx, symbol in cluster_data['atoms']:
            rdkit_idx = mol.AddAtom(Chem.Atom(symbol))
            atom_map[ovito_idx] = rdkit_idx
        
        # Add bonds
        for a_ovito, b_ovito in cluster_data['bonds']:
            if a_ovito in atom_map and b_ovito in atom_map:
                mol.AddBond(atom_map[a_ovito], atom_map[b_ovito], Chem.BondType.SINGLE)
        
        # Sanitize and finalize the molecule
        try:
            # Basic sanitization, which also calculates valency, aromaticity, etc.
            Chem.SanitizeMol(mol)
            return mol
        except Chem.rdchem.AtomValenceException:
            # A more robust solution could try to fix valency issues here,
            # for now, we just return the unsanitized molecule.
            return mol.GetMol()
        except Exception:
            return None

    def process(self):
        """
        Executes the full analysis pipeline on the LAMMPS file.
        Returns a dictionary of molecules grouped by stoichiometry.
        """
        pipeline = import_file(self.filepath)

        # Set atom types from mapping to enable OVITO's element-based modifiers
        def assign_particle_types(frame, data):
            particle_types = data.particles_.create_property('Particle Type', data=data.particles['Particle Type'])
            for type_id, symbol in self.atom_type_map.items():
                particle_types.type_by_id(type_id).name = symbol
        pipeline.modifiers.append(assign_particle_types)

        # Create bonds based on typical covalent radii cutoffs.
        # This is a general approach; specific systems might need fine-tuning.
        bonds_modifier = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
        bonds_modifier.set_pairwise_cutoff('C', 'C', 1.7)
        bonds_modifier.set_pairwise_cutoff('H', 'C', 1.3)
        bonds_modifier.set_pairwise_cutoff('H', 'H', 1.0)
        bonds_modifier.set_pairwise_cutoff('O', 'C', 1.6)
        bonds_modifier.set_pairwise_cutoff('O', 'H', 1.1)
        bonds_modifier.set_pairwise_cutoff('N', 'C', 1.6)
        bonds_modifier.set_pairwise_cutoff('N', 'H', 1.2)
        pipeline.modifiers.append(bonds_modifier)

        # Identify connected clusters of atoms
        pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=0.1, sort_by_size=True))

        data = pipeline.compute()

        # Extract cluster information
        clusters_table = data.tables['clusters']
        bonds = data.particles.bonds

        molecule_groups = defaultdict(lambda: {'count': 0, 'molecules': []})

        for cluster_index in range(len(clusters_table.clusters)):
            particle_indices = clusters_table.particle_indices_in_cluster(cluster_index)
            if len(particle_indices) <= 1: continue # Skip single atoms

            cluster_bonds = [(a, b) for (a, b) in bonds.topology if a in particle_indices and b in particle_indices]
            cluster_atoms = [(idx, self.atom_type_map[data.particles['Particle Type'][idx]]) for idx in particle_indices]
            
            rdkit_mol = self._create_rdkit_mol({'atoms': cluster_atoms, 'bonds': cluster_bonds})
            if rdkit_mol:
                mol_obj = Molecule(rdkit_mol)
                molecule_groups[mol_obj.formula]['count'] += 1
                
                # Add if it's a new isomer for this formula
                if not any(m.smiles == mol_obj.smiles for m in molecule_groups[mol_obj.formula]['molecules']):
                    molecule_groups[mol_obj.formula]['molecules'].append(mol_obj)

        # Sort by abundance
        return dict(sorted(molecule_groups.items(), key=lambda item: item[1]['count'], reverse=True))


