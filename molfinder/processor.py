from collections import defaultdict
from ovito.io import import_file
from ovito.modifiers import ClusterAnalysisModifier, CreateBondsModifier
from ovito.data import Bonds
from rdkit import Chem
from rdkit.Chem import BondType
from .molecule import Molecule
import requests

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

    def get_molecule_name(self, mol, smiles):
        """Get the best available name for a molecule"""
        try:
            name = Chem.MolToCommonName(mol)
            if name: return name
        except: pass

        try:
            url = f"https://cactus.nci.nih.gov/chemical/structure/{smiles}/iupac_name"
            response = requests.get(url, timeout=2)
            if response.ok: return response.text.strip()
        except: pass

        return smiles

    def process_clusters_to_molecules(self, data):
        """Convert OVITO clusters to RDKit molecules, handling valence issues."""
        molecules = []
        
        clusters = data.particles['Cluster']
        particle_types = data.particles['Particle Type']
        bonds = data.particles.bonds
        positions = data.particles.positions
        
        type_to_symbol = self.atom_type_map

        cluster_atoms = defaultdict(list)
        for atom_idx, (cluster_id, atom_type) in enumerate(zip(clusters, particle_types)):
            cluster_atoms[cluster_id].append((atom_idx, type_to_symbol.get(atom_type, 'X')))

        for cluster_id, atoms in cluster_atoms.items():
            if len(atoms) <= 1:
                continue

            mol = Chem.RWMol()
            atom_map = {}

            for ovito_idx, symbol in atoms:
                atom = Chem.Atom(symbol)
                rdkit_idx = mol.AddAtom(atom)
                atom_map[ovito_idx] = rdkit_idx

            if isinstance(bonds, Bonds):
                bond_array = bonds.topology
                for a_ovito, b_ovito in bond_array:
                    if a_ovito in atom_map and b_ovito in atom_map:
                        try:
                            mol.AddBond(atom_map[a_ovito], atom_map[b_ovito], BondType.SINGLE)
                        except:
                            continue

            try:
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() == 'H' and atom.GetDegree() > 1:
                          bonds_list = list(atom.GetBonds())
                          for bond in bonds_list[1:]:
                              mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                Chem.SanitizeMol(mol)
                smiles = Chem.MolToSmiles(mol, canonical=True)
                molecules.append({
                    'cluster_id': cluster_id,
                    'mol': mol,
                    'smiles': smiles,
                    'atom_count': len(atoms),
                    'positions': [positions[idx] for idx, _ in atoms]
                })
            except Exception as e:
                continue
        return molecules

    def generate_mol_report(self, cluster_molecules):
        """
        Generates a molecular report from processed cluster data.
        """
        mol_groups = defaultdict(lambda: {
            'count': 0,
            'cluster_ids': [],
            'molecules': [],
            'atom_count': None
        })

        for mol_data in cluster_molecules:
            smiles = mol_data['smiles']
            mol_groups[smiles]['count'] += 1
            mol_groups[smiles]['cluster_ids'].append(mol_data['cluster_id'])
            mol_groups[smiles]['molecules'].append(mol_data)
            mol_groups[smiles]['atom_count'] = mol_data['atom_count']

        mol_report = []
        for smiles, group_data in sorted(mol_groups.items(),
                                       key=lambda x: x[1]['count'],
                                       reverse=True):
            representative = group_data['molecules'][0]

            mol_report.append({
                'count': group_data['count'],
                'cluster_ids': sorted(group_data['cluster_ids']),
                'representative_mol': representative['mol'],
                'smiles': smiles,
                'atom_count': group_data['atom_count'],
                'name': self.get_molecule_name(representative['mol'], smiles)
            })

        return mol_report

    def process(self):
        """
        Executes the full analysis pipeline on the LAMMPS file.
        Returns a dictionary of molecules grouped by stoichiometry.
        """
        pipeline = import_file(self.filepath)

        def assign_particle_types(frame, data):
            particle_types = data.particles_.create_property('Particle Type', data=data.particles['Particle Type'])
            for type_id, symbol in self.atom_type_map.items():
                particle_types.type_by_id(type_id).name = symbol
        pipeline.modifiers.append(assign_particle_types)

        bonds_modifier = CreateBondsModifier()
        bonds_modifier.mode = CreateBondsModifier.Mode.Pairwise
        bonds_modifier.intra_molecule_only = True
        bonds_modifier.set_pairwise_cutoff('C', 'C', 1.65)
        bonds_modifier.set_pairwise_cutoff('H', 'C', 1.25)
        bonds_modifier.set_pairwise_cutoff('H', 'H', 0.95)
        pipeline.modifiers.append(bonds_modifier)

        cluster_modifier = ClusterAnalysisModifier()
        cluster_modifier.neighbor_mode = ClusterAnalysisModifier.NeighborMode.CutoffRange
        cluster_modifier.cutoff = 1.5
        cluster_modifier.sort_by_size = True
        cluster_modifier.unwrap_particles = True
        cluster_modifier.cluster_coloring = True
        pipeline.modifiers.append(cluster_modifier)

        data = pipeline.compute()

        molecules = self.process_clusters_to_molecules(data)
        mol_report = self.generate_mol_report(molecules)

        molecule_groups = defaultdict(lambda: {'count': 0, 'molecules': []})
        for entry in mol_report:
            mol_obj = Molecule(entry['representative_mol'])
            mol_obj.name = entry['name']
            molecule_groups[mol_obj.formula]['count'] += entry['count']
            
            is_new_isomer = not any(m.smiles == mol_obj.smiles for m in molecule_groups[mol_obj.formula]['molecules'])
            if is_new_isomer:
                molecule_groups[mol_obj.formula]['molecules'].append(mol_obj)

        return dict(sorted(molecule_groups.items(), key=lambda item: item[1]['count'], reverse=True))