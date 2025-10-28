import re
import io
import logging
from collections import defaultdict
from ovito import *
from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.vis import *
from ovito.data import *

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import MolToSmiles
from rdkit.Chem import BondType
from rdkit.Chem import PyMol
from .molecule import Molecule

# from IPython.display import display
# import matplotlib.pyplot as plt
# import requests
# from time import sleep





class LammpsProcessor:
    """
    Processes a LAMMPS file to identify and classify molecules.
    """
    def __init__(self, filepath, atom_mapping_str):
        """
        Initializes the processor.
        
        Args:
            filepath (str): Path to the LAMMPS data or dump file.
            atom_mapping_str (str): Comma-separated string like "C,H,O".
        """
        self.filepath = filepath
        self.user_atom_mapping_str = atom_mapping_str
        self.atom_type_map = {}
        self.log_stream = io.StringIO()

    def _setup_logging(self):
        """Sets up OVITO's logging to capture messages."""
        ovito.enable_logging(self.log_stream, level=logging.INFO)

    def _setup_particle_types(self, data):
        """
        Determines and assigns particle types based on file content or user input.
        """
        particle_types = data.particles_.create_property(ParticleType)

        # Strategy 1: Attempt to parse from LAMMPS data file 'Masses' section comments
        try:
            with open(self.filepath, 'r') as f:
                content = f.read()

            masses_section = re.search(r'^\s*Masses\s*\n\n(.*?)\n\n', content, re.DOTALL | re.MULTILINE)
            if masses_section:
                lines = masses_section.group(1).strip().split('\n')
                found_mapping = {}
                for line in lines:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    match = re.search(r'^\s*(\d+)\s+[\d.]+\s*#\s*([a-zA-Z]+)', line)
                    if match:
                        type_id = int(match.group(1))
                        symbol = match.group(2).capitalize()
                        found_mapping[type_id] = symbol

                if len(found_mapping) > 0:
                    self.atom_type_map = found_mapping
                    for type_id, symbol in self.atom_type_map.items():
                        pt = particle_types.add_type_name(symbol)
                        pt.mass = ParticleType.load_defaults()[symbol].mass
                    return # Successfully mapped from comments
        except Exception:
            pass # File might not be a data file or might be malformed

        # Strategy 2: Fallback to user-provided string (e.g., "C,H,O")
        if self.user_atom_mapping_str:
            symbols = [s.strip().capitalize() for s in self.user_atom_mapping_str.split(',')]
            for i, symbol in enumerate(symbols):
                type_id = i + 1
                self.atom_type_map[type_id] = symbol
                pt = particle_types.add_type_name(symbol)
                pt.mass = ParticleType.load_defaults()[symbol].mass
            return

        # Strategy 3: If it's a dump file with named types, infer from there
        if 'Particle Type' in data.particles and data.particles['Particle Type'].is_string_property:
            type_names = sorted(list(set(data.particles['Particle Type'])))
            for i, name in enumerate(type_names):
                type_id = i + 1
                self.atom_type_map[type_id] = name
                pt = particle_types.add_type_name(name)
                pt.mass = ParticleType.load_defaults().get(name, {}).get('mass', 1.0)
            return

        raise ValueError("Could not determine atom type mapping. Please provide a mapping string.")

    def _create_rdkit_mol(self, cluster_data):
        """Creates an RDKit molecule from a cluster of atoms and bonds."""
        mol = Chem.RWMol()
        atom_map = {}  # Maps OVITO atom index to RDKit atom index

        for ovito_idx, symbol in cluster_data['atoms']:
            rdkit_idx = mol.AddAtom(Chem.Atom(symbol))
            atom_map[ovito_idx] = rdkit_idx
        
        for a_ovito, b_ovito in cluster_data['bonds']:
            if a_ovito in atom_map and b_ovito in atom_map:
                mol.AddBond(atom_map[a_ovito], atom_map[b_ovito], Chem.BondType.SINGLE)
        
        try:
            Chem.SanitizeMol(mol)
            return mol
        except Exception:
            return mol.GetMol() # Return unsanitized if sanitization fails

    def process(self):
        """
        Executes the full analysis pipeline on the LAMMPS file.
        """
        self._setup_logging()
        pipeline = import_file(self.filepath)

        # The setup of particle types is now deferred until the pipeline is evaluated
        def setup_types_and_bonds(frame, data):
            if not self.atom_type_map: # Only run once
                self._setup_particle_types(data)

            # Now that types are set, create bonds
            bonds_modifier = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)

            # Dynamically set cutoffs from OVITO's default database
            defaults = ParticleType.load_defaults()
            type_list = list(data.particles.particle_types.types)
            for i in range(len(type_list)):
                for j in range(i, len(type_list)):
                    type_a = type_list[i]
                    type_b = type_list[j]

                    # A simple cutoff heuristic: sum of covalent radii
                    radius_a = defaults.get(type_a.name, {}).get('covalent_radius', 0.8)
                    radius_b = defaults.get(type_b.name, {}).get('covalent_radius', 0.8)
                    cutoff = radius_a + radius_b

                    bonds_modifier.set_pairwise_cutoff(type_a.name, type_b.name, cutoff)

            pipeline.modifiers.append(bonds_modifier)

        pipeline.modifiers.append(setup_types_and_bonds)

        pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=0.1, sort_by_size=True))

        data = pipeline.compute()

        clusters_table = data.tables['clusters']
        bonds = data.particles.bonds

        molecule_groups = defaultdict(lambda: {'count': 0, 'molecules': []})

        type_id_to_symbol_map = {i: ptype.name for i, ptype in data.particles.particle_types.types_by_id.items()}

        for cluster_index in range(len(clusters_table.clusters)):
            particle_indices = clusters_table.particle_indices_in_cluster(cluster_index)
            if len(particle_indices) <= 1: continue

            cluster_bonds = [(a, b) for (a, b) in bonds.topology if a in particle_indices and b in particle_indices]
            cluster_atoms = [(idx, type_id_to_symbol_map[data.particles['Particle Type'][idx]]) for idx in particle_indices]
            
            rdkit_mol = self._create_rdkit_mol({'atoms': cluster_atoms, 'bonds': cluster_bonds})
            if rdkit_mol:
                mol_obj = Molecule(rdkit_mol)
                molecule_groups[mol_obj.formula]['count'] += 1

                if not any(m.smiles == mol_obj.smiles for m in molecule_groups[mol_obj.formula]['molecules']):
                    molecule_groups[mol_obj.formula]['molecules'].append(mol_obj)

        return dict(sorted(molecule_groups.items(), key=lambda item: item[1]['count'], reverse=True))

    def get_log_messages(self):
        """Returns the captured log messages."""
        return self.log_stream.getvalue()
