# Focused version that only finds names for molecules that will be drawn
# -*- coding: utf-8 -*-

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
from collections import defaultdict
from IPython.display import display
import matplotlib.pyplot as plt
import requests
from time import sleep

# Initialize the results dictionary
results = defaultdict(lambda: {'smiles': set(), 'count': 0, 'cluster_ids': set()})

# Data import
# pipeline = import_file(['/home/dtramontina/Research/sims/23.pyrene/0.0.datas/pyrene_decompressed_ice_221_q.data',
                    #    '/home/dtramontina/Research/sims/23.pyrene/06.tracks/500ps10000K-dumps/data.track.final']) # Waiting results from serafin
# pipeline = import_file('/home/dtramontina/Research/sims/23.pyrene/0.0.datas/pyrene_decompressed_ice_221_q.data') # Big system test for initial frame
pipeline = import_file('/home/dtramontina/Research/sims/23.pyrene/0.0.datas/pyrene_decompressed_ice_q.data') # Big system test for initial frame
# pipeline = import_file('/home/dtramontina/Research/sims/23.pyrene/06.tracks/500ps19999K-dumps/data.track.final') # Little system test for 19999K

# Create bonds
mod = CreateBondsModifier()
mod.mode = CreateBondsModifier.Mode.Pairwise
mod.intra_molecule_only = True
mod.set_pairwise_cutoff('C', 'C', 1.65)
mod.set_pairwise_cutoff('H', 'C', 1.25)
mod.set_pairwise_cutoff('H', 'H', 0.95)
pipeline.modifiers.append(mod)

# Cluster analysis
mod = ClusterAnalysisModifier()
mod.neighbor_mode = ClusterAnalysisModifier.NeighborMode.CutoffRange
mod.cutoff = 1.5
mod.sort_by_size = True
mod.unwrap_particles = True
mod.cluster_coloring = True
pipeline.modifiers.append(mod)

def get_molecule_name(mol, smiles):
    """Get the best available name for a molecule"""
    # Try RDKit common name first (fastest)
    try:
        name = Chem.MolToCommonName(mol)
        if name: return name
    except: pass
    
    # Try CIR IUPAC name (requires internet)
    try:
        url = f"https://cactus.nci.nih.gov/chemical/structure/{smiles}/iupac_name"
        response = requests.get(url, timeout=2)
        if response.ok: return response.text.strip()
    except: pass
    
    # Fallback to SMILES if no name found
    return smiles

def process_clusters_to_molecules(data):
    """Convert OVITO clusters to RDKit molecules, handling valence issues."""
    molecules = []
    
    # Get cluster and bond information
    clusters = data.particles['Cluster']
    particle_types = data.particles['Particle Type']
    bonds = data.particles.bonds  # OVITO bonds object
    positions = data.particles.positions
    
    # Create atom symbol mapping (modify according to your particle types)
    type_to_symbol = {1: 'C', 2: 'H'}  # Example: type 1=Carbon, 2=Hydrogen
    
    # Group atoms by cluster
    cluster_atoms = defaultdict(list)
    for atom_idx, (cluster_id, atom_type) in enumerate(zip(clusters, particle_types)):
        cluster_atoms[cluster_id].append((atom_idx, type_to_symbol.get(atom_type, 'X')))
    
    # Process each cluster
    for cluster_id, atoms in cluster_atoms.items():
        if len(atoms) <= 1:  # Skip single-atom clusters
            continue
            
        mol = Chem.RWMol()
        atom_map = {}  # Maps OVITO atom indices to RDKit atom indices
        
        # Add atoms to the molecule
        for ovito_idx, symbol in atoms:
            atom = Chem.Atom(symbol)
            rdkit_idx = mol.AddAtom(atom)
            atom_map[ovito_idx] = rdkit_idx
        
        # Add bonds (if bonds exist in the data)
        if isinstance(bonds, Bonds):
            bond_array = bonds.topology  # Get bond pairs as Nx2 array
            for a_ovito, b_ovito in bond_array:
                # Only add bonds where both atoms are in this cluster
                if a_ovito in atom_map and b_ovito in atom_map:
                    try:
                        mol.AddBond(atom_map[a_ovito], atom_map[b_ovito], BondType.SINGLE)
                    except:
                        # Handle bond addition errors (rare)
                        continue
        
        # Try to sanitize the molecule
        try:
            # Fix common hydrogen valence issues
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == 'H' and atom.GetDegree() > 1:
                     # Keep only the first bond for this hydrogen
                      bonds = list(atom.GetBonds())
                      for bond in bonds[1:]:
                          mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            # Sanitize the molecule to check for valence issues
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
            # Handle invalid molecules (optional: try to fix before discarding)
            # print(f"Cluster {cluster_id} failed sanitization: {str(e)}")
            continue
    return molecules
def generate_mol_report(cluster_molecules):
    """
    Generates a molecular report from processed cluster data.
    
    Args:
        cluster_molecules: Output from process_clusters_to_molecules(data)
        
    Returns:
        List of dictionaries sorted by abundance, each containing:
        - 'count': Number of occurrences
        - 'cluster_ids': List of cluster IDs
        - 'representative_mol': RDKit molecule object
        - 'smiles': Canonical SMILES string
        - 'atom_count': Number of atoms
        - 'name': Common name of molecule
    """
    # Group identical molecules by SMILES
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
        mol_groups[smiles]['atom_count'] = mol_data['atom_count']  # All should be same
    
    # Generate the report sorted by abundance
    mol_report = []
    for smiles, group_data in sorted(mol_groups.items(),
                                   key=lambda x: x[1]['count'],
                                   reverse=True):
        # Use first molecule as representative
        representative = group_data['molecules'][0]
        
        mol_report.append({
            'count': group_data['count'],
            'cluster_ids': sorted(group_data['cluster_ids']),
            'representative_mol': representative['mol'],
            'smiles': smiles,
            'atom_count': group_data['atom_count'],
            'name': get_molecule_name(representative['mol'], smiles)
        })
    
    return mol_report

for frame in range(pipeline.source.num_frames):
    data = pipeline.compute(frame)
    bonds = data.particles.bonds
    clusters = data.tables['clusters']
    particle_types = data.particles['Particle Type']
    positions = data.particles.positions
    molecules = process_clusters_to_molecules(data)
    mol_report = generate_mol_report(molecules)
    for i, entry in enumerate(report[:10], 1):  # Print top 10
        print(f"{i}. {entry['name']} (Count: {entry['count']})")
        print(f"   SMILES: {entry['smiles']}")
        print(f"   Cluster IDs: {entry['cluster_ids'][:5]}...")  # Show first 5 IDs
    # Draw molecules
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(300, 300), kekulize=True)
                plt.figure(figsize=(3,3))
                plt.imshow(img)
                plt.axis('off')
                plt.title(f'{entry["name"]} ({entry["count"]})')
                display(plt.gcf())
                # plt.savefig(f"top{rank}_{formula}_mol{i}.png", 
                #         bbox_inches='tight', dpi=150) # Uncomment to save images
                plt.close()
                # print(f"    Saved: top{rank}_{formula}_mol{i}.png")
        except Exception as e:
            print(f"    Could not draw molecule {i}: {e}")
            continue
        
