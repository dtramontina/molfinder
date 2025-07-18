import requests
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors

class Molecule:
    """
    A class to represent a molecule found in the simulation.
    It holds the RDKit molecule object, its SMILES string, name, and stoichiometry.
    """
    def __init__(self, rdkit_mol):
        if not isinstance(rdkit_mol, Chem.Mol):
            raise TypeError("rdkit_mol must be an RDKit Mol object.")
            
        self.mol = rdkit_mol
        self.smiles = Chem.MolToSmiles(self.mol, canonical=True)
        self.formula = Descriptors.rdMolDescriptors.CalcMolFormula(self.mol)
        self.name = self._get_molecule_name()

    def _get_molecule_name(self):
        """
        Tries to fetch a common name for the molecule.
        Falls back to the SMILES string if no name is found.
        """
        # Try RDKit common name (fast, local)
        try:
            name = Chem.rdMolDescriptors.CalcMolFormula(self.mol)
            if name: return name
        except:
            pass

        # Try IUPAC name from NIH Cactus service (slower, requires internet)
        try:
            url = f"https://cactus.nci.nih.gov/chemical/structure/{self.smiles}/iupac_name"
            response = requests.get(url, timeout=5)
            if response.ok and response.text.strip():
                return response.text.strip().capitalize()
        except requests.exceptions.RequestException:
            # If network fails, just proceed to fallback
            pass
        
        # Fallback to SMILES
        return self.smiles

    def draw_to_file(self, filename, size=(300, 300)):
        """Generates a 2D image of the molecule and saves it to a file."""
        Draw.MolToFile(self.mol, filename, size=size, kekulize=True)

    def __repr__(self):
        return f"<Molecule name='{self.name}' smiles='{self.smiles}'>"