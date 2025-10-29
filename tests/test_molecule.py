import pytest
from rdkit import Chem
from molfinder.molecule import Molecule

@pytest.fixture
def methane_mol():
    """Fixture to provide an RDKit molecule object for methane."""
    return Chem.MolFromSmiles('C')

def test_molecule_initialization(methane_mol):
    """Tests the initialization of the Molecule class."""
    molecule = Molecule(methane_mol)
    assert molecule.mol is not None
    assert molecule.smiles == 'C'
    assert molecule.formula == 'CH4'
    assert molecule.name == 'CH4'

def test_molecule_invalid_input():
    """Tests that the Molecule class raises a TypeError for invalid input."""
    with pytest.raises(TypeError):
        Molecule("not a molecule")