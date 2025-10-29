import os
import pytest
from molfinder.processor import LammpsProcessor

@pytest.fixture
def water_data_path():
    """Fixture to provide the path to the water.data file."""
    return os.path.join(os.path.dirname(__file__), 'data', 'water.data')

def test_process_water(water_data_path):
    """
    Tests processing of a simple water molecule.
    """
    # The atom mapping is not needed here because the data file has comments
    processor = LammpsProcessor(water_data_path, atom_mapping_str="")
    result = processor.process()

    # There should be one type of molecule found: H2O
    assert len(result) == 1
    assert 'H2O' in result

    # There should be one molecule of H2O
    h2o_group = result['H2O']
    assert h2o_group['count'] == 1

    # The molecule should have the correct formula and number of atoms
    molecules = h2o_group['molecules']
    assert len(molecules) == 1
    molecule = molecules[0]
    assert molecule.formula == 'H2O'
    assert molecule.mol.GetNumAtoms() == 3