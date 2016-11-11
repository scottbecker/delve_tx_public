import random
import string
import unittest
from transcriptic_tools import utils
from autoprotocol import Unit
from transcriptic_tools.utils import ul, ml
from autoprotocol import Container
from autoprotocol.container_type import _CONTAINER_TYPES
from uuid import uuid4

def create_blank_plate(container_shortname):
    plate = Container(str(uuid4), _CONTAINER_TYPES[container_shortname], name='test plate', 
                      storage='cold_4',cover=None)   
    
    for well in plate.all_wells():
        well.volume = ul(0)
        
    return plate
        
def generate_random_transriptic_id():
    return ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(20))

class TestUtils(unittest.TestCase):
    def test_total_plate_volume(self):
        
        plate = create_blank_plate('24-deep')
        
        self.assertEqual(utils.total_plate_volume(plate),ul(0))
        
        well = plate.well(0)
        
        well.volume = ul(100)
        
        self.assertEqual(utils.total_plate_volume(plate),ul(100))
        
        
    def test_available_plate_volume(self):
        
        plate = create_blank_plate('24-deep')
                
        self.assertEqual(utils.space_available(plate),ml(10)*24)   
        
        
    def test_floor_volume(self):
        self.assertEqual(utils.floor_volume(ul(12.8)),ul(12))
        
        self.assertEqual(utils.floor_volume(ul(12.2)),ul(12))
        
        self.assertEqual(utils.floor_volume(ul(12)),ul(12))
        
        
    def test_total_plate_volume(self):
        
        plate = create_blank_plate('24-deep')
        
        plate.well(0).volume = ml(10)
        plate.well(1).volume = ml(5) + ul(30)

        self.assertEqual(utils.total_plate_volume(plate,aspiratable=True),ml(15))  
        
        self.assertEqual(utils.total_plate_volume(plate,aspiratable=False),ml(15)+ul(30)) 
        
    def test_get_volume(self):
        
        plate = create_blank_plate('24-deep')
                
        plate.well(0).volume = ml(10)
        plate.well(1).volume = ml(5) + ul(30)

        self.assertEqual(utils.get_volume(plate,aspiratable=True),ml(15))  
        
        self.assertEqual(utils.get_volume(plate,aspiratable=False),ml(15)+ul(30))  
        
        self.assertEqual(utils.get_volume(plate.well(0),aspiratable=True),ml(10)-ul(15))  
        
        
    
    def test_convert_to_wellgroup(self):
        plate1 = create_blank_plate('6-flat')
        tube1 = create_blank_plate('micro-2.0')
        wellgroup = utils.convert_to_wellgroup([plate1,tube1])
        self.assertEqual(len(wellgroup),7)
        
        wellgroup = utils.convert_to_wellgroup([plate1.well(0),tube1])
        self.assertEqual(len(wellgroup),2)        
        
        wellgroup = utils.convert_to_wellgroup(tube1.well(0))
        self.assertEqual(len(wellgroup),1)  
        
        #list of containers
        
        wellgroup = utils.convert_to_wellgroup([plate1,plate1,tube1])
        self.assertEqual(len(wellgroup),13)   
        
    def test_get_column_wells(self):
        
        plate1 = create_blank_plate('6-flat')
        
        for i in range(0,len(plate1.all_wells())):
            plate1.well(i).volume = ul(i)
        
        wells = utils.get_column_wells(plate1,0).wells
        self.assertEqual(len(wells),2)
        self.assertListEqual(wells,
                             [plate1.well(0),plate1.well(3)])
        
        
        wells = utils.get_column_wells(plate1,2).wells
        self.assertEqual(len(wells),2)
        self.assertListEqual(list(wells),
                             [plate1.well(2),plate1.well(5)])
        
        
    def test_breakup_dispense_column_volumes(self):
        
        col_vol_pairs = [{"column": 0, "volume": ml(4)},
                         {"column": 1, "volume": ml(5)},
                         {"column": 2, "volume": ml(2)}
                         ]
        self.assertListEqual([{"column": 0, "volume": ml(2)},
                              {"column": 0, "volume": ml(2)},
                              {"column": 1, "volume": ml(2.5)},
                              {"column": 1, "volume": ml(2.5)},      
                              {"column": 2, "volume": ml(2)}

                              ],
                             utils.breakup_dispense_column_volumes(col_vol_pairs))
        
    
        
        
        
    def test_ceil_volume(self):
        self.assertEqual(utils.ceil_volume(ul(0.54),1),ul(0.6)) 
        self.assertEqual(utils.ceil_volume(ul(0.54)),ul(1)) 
        self.assertEqual(utils.ceil_volume(ul(1.0),1),ul(1)) 
        
    def test_round_volume(self):
        self.assertEqual(utils.round_volume(ul(0.54),1),ul(0.5))
        
        
    def test_convert_mass_to_volume(self):
        plate = create_blank_plate('24-deep')
        dna_well = plate.well(0)
        dna_well.volume = ml(10)
        
        dna_well.properties['Concentration (DNA)'] = '4:nanograms/microliter'
        
        self.assertEqual(utils.convert_mass_to_volume(utils.ng(20), 
                                                     dna_well),
                         ul(5))
        
        
    def test_convert_moles_to_volume(self):
        plate = create_blank_plate('24-deep')
        dna_well = plate.well(0)
        dna_well.volume = ml(10)
    
        dna_well.properties['Concentration (DNA)'] = '649:nanograms/microliter'
        dna_well.properties['dna_length'] = 1000
        
        self.assertEqual(utils.convert_moles_to_volume(utils.pmol(1), 
                                                       dna_well),
                         ul(1))

    def test_convert_moles_to_volume2(self):
        plate = create_blank_plate('24-deep')
        dna_well = plate.well(0)
        dna_well.volume = ml(10)
    
        dna_well.properties['Concentration (DNA)'] = '59.6:nanograms/microliter'
        dna_well.properties['dna_length'] = 5329
    
        self.assertEqual(utils.convert_moles_to_volume(utils.pmol(0.1), 
                                                       dna_well),
                         ul(5.81))    
    
        
    def test_convert_stamp_shape_to_wells(self):
        plate1 = create_blank_plate('96-flat')
        plate2 = create_blank_plate('96-flat')
        
        for well in plate1.all_wells():
            well.volume = ul(300)
        
        source_wells, dest_wells = utils.convert_stamp_shape_to_wells(plate1.well(0), plate2.well(1), 
                                           {'rows':8,'columns':1})
        
        
        self.assertListEqual(source_wells, utils.get_column_wells(plate1,0).wells )
        self.assertListEqual(dest_wells, utils.get_column_wells(plate2,1).wells )
        
        
        
    def test_calculate_dilution_volume(self):
        
        self.assertEqual(utils.calculate_dilution_volume(utils.mM(100), 
                                                        utils.uM(500), 
                                                        ul(440)), ul(2.2))