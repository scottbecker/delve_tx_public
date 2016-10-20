import random
import string
import unittest
from transcriptic_tools import CustomProtocol as Protocol
from transcriptic_tools import bacteria_utils
from transcriptic_tools.enums import Antibiotic
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

class TestBacteriaUtils(unittest.TestCase):
    pass    