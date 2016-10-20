from __future__ import print_function
from transcriptic_tools.utils import ml, ul
from transcriptic_tools.harness import run
from transcriptic_tools.custom_protocol import CustomProtocol as Protocol
from autoprotocol.protocol import Container

def main(p, params):    
    assert isinstance(p, Protocol)
    
    #create a container 
    
    container = p.ref(params['container_name'], cont_type=params['container_type'],
                 storage=params['storage_conditions'], discard=False)
    
    p.provision_by_name('water',container.all_wells(), ul(params['volume_ul']))
    
if __name__ == '__main__':
    run(main, "CreateWaterContainer")
