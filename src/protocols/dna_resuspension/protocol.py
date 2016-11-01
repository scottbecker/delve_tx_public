from __future__ import print_function
from transcriptic_tools.utils import ul, ml, set_property
from transcriptic_tools.harness import run
from transcriptic_tools import CustomProtocol as Protocol
from autoprotocol.protocol import Container, Well
from transcriptic_tools.enums import Temperature


PROPERTIES_TO_COPY = ['Concentration (DNA)',
                      'A260/A280']

def main(p, params):    
    assert isinstance(p, Protocol)
    
    dna_tube = params["dna_tube"]
    
    assert isinstance(dna_tube,Container)
    
    dna_tube.set_storage(Temperature.cold_80.name)
    
    #Centrifuge the tube for 3-5 sec at a minimum of 3000 x g to pellet the material to the bottom of the tube.
    p.spin(dna_tube, '3000:g', '5:second')
    
    #Add 500uL TE (to reach 1 ng/uL)
    p.provision_by_name('te', dna_tube.well(0), ul(params.get('resuspension_volume_uL',500)))
    
    #Briefly vortex (here we have to use mix) and centrifuge
    
    p.mix(dna_tube.well(0))
    
    p.spin(dna_tube, '3000:g', '3:second')    
    
    p.mix(dna_tube.well(0))
    
    if params.get('measure_concentration'):
        p.measure_concentration(dna_tube.well(0), 'dna_concentration', 'DNA', ul(2))
    
    if params.get('split_tubes'):
        
        #create a separate container and take half to be used for other experiments
        
        new_dna_tube = p.ref('%s_2'%dna_tube.name,cont_type=dna_tube.container_type,
              storage=Temperature.cold_4.name).well(0)
        
        p.transfer(dna_tube.well(0), new_dna_tube, dna_tube.well(0).volume / 2)
        
        for property_name in PROPERTIES_TO_COPY:
            if dna_tube.well(0).properties.get(property_name):
                set_property(new_dna_tube, property_name, dna_tube.well(0).properties[property_name])
    
    
if __name__ == '__main__':
    run(main, "DNAResuspension")
