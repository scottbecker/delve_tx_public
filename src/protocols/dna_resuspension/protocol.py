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
    
    dna_wells = params["dna_wells"]
    
    if params.get('send_to_80C',False):
        for dna_well in dna_wells:
            dna_well.container.set_storage(Temperature.cold_80.name)
    
    containers = set([well.container for well in dna_wells])
    
    assert len(containers)==1, 'we allow wells from one container at a time'
    container = list(containers)[0]
    
    #Centrifuge the tube for 3-5 sec at a minimum of 3000 x g to pellet the material to the bottom of the tube.
    p.spin(container, '3000:g', '10:second')
    
    #Add 500uL TE (to reach 1 ng/uL)
    p.provision_by_name('te', dna_wells, ul(params.get('resuspension_volume_uL',500)))
    
    #Briefly vortex (here we have to use mix) and centrifuge
    
    p.mix(dna_wells)
    
    p.spin(container, '3000:g', '3:second')    
    
    p.mix(dna_wells)
    
    if params.get('measure_concentration'):
        p.measure_concentration(dna_wells, 'dna_concentration', 'DNA', ul(2))
    
    if params.get('split_tubes'):
        
        #create a separate container and take half to be used for other experiments

        new_container = p.ref('%s_2'%container.name,cont_type=container.container_type,
                              storage=Temperature.cold_4.name)
        

        for dna_well in dna_wells:
        
            new_dna_well = new_container.well(dna_well.index)
            new_dna_well.name = dna_well.name
            
            p.transfer(dna_well, new_dna_well, dna_well.volume / 2)
            
            for property_name in PROPERTIES_TO_COPY:
                if dna_well.properties.get(property_name):
                    set_property(new_dna_well, property_name, dna_well.properties[property_name])
    
    
if __name__ == '__main__':
    run(main, "DNAResuspension")
