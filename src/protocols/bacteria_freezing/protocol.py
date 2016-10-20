from __future__ import print_function
from transcriptic_tools.utils import (ul, ml, get_cell_line_name, get_well_max_volume,
                                      copy_cell_line_name, set_property, get_volume)

from transcriptic_tools.harness import run
from transcriptic_tools import CustomProtocol as Protocol
from transcriptic_tools.enums import Antibiotic, Reagent

def amplify_and_freeze_bacteria(p, source_bacteria_well,
                                antibiotic=None,
                                destroy_source_tube=False):
    
    cell_line_name = get_cell_line_name(source_bacteria_well)
    
    assert isinstance(p,Protocol)
    
    if not antibiotic:
        #get the antibiotic from the source well
        well_antibiotic_str = source_bacteria_well.properties.get('antibiotic')
        if well_antibiotic_str:
            antibiotic = Antibiotic.from_string(well_antibiotic_str)
        else:
            raise Exception('Source Well must have property \'antibiotic\' set if antibiotic not specified')
    
    assert isinstance(antibiotic,Antibiotic)
    
    # Tubes and plates
    growth_plate = p.ref('growth_plate', cont_type="96-flat", discard=True)
    growth_wells = growth_plate.wells(['A1','A2'])
    
    mix_well =  p.ref('temp_mix_tube', cont_type='micro-1.5',
                      discard=True).well(0)
    
    p.transfer(source_bacteria_well,mix_well,ul(25),mix_before=True)
    
    #adding antibiotic mixes after by default
    p.add_antibiotic(mix_well, antibiotic, broth_volume=ul(650))
    
    p.distribute(mix_well, growth_wells, ul(325), allow_carryover=True)
    
    p.cover(growth_plate)
    
    #grow bacteria until they are in their log phase of growth
    #https://www.qiagen.com/us/resources/technologies/plasmid-resource-center/growth%20of%20bacterial%20cultures/
    p.incubate(growth_plate, "warm_37", "{}:hour".format(15), shaking=True)
    
    p.uncover(growth_plate)
    
    p.absorbance(growth_plate, growth_wells.indices(),
                 wavelength="600:nanometer",
                 dataref='cell_density_600nm_absorbance', flashes=25)
    
    mix_well2 =  p.ref('temp_mix_tube2', cont_type='micro-1.5',
                       discard=True).well(0)    
    
    #add glycerol
    p.provision_by_name(Reagent.glycerol,mix_well2,ul(600))
    
    p.consolidate(growth_wells, mix_well2, get_volume(growth_wells[0],aspiratable=True),
                  mix_after=True)
    
    
    
    cryo_vial_wells = []
    for x in range(0,10):
        cryo_vial_wells.append(p.ref("cryo_frozen_%s_cells_%s"%(cell_line_name,x), 
                                     cont_type="micro-1.5", 
                                     storage="cold_80",
                                     cell_line_name=cell_line_name).well(0))
        
    copy_cell_line_name(source_bacteria_well, cryo_vial_wells)
    set_property(cryo_vial_wells,'antibiotic',antibiotic.name)
    
    p.distribute(mix_well2, cryo_vial_wells, ul(115), allow_carryover=True)
    
    if destroy_source_tube:
        
        source_bacteria_well.container.discard()
    
    

def main(p, params):    
    """This protocol takes a tube of bacteria, amplifies it, and creates/stores new tubes of frozen cells
    """
    #bacterial protocol
    p.mammalian_cell_mode = False
    
    amplify_and_freeze_bacteria(p, params['bacteria_well'],
                                Antibiotic.from_string(params['antibiotic']) if params['antibiotic'] != 'cell_line' else None,
                                params['destroy_source_tube']
                                )
    
if __name__ == '__main__':
    run(main, "BacteriaFreezing")
