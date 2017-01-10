from __future__ import print_function
from transcriptic_tools.utils import (ul, ml, get_cell_line_name, get_well_max_volume,
                                      copy_cell_line_name, set_property, get_column_wells,
                                      get_volume, get_well_dead_volume, floor_volume,
                                      round_volume, get_well_safe_volume)

from transcriptic_tools.harness import run
from transcriptic_tools import CustomProtocol as Protocol
from transcriptic_tools.enums import Antibiotic, Reagent

def pellet_bacteria(p, source_bacteria_well,
                    antibiotic=None,
                    destroy_source_tube=False,
                    induce_high_copy_number=False):
    
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
    growth_plate = p.ref('growth_plate', cont_type="96-deep", discard=True)
    
    if induce_high_copy_number:
        overnight_volume = ul(550)
        growth_wells = growth_plate.wells_from(0,3,columnwise=True)
    else:
        #only use 50% to maximize growth
        overnight_volume = ul(1000)
        growth_wells = get_column_wells(growth_plate, range(0,4))
        
    
    if destroy_source_tube:
        source_volume = get_volume(source_bacteria_well) - get_well_safe_volume(source_bacteria_well)
    else:
        source_volume = ul(50)
    
    
    #adding antibiotic mixes after by default
    p.add_antibiotic(growth_wells, antibiotic, broth_volume=overnight_volume)
    
    p.distribute(source_bacteria_well, growth_wells, round_volume(source_volume / len(growth_wells) * 1.0,2), 
                 allow_carryover=True)
    
    p.cover(growth_plate)
    
    #grow bacteria until they are in their log phase of growth
    #https://www.qiagen.com/us/resources/technologies/plasmid-resource-center/growth%20of%20bacterial%20cultures/
    p.incubate(growth_plate, "warm_37", "{}:hour".format(16), shaking=True)
    
    p.uncover(growth_plate)
    

    
    
    if induce_high_copy_number:
        induction_wells = get_column_wells(growth_plate, 1)
        
        mix_well2 =  p.ref('temp_mix_tube2', cont_type='micro-2.0',
                           discard=True).well(0)
        
        
        overnight_bacteria_volume = ul(400)
        
        p.consolidate(growth_wells,mix_well2,floor_volume(overnight_bacteria_volume / len(growth_wells)))
        
        copy_control_volume = ul(2)
        
        #adding antibiotic mixes after by default
        p.add_antibiotic(mix_well2, antibiotic, total_volume_to_add_including_broth=ml(2)-overnight_bacteria_volume-\
                         copy_control_volume)
        
        p.provision_by_name(Reagent.copycontrol_induction_solution, mix_well2, copy_control_volume)
        
        p.distribute(mix_well2, induction_wells, floor_volume(get_volume(mix_well2,aspiratable=True) / len(induction_wells)), 
                     allow_carryover=True)

        p.incubate(growth_plate, "warm_37", "{}:hour".format(5), shaking=True)
        
        p.measure_bacterial_density(induction_wells)
        
        final_growth_wells = induction_wells
        
    else:
        final_growth_wells = growth_wells
        
        top_row_well_indices = range(0,8*4,8)
        top_wells = [growth_wells[i] for i in top_row_well_indices]    
        
        non_top_wells = list(set(final_growth_wells) - set(top_wells))
            
        p.measure_bacterial_density(top_wells, one_tip=True)
        
        #spin down each column of cells
        
        p.spin(growth_plate, '4000:g', '5:minute')
        
        trash_wells = growth_plate.wells_from(4,len(final_growth_wells),columnwise=True)
        
        consolidate_volume = ul(250)
        
        
        #trash supernatent except consolidate_volume in non top row wells
        
        p.transfer(non_top_wells, trash_wells, get_volume(non_top_wells[0],aspiratable=True)-consolidate_volume,
                   one_tip=True)
        
        
        #trash all volume in the top row wells
        
        trash_wells = growth_plate.wells_from(8,4,columnwise=True)
        p.transfer(top_wells, trash_wells, get_volume(top_wells[0],aspiratable=True),
                   one_tip=True)     
        
        #break up pellet in all wells
        p.mix(non_top_wells, repetitions=20, one_tip=True)
        
        #consolidate into the top row
        for top_well in top_wells:
            non_top_wells = set(growth_plate.wells_from(top_well.index,8,columnwise=True))
            non_top_wells.remove(top_well)
            non_top_wells = list(non_top_wells)
            
            p.consolidate(non_top_wells, top_well, consolidate_volume)
            
        #spin down again and trash all but dead_volume+250ul
        p.spin(growth_plate, '4000:g', '5:minute')
        trash_wells = growth_plate.wells_from(9,4,columnwise=True)
        p.transfer(top_wells, trash_wells, get_volume(growth_wells[0],aspiratable=True)-consolidate_volume*2,
                   one_tip=True)        
    
        #break up the pellet created by spinning
        p.mix(top_wells, repetitions=20, one_tip=True)
        
        final_growth_wells = top_wells
   
    pellet_well = p.ref("pelleted_%s_cells"%(cell_line_name), 
                        cont_type="micro-2.0", 
                        storage="cold_80",
                        cell_line_name=cell_line_name).well(0)
        
    copy_cell_line_name(source_bacteria_well, pellet_well)
    
    p.consolidate(final_growth_wells, pellet_well, get_volume(final_growth_wells[0],aspiratable=True), 
                  allow_carryover=True)
    
    p.measure_bacterial_density(pellet_well)
    
    p.spin(pellet_well.container, '4000:g', '10:minute')
    
    trash_well = p.ref('trash_tube', cont_type='micro-2.0',
                       discard=True).well(0)    
    
    
    #keep ul(25) to make the safe volume
    p.transfer(pellet_well, trash_well, pellet_well.volume-ul(25), one_source=True, one_tip=True)
    p.provision_by_name(Reagent.glycerol, pellet_well, ul(15),mix_after=True)

    #remove max (this means the remaining volume is 20% glycerol)
    p.transfer(pellet_well, trash_well, get_volume(pellet_well,aspiratable=True), one_source=True, one_tip=True)
    
    if destroy_source_tube:
        
        source_bacteria_well.container.discard()
    
    

def main(p, params):    
    """This protocol takes a tube of bacteria, amplifies it, pellets the results. It can optionally induce additional BAC plasmids with IPTG
    """
    #bacterial protocol
    p.mammalian_cell_mode = False
    
    pellet_bacteria(p, params['bacteria_well'],
                    Antibiotic.from_string(params['antibiotic']) if params['antibiotic'] != 'cell_line' else None,
                    params['destroy_source_tube'],
                    params['induce_high_copy_number']
                    )
    
if __name__ == '__main__':
    run(main, "BacteriaPelleting")
