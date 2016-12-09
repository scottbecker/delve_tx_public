from __future__ import print_function
from transcriptic_tools.utils import (ul, ml, get_cell_line_name, get_well_max_volume,
                                      get_column_wells,
                                      copy_cell_line_name, set_property, get_volume,
                                      set_name)

from transcriptic_tools.harness import run
from transcriptic_tools import CustomProtocol as Protocol
from transcriptic_tools.enums import Antibiotic, Reagent, Temperature

def miniprep(p, source_bacteria_well,
             antibiotic=None,
             destroy_source_tube=False,
             incubate_hours=24,
             measure_od=True,
             induce_high_copy_number=False,
             slow_growth=False
             ):
    
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
    if slow_growth:
        growth_wells = get_column_wells(growth_plate,0)
    else:
        growth_wells = [growth_plate.well(0)]
    
    
    miniprep_result_plate = p.ref('%s_miniprep'%cell_line_name, cont_type="96-pcr", storage=Temperature.cold_4)
    miniprep_result_well = miniprep_result_plate.well(0)
    
    #adding antibiotic mixes after by default
    p.add_antibiotic(growth_wells, antibiotic, broth_volume=ul(1930))    
    source_bacteria_volume = ul(50)/len(growth_wells)
    p.transfer(source_bacteria_well,growth_wells,source_bacteria_volume,mix_before=True,
               mix_after=source_bacteria_volume<=ul(10))
    
    interval_hr = 24
    
    if incubate_hours>=24:
        for t in range(0,incubate_hours/interval_hr):
            p.incubate(growth_plate, "warm_37", "{}:hour".format(interval_hr), shaking=True)
               
            if measure_od:
                p.measure_bacterial_density(growth_wells,
                                            blanking_antibiotic=antibiotic,
                                            one_tip=True)
            
    if incubate_hours%24:
        p.incubate(growth_plate, "warm_37", "{}:hour".format(incubate_hours%24), shaking=True)

        if measure_od:
            p.measure_bacterial_density(growth_wells,
                                        blanking_antibiotic=antibiotic,
                                        one_tip=True)


    if induce_high_copy_number:
        induction_wells = growth_plate.wells_from(1,len(growth_wells),columnwise=True)
    
        overnight_bacteria_volume_to_use = ul(400)
    
        p.transfer(growth_wells, induction_wells, overnight_bacteria_volume_to_use,
                   mix_before=True)
    
        #adding antibiotic mixes after by default
        p.add_antibiotic(induction_wells, antibiotic, 
                         total_volume_to_add_including_broth=ul(1998)-overnight_bacteria_volume_to_use)
    
        p.provision_by_name(Reagent.copycontrol_induction_solution, induction_wells, ul(2))

        p.incubate(growth_plate, "warm_37", "{}:hour".format(5), shaking=True)
    
        if measure_od:
            p.measure_bacterial_density(induction_wells,
                                        blanking_antibiotic=antibiotic,
                                        one_tip=True)
    
        growth_wells = induction_wells
 


    if slow_growth:

        #consolidate bacteria

        trash_wells = get_column_wells(growth_plate,11)
        
        set_name(trash_wells,'trash_well')

        p.spin(growth_plate, '4000:g', '5:minute')
        
        consolidation_well = growth_wells[0]

        trash_volumes = [get_volume(consolidation_well, aspiratable=True)]+\
            [get_volume(consolidation_well, aspiratable=True)-ul(298)]*(len(growth_wells)-1)
        p.transfer(growth_wells, trash_wells, trash_volumes)
        
        p.mix(growth_wells[1:],repetitions=20)

        p.consolidate(growth_wells[1:], consolidation_well,ul(283),allow_carryover=True,
                      mix_after=True)

        if measure_od:
            p.measure_bacterial_density(consolidation_well,
                                        blanking_antibiotic=antibiotic,
                                        one_tip=True)
  
    final_growth_well = growth_wells[0]
        
    final_growth_well.name = miniprep_result_plate.name
        
    p.miniprep(final_growth_well,miniprep_result_well)
    
    p.measure_concentration(miniprep_result_well, "dna concentration", 'DNA')
    
    if destroy_source_tube:
        
        source_bacteria_well.container.discard()
    
    

def main(p, params):    
    """This protocol takes a tube of bacteria, amplifies it, and creates/stores new tubes of frozen cells
    """
    #bacterial protocol
    p.mammalian_cell_mode = False
    
    miniprep(p, params['bacteria_well'],
             antibiotic=Antibiotic.from_string(params['antibiotic']) if params['antibiotic'] != 'cell_line' else None,
             destroy_source_tube=params['destroy_source_tube'],
             incubate_hours=params['incubate_hours'],
             measure_od=params['measure_od'],
             induce_high_copy_number=params['induce_high_copy_number'],
             slow_growth=params['slow_growth']
             )
    
if __name__ == '__main__':
    run(main, "Miniprep")
