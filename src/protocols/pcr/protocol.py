from __future__ import print_function
from transcriptic_tools.utils import (ul, ml, ng,
                                      convert_mass_to_volume
                                      )

from transcriptic_tools.harness import run
from transcriptic_tools import CustomProtocol as Protocol
from transcriptic_tools.inventory import get_reagent_by_transcritpic_resource_id
from transcriptic_tools.enums import Antibiotic, Reagent, Temperature


def main(p, params):    
    """This protocol takes a tube of bacteria, amplifies it, and creates/stores new tubes of frozen cells
    """
    #bacterial protocol
    p.mammalian_cell_mode = False
    
    
    assert isinstance(p,Protocol)
    
    #determine if either primer is a reagent
    if params['primer_type1']['value'] == 'resource':
        
        resource_id = params['primer_type1']['inputs']['resource']['resource_id']
        
        primer1_well_or_reagent = get_reagent_by_transcritpic_resource_id(resource_id, fail_if_not_found=True)
        
    else:
        primer1_well_or_reagent = params['primer_type1']['inputs']['custom_inventory']['primer_well']
        
    if params['primer_type2']['value'] == 'stock':
    
        resource_id = params['primer_type2']['inputs']['stock']['resource_id']
    
        primer2_well_or_reagent = get_reagent_by_transcritpic_resource_id(resource_id, fail_if_not_found=True)
    
    else:
        primer2_well_or_reagent = params['primer_type2']['inputs']['custom_inventory']['primer_well']    
        
    
    p.pcr(template_well=params['template_well'],
          primer1_well_or_reagent = primer1_well_or_reagent,
          primer2_well_or_reagent = primer2_well_or_reagent,
          annealing_temp_c = params['annealing_temp_c'],
          greater_than_65_percent_gc_primer = params['greater_than_65_percent_gc'],
          product_length = params['product_length'],
          product_name = params['product_name']
            )

if __name__ == '__main__':
    run(main, "PCR")
