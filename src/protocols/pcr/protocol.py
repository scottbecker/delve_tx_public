from __future__ import print_function
from transcriptic_tools.utils import (ul, ml, ng,
                                      convert_mass_to_volume
                                      )

from transcriptic_tools.harness import run
from transcriptic_tools import CustomProtocol as Protocol
from transcriptic_tools.enums import Antibiotic, Reagent, Temperature

def run_pcr(p, oligo1_tube, oligo2_tube,
            template_well,
            annealing_temp_c,
            greater_than_65_percent_gc_primer,
            product_length,
            product_name,
            negative_control=True,
            ):
    """
    
    
    """
    
    
    #@TODO: update to dynamically detect concentration of primer tubes (currently assumes its 100uM)
    #@TODO: include option to proceed immediately to spread/transform/pick
    
    assert isinstance(p,Protocol)
   
    p.pcr(template_well, oligo1_tube.well(0), oligo2_tube.well(0), 
         annealing_temp_c, greater_than_65_percent_gc_primer,product_length,
         product_name=product_name)
    

def main(p, params):    
    """This protocol takes a tube of bacteria, amplifies it, and creates/stores new tubes of frozen cells
    """
    #bacterial protocol
    p.mammalian_cell_mode = False
    
    run_pcr(p,params['oligo1_tube'],params['oligo2_tube'],params['template_well'],
                                       params['annealing_temp_c'],
                                       params['greater_than_65_percent_gc'],
                                       params['product_length'],
                                       params['product_name']
                                       )

if __name__ == '__main__':
    run(main, "PCR")
