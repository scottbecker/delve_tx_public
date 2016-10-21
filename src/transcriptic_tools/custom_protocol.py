import datetime
import math
from enum import Enum
from collections import defaultdict
from autoprotocol import Unit
from autoprotocol import container
from autoprotocol.container import COVER_TYPES
from autoprotocol.pipette_tools import depth, dispense_target
from autoprotocol.container_type import _CONTAINER_TYPES
from autoprotocol.protocol import (Protocol, Well, Container, is_valid_well,
                                   WellGroup, Uncover, Ref)
from transcriptic_tools.utils import (ul, get_well_dead_volume,
                                      get_well_max_volume, assert_valid_volume, ml,
                                      space_available, get_well_safe_volume, ensure_list,
                                      init_inventory_container, total_plate_volume,
                                      assert_non_negative_well, convert_to_wellgroup,
                                      get_volume, floor_volume, hours, total_plate_available_volume,
                                      breakup_dispense_column_volumes, get_column_wells,
                                      set_property, ceil_volume, init_inventory_well,
                                      touchdown_pcr, convert_stamp_shape_to_wells,
                                      convert_mass_to_volume, ug)
from lib import lists_intersect, get_dict_optional_value, get_melting_temp
from transcriptic_tools.enums import Reagent, Antibiotic, Temperature
from instruction import MiniPrep
from Bio.SeqUtils import GC
from transcriptic_tools.inventory import get_transcriptic_inventory
from autoprotocol.util import make_gel_extract_params, make_band_param


MAX_PIPETTE_TIP_VOLUME = ul(900)

class CustomProtocol(Protocol):
   
    def __init__(self,parent_protocol=None,
                 mammalian_cell_mode=True):
        #catalog inventory
        self.transcriptic_inv = get_transcriptic_inventory()
        
        
        #default to test mode (users need to change this after creating the protocol before using it)
        self.set_test_mode(True)
        
        #self.trash_tube_count = 0
        
        self._trash_plate = None
        self._trash_count = 0
        self._reagant_plate_cache = {}
        
        self.parent_protocol = parent_protocol
        self._last_incubated = []
        
        def get_zero():
            return ul(0)
        
        
        #changes default incubation w/ CO2
        self.mammalian_cell_mode = mammalian_cell_mode
        
        self.unprovisioned_stock_reagent_volumes = defaultdict(get_zero)
        
        self.absorbance_measurement_count= 0
        
        #used for measuring optical density
        self.absorbance_plate = None
        self.next_absorbance_plate_index = 0
        
        super(CustomProtocol,self).__init__()
 
    def set_test_mode(self, test_mode_on):
        self.is_test_mode = test_mode_on
        
        #@TODO: this should be managed by an external database or via the transcriptic api itself
        #@TODO: update to be smart about test mode (currently only has test reagents)
        
        #circular dependency prevents importing at the top
        from transcriptic_tools.inventory import get_our_inventory
        
        self.our_inventory = get_our_inventory(test_mode_on)    
            
               

    def use_shared_culture_medium(self):
        raise Exception('deprecated')

    def get_10ml_well(self, name, discard=False):
        """

        This method will return a new well of a shared 24-deep plate.
        Meant for materials that are meant to be incubated for a short period of time
        
        """
        curr_time = datetime.datetime.now().strftime("%m_%d_%Y_%H_%M_%S")
        ten_ml_plate = self.ref("%s_%s"%(name,curr_time), cont_type="24-deep", 
                                      storage="warm_37" if not discard else None,
                                      discard=discard)     
        return ten_ml_plate.well(0)

    def _ensure_trash_tube_has_space(self, volume_to_trash):
        """Makes sure that the trash tube is available with enough space available for volume_to_trash
        This will generate a new trash tube if there isn't enough volume available"""
        
        #if the trash tube is full, transfer to itself and use max disposal volume
        
        if not self._trash_plate:
            self._trash_plate = self.ref("trash_plate", cont_type="24-deep", discard=True) 

        if space_available(self._trash_plate) < volume_to_trash:
            self._trash_count+=1
            self._trash_plate = self.ref("trash_plate_%s"%self._trash_count,
                                         cont_type="24-deep", discard=True)   

    def transfer_from_well_to_plate(self, src_well, plate, volume, one_tip=False,
                                    first_well_index=0,
                                    new_group=False,
                                    ignore_low_volume=False,
                                    dispense_target=None):
        """ Finds available wells in plate and moves volume from well to them"""
        
        if space_available(plate, first_well_index) < volume:
            raise Exception('not enough volume available in %s for %s transfer'%(plate,volume))      
        
        remaining_volume_to_xfer = volume
        dst_wells = []
        src_wells = []
        xfer_volumes = []
        
        for dst_well in plate.all_wells():
            available_volume = space_available(dst_well)
            
            xfer_volume = min(available_volume,remaining_volume_to_xfer)
            
            if xfer_volume:
                dst_wells.append(dst_well)
                src_wells.append(src_well)
                xfer_volumes.append(xfer_volume)
                
                remaining_volume_to_xfer-=xfer_volume
            
            if not remaining_volume_to_xfer:
                break

            
        self.transfer(src_wells, dst_wells, xfer_volumes,one_tip=one_tip, new_group=new_group,
                      ignore_low_volume=ignore_low_volume,dispense_target=dispense_target)

    def _consolidate_last_pipette_operation(self):
        """
        
        Moves all transfer operations into a single group (to use the same tip)
        
        
        """
        
        last_instruction = self.instructions[-1]
        
        new_groups = []
        
        current_xfer_group = None
        for group in last_instruction.data['groups']:
            
            if len(group)>1 or group.keys()[0]!='transfer':
                if current_xfer_group:
                    new_groups.append(current_xfer_group)
                new_groups.append(group)
                current_xfer_group = None
                continue
                
            
            if not current_xfer_group:
                current_xfer_group = group
                continue
            
            current_xfer_group['transfer']+=group['transfer']
            
        if current_xfer_group:
            new_groups.append(current_xfer_group)     
            
        
        #for some reason we need to update both of these
        last_instruction.data['groups'] = new_groups 
        last_instruction.groups = new_groups  
        


    def set_container_name(self, ref, new_name):
        
        old_name = ref.name
        ref.name = new_name
        
        if old_name in self.refs:
            ref = self.refs[old_name]
        
            del self.refs[old_name]
            self.refs[new_name] = ref

    #Needed to fix this bug
    #https://github.com/autoprotocol/autoprotocol-python/issues/133    
    #@TODO: make this smart enough to also work if there are Mix, Consolidate, and Distribute steps involved
    #@TODO: also consider distribute, consolidate, and mix in the keys that are present
    def _remove_redundant_mixing(self, pipette_instruction_groups):

        transfer_groups = []
        all_transfers = True
        for xfer_group in pipette_instruction_groups:
            if xfer_group.keys() == ['transfer']:
                transfer_groups+=xfer_group['transfer']
            else:
                all_transfers = False
                
        if not transfer_groups:
            return 
        
        if not all_transfers:
            #cleanup the transfer operations in blocks to prevent 
            # distribute and consolidate from causing bad behavior
            serial_transfer_groups = []
            for xfer_group in pipette_instruction_groups:
                if xfer_group.keys() == ['transfer']:
                    serial_transfer_groups.append(xfer_group)
                elif serial_transfer_groups:
                    self._remove_redundant_mixing(serial_transfer_groups)
                    serial_transfer_groups = []
            if serial_transfer_groups:
                self._remove_redundant_mixing(serial_transfer_groups)
                serial_transfer_groups = []            
                
            return
        
        #check that sources are always sources and dests are always dests (this means that we don't need to remix)
        sources = set()
        dests = set()
        for xfer in transfer_groups:
            sources.add(self._refify(xfer['from']))
            dests.add(self._refify(xfer['to']))
            
        #we can't optimize a set of transfer groups if a source becomes a dest or visa versa right now
        if lists_intersect(sources,dests):
            return
        
        #cleanup mix_before
        dest_source_cache = {}
        for xfer in transfer_groups:
            mix_key = self._refify(xfer['from'])

            if mix_key in dest_source_cache and \
               'mix_before' in xfer:
                del xfer['mix_before']

            dest_source_cache[mix_key] = True

        #cleanup mix_after
        dest_source_cache = {}
        for xfer in transfer_groups:
            mix_key = self._refify(xfer['to'])

            if 'mix_after' in xfer:       
                dest_source_cache[mix_key] = xfer

        for xfer in transfer_groups:
            if 'mix_after' in xfer and xfer not in dest_source_cache.values() and \
            ul(xfer['volume']) >= ul(10):
                del xfer['mix_after']  
        
    def transfer_column(self,source_plate,source_column_index,dest_plate,dest_column_index,volume):
        
        assert isinstance(source_plate,Container)
        
        num_rows = source_plate.container_type.row_count() 
        
        self.stamp(get_column_wells(source_plate,source_column_index)[0], 
                get_column_wells(dest_plate,dest_column_index)[0],
                volume,
                shape={'rows':num_rows,'columns':1},
                mix_before=True)        
        
    def transfer_all_volume_evenly(self,sources,dests,**mix_kwargs):
        """
        
        This function moves all volume in source_wellsorcontainer to dest_wellsorcontainer 
        resulting in an even distribution of volume
        
        """
        sources = convert_to_wellgroup(sources)      
        dests = convert_to_wellgroup(dests)
        
        xfer_sources = WellGroup([])
        xfer_dests = WellGroup([])
        xfer_volumes = []
        
        source_volumes = [get_volume(well, aspiratable=True) for well in sources.wells]
        
        ul_per_dest = floor_volume(get_volume(sources, aspiratable=True) / len(dests))
        
        last_source_index = 0
        
        for dest in dests:
            
            ul_per_dest_remaining = ul_per_dest

            for source_index in range(last_source_index,len(sources)):
                
                volume_to_xfer_for_source_remaining = min(source_volumes[source_index],
                                                          ul_per_dest_remaining) 
                while volume_to_xfer_for_source_remaining:
                    volume_to_xfer = min(MAX_PIPETTE_TIP_VOLUME,volume_to_xfer_for_source_remaining)
                    xfer_sources.append(sources[source_index])
                    xfer_dests.append(dest)
                    xfer_volumes.append(volume_to_xfer)
                    ul_per_dest_remaining-=volume_to_xfer
                    source_volumes[source_index]-=volume_to_xfer
                    volume_to_xfer_for_source_remaining-=volume_to_xfer
                
                if not ul_per_dest_remaining:
                    break
                    
            last_source_index = source_index
            
            
        #using a new_group allows us to use simple_mode in the remove redundant mixing 
        #(which has better optimization)
        self.transfer(xfer_sources, xfer_dests, xfer_volumes,new_group=True,**mix_kwargs)
        
        
    def trash_volume(self,wellsorcontainer,volume=None,sterile=False):
        """Pipette out as much volume as possible to a trash tube.
                
            if sterile is set, we will use a different pipette tip for each transfer operation to the trash
        """        
        
        wells = convert_to_wellgroup(wellsorcontainer)
        
        wells = [well for well in wells \
                 if well.volume > get_well_dead_volume(well)]
        
        if not wells:
            return
       
        new_group_made = False
        for well in wells:
            if not volume:
                volume_to_trash = well.volume - get_well_dead_volume(well)
            else:
                volume_to_trash = volume
                
            self._ensure_trash_tube_has_space(volume_to_trash)
            self.transfer_from_well_to_plate(well, self._trash_plate, volume_to_trash,
                                             one_tip=not sterile,
                                             new_group=not new_group_made,
                                             ignore_low_volume=True,
                                             dispense_target=dispense_target(depth("ll_top", distance="0.001:meter"),
                                                                             {'start':'900:microliter/second',
                                                                              'max':'900:microliter/second'},
                                                                             ))
            new_group_made = True
    
        if not sterile:
            self._consolidate_last_pipette_operation()        
        
    
    def trash_max_volume(self, wellsorcontainer, sterile=False):
        self.trash_volume(wellsorcontainer, volume=None, sterile=sterile)
        
        
    def _generate_refill_reagent_function(self,reagent_name):
        #mostly a copy of prepare_stock_culture_medium.protocol.refill_culture_medium
        def refill_reagent_function(p, reagent_plate,volume=None):    
            assert isinstance(p, Protocol)
            assert isinstance(reagent_plate,Container)
            
            max_volume = total_plate_available_volume(reagent_plate) - \
                len(reagent_plate.all_wells())*get_well_dead_volume(reagent_plate)    
            if not volume:
                volume = max_volume
            
            if volume > max_volume:
                raise Exception('Can\'t put more than %s in the reagent plate'%max_volume)
        
            wells_to_fill = []
            reagent_volumes = [] 
            
            ul_remaining_to_fill = volume.to('microliter').magnitude
            
            for well in reagent_plate.all_wells():
                well_volume_to_dispense = min(space_available(well).magnitude,
                                              ul_remaining_to_fill +
                                              get_well_dead_volume(well).magnitude)
                
                wells_to_fill.append(well)
                
                reagent_volume = ul(well_volume_to_dispense)
                
                reagent_volumes.append(reagent_volume)
                
                ul_remaining_to_fill-=well_volume_to_dispense-get_well_dead_volume(well).magnitude
                
                if ul_remaining_to_fill == 0:
                    break
                
            p.provision_by_name(reagent_name, wells_to_fill, reagent_volumes)
            
        return refill_reagent_function
        
        
    def dispense_by_name(self, reagent, plate_or_plates, col_volumes_or_volume,
                         speed_percentage=None):
        """
        col_volumes_or_volume is a volume, we will dispense this to all wells in the plate
        """
                
        
        plate_or_plates = ensure_list(plate_or_plates)
        
        if any([plate.container_type.shortname.startswith('6-flat') for plate in plate_or_plates]):
            raise ValueError('can\'t dispense into 6-flat plates')
        
        for plate in plate_or_plates:
            
            
            if isinstance(col_volumes_or_volume,Unit):
                self.dispense_full_plate(plate,self.transcriptic_inv[reagent],col_volumes_or_volume,
                                                               speed_percentage=speed_percentage,is_resource_id=True
                                                               )   
            else:
                #dispense dpesn't accept duplicate columns
                #col_volumes_or_volume = breakup_dispense_column_volumes(col_volumes_or_volume)
                self.dispense(plate,self.transcriptic_inv[reagent],col_volumes_or_volume,
                                                    speed_percentage=speed_percentage,is_resource_id=True
                                                    )
        
        
    def _dispense_if_possible(self, reagent, dests, volume):
        """
        Attempts to convert a set of dests and volumes into matching dispense calls.  Returns True if it succeeds
        #@TODO: return a new dests and volumes that still need to be completed
        """
        if not reagent.is_dispensable:
            return False
        
        containers = set([dest.container for dest in dests])
        for container in containers:
            if 'dispense' not in container.container_type.capabilities:
                return False
            
        
        if isinstance(volume, list) and len(volume) == len(dests)\
           and all([v==volume[0] for v in volume]):
            volume = volume[0]
        
        #we can't handle different volumes
        if not isinstance(volume, Unit):
            return False
                
        if volume.to('microliter').magnitude % 20 != 0:
            return False       
        
        #not able to handle 6-flat
        if dests[0].container.container_type.shortname.startswith('6-flat'):
            return False

        #we can only deal with wells from a single contianer
        if len(containers)>1:
            return False
        
        container = dests[0].container
        
        #are there any entire plates?
        
        if dests.wells == dests[0].container.all_wells().wells:
            self.dispense_by_name(reagent, dests[0].container, 
                                 volume)
            return True
        
        
        full_column_indexes = []
        #try to find entire columns
        for col_index in range(0,container.container_type.col_count):
            if set(get_column_wells(container, col_index)).issubset(set(dests)):
                full_column_indexes.append(col_index)
            
        #if we only have full columns, dispense
        if len(dests) == len(full_column_indexes)*container.container_type.row_count():
            col_volumes = []
            for col_index in full_column_indexes:
                col_volumes.append({'column':col_index, 'volume': volume})
            self.dispense_by_name(reagent, dests[0].container, 
                                  col_volumes)
            return True            
        
        #@TODO: are there any entire columns belonging to a single plate?
        #if dests in column_sets from dests[0].container
        
        return False
        
    def provision_by_name(self, reagent_name, dests, volumes, pre_warm_minutes=None,
                          mix_after=False, **mix_kwargs):
        """ This function allows provisioning reagents by their name.
        It will handle using a pre-made stock of a reagent if we have it in stock.
        If more reagent is needed, it will handle automatically privisioning the optimal amount needed.
        
        Note: its not that smart right now but will get smarter over time.
        
        pre_warm_minutes can be the number of minutes to warm a reagent at 37C before provisioning
        """
        
        if isinstance(reagent_name, str):
            reagent_name = reagent_name.lower()
            
        if reagent_name in Reagent._member_map_:
            reagent = Reagent[reagent_name]
            
        if isinstance(reagent_name,Reagent):
            reagent = reagent_name

        dests = convert_to_wellgroup(dests)

        #----- copied from parent provision ----- 

        # Check valid well inputs
        if not is_valid_well(dests):
            raise TypeError("Destinations (dests) must be of type Well, "
                            "list of Wells, or WellGroup.")
        dests = WellGroup(dests)
        if not isinstance(reagent_name, basestring) and not isinstance(reagent_name,Reagent):
            raise TypeError("reagent_name must be a string or Reagent.")
        if not isinstance(volumes, list):
            volumes = [Unit.fromstring(volumes)] * len(dests)
        else:
            if len(volumes) != len(dests):
                raise RuntimeError("To provision a resource into multiple "
                                   "destinations with multiple volumes, the  "
                                   "list of volumes must correspond with the "
                                   "destinations in length and in order.")
            volumes = [Unit.fromstring(v) for v in volumes]
        for v in volumes:
            if not isinstance(v, (basestring, Unit)):
                raise TypeError("Volume must be a string or Unit.")
              
    
        #do we worry about getting dispense volumes perfect?
        ignore_low_volume = min(volumes) >= ul(20)
            
        #------ end of copy from parent ---- 
        
        #is this a commercial reagent that Transcriptic provides?
        if reagent in self.transcriptic_inv and not pre_warm_minutes:
            
            #@TODO: intelligently convert these to dispense where possible
            
            
            dispenced = False
            if reagent.is_dispensable:
                dispenced = self._dispense_if_possible(reagent,dests,volumes)
            
            if not dispenced:
                super(CustomProtocol,self).provision(self.transcriptic_inv[reagent],dests,volumes)
            
            if mix_after:
                self.mix(dests, volume=mix_kwargs.get('mix_vol'),
                         repetitions=mix_kwargs.get('repetitions'),
                         speed=mix_kwargs.get('flowrate'))
                
            #@TODO: move this to custom dispense and provision methods
            self._assert_valid_add_volume(dests)
            
            return
            
        #must be a custom reagent that we manage            
        if reagent not in self.our_inventory and reagent not in self.transcriptic_inv:
            raise Exception('Reagent %s not found in inventory'%reagent)
        
        if reagent not in self.our_inventory:
            reagent_info = None
            #need to have a dynamic class that automatically handles filling the reagent as needed
            #the hard part is knowing what container size to use, we may need to have this pre-set or given to us
            #how do we also ensure that this warming/cooling is right before any pipette instructions we issue
            #(applies to our reagents as well)
       
            reagent_info = {
                #@TODO: make this more intelligent to resize if more volume is necessary
                'cont_type':'96-deep',
                'refill_function':self._generate_refill_reagent_function(reagent)
            }
            
            self.our_inventory[reagent] = reagent_info
        
        else:
            reagent_info = self.our_inventory[reagent]

        #@TODO: update this to auto-detect the countainer type and storage from the reference id (or from the information stored above)
        #See https://www.evernote.com/shard/s49/nl/5433404/9b8df54d-fda1-434d-836e-0dafcb6f93e3

        
        if reagent not in self._reagant_plate_cache:
            
            if not reagent_info.get('id'):
                reagent_plate = self.ref(reagent.name,
                                         cont_type=reagent_info['cont_type'],
                                         discard=True)
                
                #we set the volume to max, we will later provision volume into these
                for well in reagent_plate.all_wells():
                    well.volume = get_well_max_volume(well)
            else:
                reagent_plate = self.ref(reagent.name,
                                         id=reagent_info['id'],
                                         cont_type=reagent_info['cont_type'],
                                         storage=reagent_info['storage'])
            

            self._reagant_plate_cache[reagent] = reagent_plate
        else:
            reagent_plate = self._reagant_plate_cache[reagent]
            
            
        #if its an unprovisioned reagent
        if not reagent_info.get('id'):
            #remember how much we need to provision 
            self.unprovisioned_stock_reagent_volumes[reagent] += sum(volumes)
        
        #refill if we don't have enough reagent for the operation
        if reagent_info.get('id') and \
           total_plate_volume(reagent_plate, aspiratable=True) < sum(volumes) \
           and 'refill_function' in reagent_info:
            reagent_info['refill_function'](self, reagent_plate)
            self._last_incubated = None
            
        #incubate as needed
        if reagent == Reagent.culture_medium and pre_warm_minutes==None:
            pre_warm_minutes = 5
                
                
        incubate_instruction_index = None
        #don't re-incubate the same plate twice
        if self._last_incubated!=reagent_plate and pre_warm_minutes:
            #we don't need co2 for reagent plates
            self.incubate(reagent_plate, 'warm_37', '%s:minute'%pre_warm_minutes, co2=0)
            incubate_instruction_index = self.get_instruction_index()
            
             
        
        requires_mixing = reagent_info.get('requires_mixing',False)
        
        #create transfer instructions to satisfy the provision request
        
        first_pippette_instruction_index = None
        mixed_well_indexes = []
        for dest, volume in zip(dests, volumes):  
            #get next well for the required volume
            #@TODO: make this smart enough to take from multiple wells
            well_volumes = self._find_well_with_volume(reagent_plate, volume)      
            
            #@TODO: update this to be smart about wells that don't have enough volume        
            #check if we have enough available
            #if we don't have enough available, provision more            
            if not well_volumes:
                raise Exception('not enough volume available of %s and unable to refill, please provision more'%reagent.name)
            
            for well,volume_to_take in well_volumes:
                
                mix_before = False
                if (requires_mixing and well.index not in mixed_well_indexes):
                    mix_before = True
                
                if volume_to_take<=ul(10):
                    mix_after=True
                
                #only use the same tip if the destination starts as empty 
                one_tip = dest.volume < get_well_dead_volume(dest)
                
                self.transfer(well, dest, volume_to_take, 
                              ignore_low_volume=ignore_low_volume,
                              mix_before = mix_before,
                              one_tip=one_tip,
                              mix_after = mix_after,
                              **mix_kwargs)
                
                if first_pippette_instruction_index==None:
                    first_pippette_instruction_index = self.get_instruction_index()
            
            
        if incubate_instruction_index!=None:
            self.add_time_constraint({"mark": incubate_instruction_index, "state": "end"},
                                     {"mark": first_pippette_instruction_index, "state": "start"},
                                     '1.0:hour')              
        
        assert_valid_volume(dests)
            
       
    def _get_mix_args(self, well, volume=None, speed=None,
                      repetitions=5,
                      kwargs_format=False,
                      mix_after=False,
                      mix_before=False,
                      mix_seconds=2,
                      mix_percent=50):
        """
        mix_seconds is the total mix duration (both pipette up and down)
        
        """
        
        assert mix_seconds>=2,'We can\'t mix faster than 1 second for half of a mix operation'
    
        if isinstance(volume,str):
            volume = Unit(volume)
    
        #default volume is 50% of well volume
        if not volume:
            volume = mix_percent / 100.0 * well.volume.to('microliter')
    
        #max allowed volume is 900uL
        volume = min(MAX_PIPETTE_TIP_VOLUME,volume)
    
        if not speed:
            speed = '%s/second'%(floor_volume(volume.to('microliter') / (mix_seconds / 2.0)))            
    
        #prevent bubbles forming by keeping the pipette time to 1 second minimum            
        if ul(Unit(speed).to('microliter/second').magnitude) > volume:
            speed = '%s/second'%volume
    
        #this is the max speed
        if ul(Unit(speed).to('microliter/second').magnitude) > ul(1000):
            speed = '1000:microliter/second'
    
        if ul(Unit(speed).to('microliter/second').magnitude) < ul(2.5):
            if volume < ul(2.5):
                raise Exception('can\'t mix less than 2.5uL of volume without creating bubbles')
    
            speed = '2.5:microliter/second'    
            
        if not kwargs_format:
            return volume, speed, repetitions
        else:
            kwargs = {'mix_vol':volume, 
                      'repetitions': repetitions,
                      'flowrate': speed
                      }
            
            postfix = ''
            if mix_after:
                postfix = '_a'
            elif mix_before:
                postfix = '_b'
            
            if postfix:
                new_kwargs = {}
                for key in kwargs.keys():
                    new_kwargs[key+postfix] = kwargs[key]
                    
                kwargs = new_kwargs
                
            return kwargs 
                    
         
    def mix(self, well_well_group_or_plate, volume=None, speed=None,
                repetitions=5, one_tip=True, mix_seconds=2, mix_percent=50):
            """
            Mix specified well using a new pipette tip
            
            volume defaults to 50%
            speed defaults to complete the mix in 1 second
        
            Example Usage:
        
            .. code-block:: python
        
                p = Protocol()
                sample_source = p.ref("sample_source",
                                      None,
                                      "micro-1.5",
                                      storage="cold_20")
        
                p.mix(sample_source.well(0), volume="200:microliter",
                      repetitions=25)
        
            Autoprotocol Output:
        
            .. code-block:: json
        
                "instructions": [
                    {
                      "groups": [
                        {
                          "mix": [
                            {
                              "volume": "200:microliter",
                              "well": "sample_source/0",
                              "repetitions": 25,
                              "speed": "100:microliter/second"
                            }
                          ]
                        }
                      ],
                      "op": "pipette"
                    }
                  ]
                }
        
        
            Parameters
            ----------
            well : Well, WellGroup, list of Wells
                Well(s) to be mixed. If a WellGroup is passed, each well in the
                group will be mixed using the specified parameters.
            volume : str, Unit, optional
                volume of liquid to be aspirated and expelled during mixing
            speed : str, Unit, optional
                flowrate of liquid during mixing
            repetitions : int, optional
                number of times to aspirate and expell liquid during mixing
            one_tip : bool
                mix all wells with a single tip
        
            """   
            
            #prevent repetitions from being 0 or None
            if not repetitions:
                repetitions = 5
            
            #if all wells have the same volume, we can use an optimized version
            wells = convert_to_wellgroup(well_well_group_or_plate)
            first_well_volume = wells[0].volume
            all_wells_same_volume = True
            for well in wells:
                if well.volume!=first_well_volume:
                    all_wells_same_volume=False
                    break
            
            if all_wells_same_volume:
                volume, speed, repetitions = self._get_mix_args(wells[0], volume, 
                                                                speed, repetitions,
                                                                mix_percent=mix_percent,
                                                                mix_seconds=mix_seconds)
            
                super(CustomProtocol,self).mix(wells, volume, speed, repetitions, one_tip)
            else:
                for well in wells:
                    self.mix(well, volume, 
                            speed, 
                            repetitions, 
                            one_tip, 
                            mix_seconds, 
                            mix_percent)
            
        
    def _find_well_with_volume(self, plate, volume):
        well_volumes = []
        
        remaining_volume = volume
        
        for well in plate.all_wells(columnwise=True):
            assert isinstance(well,Well)

            #we can't draw from a well that is under safe volume
            if well.volume < get_well_safe_volume(well):
                continue
            
            #don't bring the well under its dead volume
            volume_to_take = min(well.volume - get_well_dead_volume(well),remaining_volume)
            
            well_volumes.append([well,volume_to_take])
            remaining_volume-=volume_to_take
            
            if not remaining_volume: 
                return well_volumes
                
        if remaining_volume:
            return None
        
        return well_volumes
       
    def ref(self, name, id=None, cont_type=None, storage=None, discard=None, cover=None,
            **properties):
        """
        
        Add a Ref object to the dictionary of Refs associated with this protocol
        and return a Container with the id, container type and storage or
        discard conditions specified.
        
        This custom version allows discard=True to overide storage
    
        Example Usage:
    
        .. code-block:: python
    
            p = Protocol()
    
            # ref a new container (no id specified)
            sample_ref_1 = p.ref("sample_plate_1",
                                 cont_type="96-pcr",
                                 discard=True)
    
            # ref an existing container with a known id
            sample_ref_2 = p.ref("sample_plate_2",
                                 id="ct1cxae33lkj",
                                 cont_type="96-pcr",
                                 storage="ambient")
    
        Autoprotocol Output:
    
        .. code-block:: json
    
            {
              "refs": {
                "sample_plate_1": {
                  "new": "96-pcr",
                  "discard": true
                },
                "sample_plate_2": {
                  "id": "ct1cxae33lkj",
                  "store": {
                    "where": "ambient"
                  }
                }
              },
              "instructions": []
            }
    
        Parameters
        ----------
        name : str
            name of the container/ref being created.
        id : str
            id of the container being created, from your organization's
            inventory on http://secure.transcriptic.com.  Strings representing
            ids begin with "ct".
        cont_type : str, ContainerType
            container type of the Container object that will be generated.
        storage : {"ambient", "cold_20", "cold_4", "warm_37"}, optional
            temperature the container being referenced should be stored at
            after a run is completed.  Either a storage condition must be
            specified or discard must be set to True.
        discard : bool, optional
            if no storage condition is specified and discard is set to True,
            the container being referenced will be discarded after a run.
    
        Returns
        -------
        container : Container
            Container object generated from the id and container type
             provided.
    
        Raises
        ------
        RuntimeError
            If a container previously referenced in this protocol (existant
            in refs section) has the same name as the one specified.
        RuntimeError
            If no container type is specified.
        RuntimeError
            If no valid storage or discard condition is specified.
    
        """    
        
        if not isinstance(name,basestring):
            raise ValueError('invalid name: %s'%name)
        
        if storage:
            discard=False
        
        if discard:
            storage=None
        
        if isinstance(storage, Temperature):
            storage = storage.name        
            
        ref = super(CustomProtocol,self).ref(name, id=id, cont_type=cont_type, 
                                              storage=storage, discard=discard,
                                              cover=cover)
        
        assert isinstance(ref,Container)
        
        #this is a new container so we should initialize all well's to 0 volume
        if not id:
            for well in ref.all_wells():
                well.volume = ul(0)
        else:
            init_inventory_container(ref)
        
        for property_name, property_value in properties.items():
            set_property(ref, property_name, property_value)
        
        return ref
    
    
    def transfer(self, source, dest, volume, one_source=False, one_tip=True, 
                aspirate_speed=None, dispense_speed=None, 
                aspirate_source=None, dispense_target=None, 
                pre_buffer=None, disposal_vol=None, 
                transit_vol=None, blowout_buffer=None, 
                tip_type=None, new_group=False, ignore_low_volume=False,
                **mix_kwargs):
        """
        Transfer liquid from one specific well to another.  A new pipette tip
        is used between each transfer step unless the "one_tip" parameter
        is set to True.

        Example Usage:

        .. code-block:: python

            p = Protocol()
            sample_plate = p.ref("sample_plate",
                                 ct32kj234l21g,
                                 "96-flat",
                                 storage="warm_37")


            # a basic one-to-one transfer:
            p.transfer(sample_plate.well("B3"),
                       sample_plate.well("C3"),
                       "20:microliter")

            # using a basic transfer in a loop:
            for i in xrange(1, 12):
              p.transfer(sample_plate.well(i-1),
                         sample_plate.well(i),
                         "10:microliter")

            # transfer liquid from each well in the first column of a 96-well
            # plate to each well of the second column using a new tip and
            # a different volume each time:
            volumes = ["5:microliter", "10:microliter", "15:microliter",
                       "20:microliter", "25:microliter", "30:microliter",
                       "35:microliter", "40:microliter"]

            p.transfer(sample_plate.wells_from(0,8,columnwise=True),
                       sample_plate.wells_from(1,8,columnwise=True),
                       volumes)

            # transfer liquid from wells A1 and A2 (which both contain the same
            # source) into each of the following 10 wells:
            p.transfer(sample_plate.wells_from("A1", 2),
                       sample_plate.wells_from("A3", 10),
                       "10:microliter",
                       one_source=True)

            # transfer liquid from wells containing the same source to multiple
            # other wells without discarding the tip in between:
            p.transfer(sample_plate.wells_from("A1", 2),
                       sample_plate.wells_from("A3", 10),
                       "10:microliter",
                       one_source=True,
                       one_tip=True)


        Parameters
        ----------
        source : Well, WellGroup
            Well or wells to transfer liquid from.  If multiple source wells
            are supplied and one_source is set to True, liquid will be
            transfered from each source well specified as long as it contains
            sufficient volume. Otherwise, the number of source wells specified
            must match the number of destination wells specified and liquid
            will be transfered from each source well to its corresponding
            destination well.
        dest : Well, WellGroup
            Well or WellGroup to which to transfer liquid.  The number of
            destination wells must match the number of source wells specified
            unless one_source is set to True.
            You can have 1 dest well and multiple source wells (equivalent to a consolidate)
        volume : str, Unit, list
            The volume(s) of liquid to be transferred from source wells to
            destination wells.  Volume can be specified as a single string or
            Unit, or can be given as a list of volumes.  The length of a list
            of volumes must match the number of destination wells given unless
            the same volume is to be transferred to each destination well.
        one_source : bool, optional
            Specify whether liquid is to be transferred to destination wells
            from a group of wells all containing the same substance.
        one_tip : bool, optional
            Specify whether all transfer steps will use the same tip or not.
        mix_after : bool, optional
            Specify whether to mix the liquid in the destination well after
            liquid is transferred.
        mix_before : bool, optional
            Specify whether to mix the liquid in the source well before
            liquid is transferred.
        mix_vol : str, Unit, optional
            Volume to aspirate and dispense in order to mix liquid in a wells
            before and/or after each transfer step.
        repetitions : int, optional
            Number of times to aspirate and dispense in order to mix
            liquid in well before and/or after each transfer step.
        flowrate : str, Unit, optional
            Speed at which to mix liquid in well before and/or after each
            transfer step.
        aspirate speed : str, Unit, optional
            Speed at which to aspirate liquid from source well.  May not be
            specified if aspirate_source is also specified. By default this is
            the maximum aspiration speed, with the start speed being half of
            the speed specified.
        dispense_speed : str, Unit, optional
            Speed at which to dispense liquid into the destination well.  May
            not be specified if dispense_target is also specified.
        aspirate_source : fn, optional
            Can't be specified if aspirate_speed is also specified.
        dispense_target : fn, optional
            Same but opposite of  aspirate_source.
        pre_buffer : str, Unit, optional
            Volume of air aspirated before aspirating liquid.
        disposal_vol : str, Unit, optional
            Volume of extra liquid to aspirate that will be dispensed into
            trash afterwards.
        transit_vol : str, Unit, optional
            Volume of air aspirated after aspirating liquid to reduce presence
            of bubbles at pipette tip.
        blowout_buffer : bool, optional
            If true the operation will dispense the pre_buffer along with the
            dispense volume. Cannot be true if disposal_vol is specified.
        tip_type : str, optional
            Type of tip to be used for the transfer operation.
        new_group : bool, optional

        Raises
        ------
        RuntimeError
            If more than one volume is specified as a list but the list length
            does not match the number of destination wells given.
        RuntimeError
            If transferring from WellGroup to WellGroup that have different
            number of wells and one_source is not True.

        """
        
        source = convert_to_wellgroup(source)
            
        dest = convert_to_wellgroup(dest)
        
        if len(dest) == 1 and len(source) > 1:
            dest = WellGroup([dest[0]*len(source)])
        
        if type(volume) == list:
            min_volume = min(volume)
        else:
            min_volume = volume        
        
        if not ignore_low_volume and min_volume < ul(10) and not mix_kwargs.get('mix_after') \
           and not mix_kwargs.get('ignore_mix_after_warning'):
            raise Exception('mix_after required for <10uL of solution to ensure complete transfer. \n'
                            'Ensure you have are pipetting into something big enough and set this')
            
        if 'ignore_mix_after_warning' in mix_kwargs:
            del mix_kwargs['ignore_mix_after_warning']
        
        self._complete_mix_kwargs(mix_kwargs, source, dest, volume,
                                  allow_different_mix_options=True)
            
        
        self._assert_safe_to_pipette_from(source)
            
        super(CustomProtocol,self).transfer(source, dest, volume, one_source=one_source, one_tip=one_tip, 
              aspirate_speed=aspirate_speed, dispense_speed=dispense_speed, 
              aspirate_source=aspirate_source, dispense_target=dispense_target, 
              pre_buffer=pre_buffer, disposal_vol=disposal_vol, 
              transit_vol=transit_vol, blowout_buffer=blowout_buffer, 
              tip_type=tip_type, new_group=new_group,**mix_kwargs)
        
        self._remove_redundant_mixing(self.instructions[-1].data['groups'])
        
        self._assert_valid_transfer(source, dest)
        
        #if len(self.instructions)>=155:
            #raise Exception('here')
        
    def _complete_mix_kwargs(self, mix_kwargs, sources, dests, volume, one_source=False,
                             allow_different_mix_options=False
                             ):
        """
        This function fills out the remaining mix_kwargs if only mix_before or mix_after have been specified
        """
        
        original_mix_kwargs = mix_kwargs.copy()
        
        dests = ensure_list(dests)
        sources = ensure_list(sources)
        
        sources = WellGroup(sources)
        dests = WellGroup(dests)        
        len_source = len(sources.wells)
        len_dest = len(dests.wells)        
        
        # Auto-generate list from single volume, check if list length matches
        if isinstance(volume, basestring) or isinstance(volume, Unit):
            if len_dest == 1 and not one_source:
                volume = [Unit.fromstring(volume).to("ul")] * len_source
            else:
                volume = [Unit.fromstring(volume).to("ul")] * len_dest
        elif isinstance(volume, list) and len(volume) == len_dest:
            volume = list(
                map(lambda x: Unit.fromstring(x).to("ul"), volume))
        else:
            raise RuntimeError("Unless the same volume of liquid is being "
                               "transferred to each destination well, each "
                               "destination well must have a corresponding "
                               "volume in the form of a list.")        
        
        

        
        if 'mix_after' in mix_kwargs and \
            not lists_intersect(['mix_vol','flowrate',
                                'repetitions'],
                               mix_kwargs.keys()) :
            
            #needed to calculate proper mixing speed
            future_dest_well = Well(None, 0)
            future_dest_well.volume = ul(0)
            
            #we use the first destination to determine mix volume
            
            last_dest = None
            for well, volume in zip(dests,volume):
                if last_dest and last_dest != well:
                    break
                if not last_dest:
                    future_dest_well.volume+=well.volume
                
                future_dest_well.volume += volume
                
            mix_kwargs.update(self._get_mix_args(future_dest_well, kwargs_format=True,
                                                 mix_after=allow_different_mix_options,
                                                 mix_seconds=get_dict_optional_value(mix_kwargs,
                                                                                          ['mix_seconds','mix_seconds_a'],
                                                                                          2),
                                                 mix_percent=get_dict_optional_value(mix_kwargs,
                                                                                     ['mix_percent','mix_percent_a'],
                                                                                     50)
                                                 ))
            
        if 'mix_before' in mix_kwargs and \
            not lists_intersect(['mix_vol','flowrate',
                                'repetitions'],
                               mix_kwargs.keys()):
            mix_kwargs.update(self._get_mix_args(sources[0], kwargs_format=True,
                                                 mix_before=allow_different_mix_options,
                                                 mix_seconds=get_dict_optional_value(mix_kwargs,
                                                                                     ['mix_seconds','mix_seconds_b'],
                                                                                     2),
                                                 mix_percent=get_dict_optional_value(mix_kwargs,
                                                                                     ['mix_percent','mix_percent_b'],
                                                                                     50)                                                 
                                                 ))
            
        #this fills in any parameters that we specified
        mix_kwargs.update(original_mix_kwargs)        
        
        #delete our extra mix_kwargs that Transcriptic doesn't understand
        keys_to_remove = ['mix_seconds','mix_percent','mix_percent_a',
                          'mix_seconds_a','mix_seconds_b','mix_percent_b']
        
        for key_to_remove in keys_to_remove:
            if key_to_remove in mix_kwargs:
                del mix_kwargs[key_to_remove]
        
    def distribute(self, source, dest, volume, allow_carryover=False,
                   aspirate_speed=None,
                   aspirate_source=None, dispense_speed=None,
                   distribute_target=None, pre_buffer=None, disposal_vol=None,
                   transit_vol=None, blowout_buffer=None, tip_type=None,
                   ignore_mix_after_warning = False,
                   new_group=False,**mix_kwargs):
        """
        Distribute liquid from source well(s) to destination wells(s).


        Example Usage:

        .. code-block:: python

            p = Protocol()
            sample_plate = p.ref("sample_plate",
                                 None,
                                 "96-flat",
                                 storage="warm_37")
            sample_source = p.ref("sample_source",
                                  "ct32kj234l21g",
                                  "micro-1.5",
                                  storage="cold_20")

            p.distribute(sample_source.well(0),
                         sample_plate.wells_from(0,8,columnwise=True),
                         "200:microliter",
                         mix_before=True,
                         mix_vol="500:microliter",
                         repetitions=20)

        Autoprotocol Output:

        .. code-block:: json

            "instructions": [
              {
                "groups": [
                  {
                    "distribute": {
                      "to": [
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/0"
                        },
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/12"
                        },
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/24"
                        },
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/36"
                        },
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/48"
                        },
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/60"
                        },
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/72"
                        },
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/84"
                        }
                      ],
                      "from": "sample_source/0",
                      "mix_before": {
                        "volume": "500:microliter",
                        "repetitions": 20,
                        "speed": "100:microliter/second"
                      }
                    }
                  }
                ],
                "op": "pipette"
              }
            ]

        Parameters
        ----------
        source : Well, WellGroup
            Well or wells to distribute liquid from.  If passed as a WellGroup
            with set_volume() called on it, liquid will be automatically be
            drawn from the wells specified using the fill_wells function.
        dest : Well, WellGroup
            Well or wells to distribute liquid to.
        volume : str, Unit, list
            Volume of liquid to be distributed to each destination well.  If a
            single string or unit is passed to represent the volume, that
            volume will be distributed to each destination well.  If a list of
            volumes is provided, that volume will be distributed to the
            corresponding well in the WellGroup provided. The length of the
            volumes list must therefore match the number of wells in the
            destination WellGroup if destination wells are recieving different
            volumes.
        allow_carryover : bool, optional
            specify whether the same pipette tip can be used to aspirate more
            liquid from source wells after the previous volume aspirated has
            been depleted.
        mix_before : bool, optional
            Specify whether to mix the liquid in the destination well before
            liquid is transferred.
        mix_vol : str, Unit, optional
            Volume to aspirate and dispense in order to mix liquid in a wells
            before liquid is distributed.
        repetitions : int, optional
            Number of times to aspirate and dispense in order to mix
            liquid in a well before liquid is distributed.
        flowrate : str, Unit, optional
            Speed at which to mix liquid in well before liquid is distributed.
        aspirate_speed : str, Unit, optional
            Speed at which to aspirate liquid from source well.  May not be
            specified if aspirate_source is also specified. By default this is
            the maximum aspiration speed, with the start speed being half of
            the speed specified.
        aspirate_source : fn, optional
            Can't be specified if aspirate_speed is also specified.
        dispense_speed : str, Unit, optional
            Speed at which to dispense liquid into the destination well.  May
            not be specified if dispense_target is also specified.
        distribute_target : fn, optional
            A function that contains additional parameters for distributing to
            target wells including depth, dispense_speed, and calibrated
            volume.
            If this parameter is specified, the same parameters will be
            applied to every destination well.
            Will supersede dispense_speed parameters if also specified.
        pre_buffer : str, Unit, optional
            Volume of air aspirated before aspirating liquid.
        disposal_vol : str, Unit, optional
            Volume of extra liquid to aspirate that will be dispensed into
            trash afterwards.
        transit_vol : str, Unit, optional
            Volume of air aspirated after aspirating liquid to reduce presence
            of bubbles at pipette tip.
        blowout_buffer : bool, optional
            If true the operation will dispense the pre_buffer along with the
            dispense volume.
            Cannot be true if disposal_vol is specified.

        Raises
        ------
        RuntimeError
            If no mix volume is specified for the mix_before instruction.
        ValueError
            If source and destination well(s) is/are not expressed as either
            Wells or WellGroups.

        """
        
        if isinstance(source,Container):
            if len(source.all_wells())!=1:
                raise Exception('distribute only accepts containers with one well')
            
            source = source.well(0)
            
        if isinstance(dest,Container):
            dest = dest.all_wells()        
        
        #bug: https://github.com/autoprotocol/autoprotocol-python/issues/135
        if mix_kwargs.get('mix_before') and mix_kwargs.get('mix_after') and not mix_kwargs.get('mix_vol'):
            raise Exception('can\'t specify mix_before and mix_after set True to distribute without also'+\
                            ' setting mix_vol')
        
        self._complete_mix_kwargs(mix_kwargs,source,dest,volume)
        self._assert_safe_to_pipette_from(source)

        if volume < ul(10):
            for dest_well in dest:
                self.transfer(source,dest_well,volume,
                              aspirate_speed=aspirate_speed, 
                              aspirate_source=aspirate_source, 
                              pre_buffer=pre_buffer, disposal_vol=disposal_vol, 
                              transit_vol=transit_vol, blowout_buffer=blowout_buffer,
                              tip_type=tip_type, new_group=new_group,
                              ignore_mix_after_warning=ignore_mix_after_warning,
                              **mix_kwargs)                               
        else:
            
            if volume>ul(900):
                repeat_count = int(math.ceil(volume.to('microliter').magnitude / 900.0))
                volume = floor_volume(volume/repeat_count)
                 
                dest = repeat_count*ensure_list(dest)
                
                
            super(CustomProtocol,self).distribute(source, dest, volume, allow_carryover=allow_carryover,
                                                  aspirate_speed=aspirate_speed,
                                                  aspirate_source=aspirate_source, 
                                                  dispense_speed=dispense_speed,
                                                  distribute_target=distribute_target,
                                                  pre_buffer=pre_buffer, disposal_vol=disposal_vol, 
                                                  transit_vol=transit_vol,
                                                  blowout_buffer=blowout_buffer, tip_type=tip_type, 
                                                  new_group=new_group,
                                                  **mix_kwargs)
            
        self._assert_valid_transfer(source, dest)
            
    def incubate(self, ref, where, duration, shaking=None, co2=None):
        """
                Move plate to designated thermoisolater or ambient area for incubation
                for specified duration.
        
                Example Usage:
        
                .. code-block:: python
        
                    p = Protocol()
                    sample_plate = p.ref("sample_plate",
                                         None,
                                         "96-pcr",
                                         storage="warm_37")
        
                    # a plate must be sealed/covered before it can be incubated
                    p.seal(sample_plate)
                    p.incubate(sample_plate, "warm_37", "1:hour", shaking=True)
        
                Autoprotocol Output:
        
                .. code-block:: json
        
                    "instructions": [
                        {
                          "object": "sample_plate",
                          "op": "seal"
                        },
                        {
                          "duration": "1:hour",
                          "where": "warm_37",
                          "object": "sample_plate",
                          "shaking": true,
                          "op": "incubate",
                          "co2_percent": 5
                        }
                      ]
        
                """   
        if self.mammalian_cell_mode:
            if co2 == None:
                co2 = 5
        else:
            if shaking==None:
                shaking = True
        
        if not shaking:
            shaking=False
        
        if isinstance(ref,list):
            refs = ref
            for ref in refs:
                self.incubate(ref, where, duration, shaking, co2)
                
            return
        

        if isinstance(where, Temperature):
            where = where.name

        self._last_incubated = ref
        
        super(CustomProtocol,self).incubate(ref,where,duration,shaking,co2)
        
        
    def _assert_safe_to_pipette_from(self, sources):
        sources = ensure_list(sources)
        if not all([well.volume >= get_well_safe_volume(well) for well in sources]):
            raise Exception('Not safe to pipette from %s'%sources)
        
    def _assert_valid_add_volume(self, dests):
        dests = ensure_list(dests)
        if not all([well.volume <= get_well_max_volume(well) for well in dests]):
            raise Exception('Too much liquid added to dests: %s'%dests)            

    def _assert_valid_transfer(self, sources, dests):
        sources = ensure_list(sources)
        if not all([well.volume >= get_well_dead_volume(well) for well in sources]):
            raise Exception('Too much liquid drawn from sources: %s'%sources)
        self._assert_valid_add_volume(dests) 

    def stamp(self, source_origin, dest_origin, volume, shape=dict(rows=8,
                                                                   columns=12),
             
              aspirate_speed=None, dispense_speed=None, aspirate_source=None,
              dispense_target=None, pre_buffer=None, disposal_vol=None,
              transit_vol=None, blowout_buffer=None, one_source=False,
              one_tip=False, new_group=False,
              **mix_kwargs):
        """
        **Note: the way this method now works is significantly different to the
        way it has in previous versions, please make sure to read the
        documentation below and adjust existing scripts utilizing stamp()
        accordingly**
    
        A stamp instruction consists of a list of groups of transfers, each of
        which specifies from and to well references (ref/well_index)
        representing the top-left well or origin of a specified shape.
    
        The volume field defines the volume of liquid that will be aspirated
        from every well of the shape specified starting at the from field and
        dispensed into the corresponding wells starting at the to field.
    
        Currently, the shape field may only be a rectangle object defined by
        rows and columns attributes representing the number of contiguous tip
        rows and columns to transfer.
    
        The shape parameter is optional and will default to a full 8 rows by
        12 columns. The tip_layout field refers to the SBS compliant layout of
        tips, is optional, and will default to the layout of a 96 tip box.
    
        The following plate types are currently supported: 96 and 384.
    
    
        Example Usage:
    
        .. code-block:: python
    
            p = Protocol()
    
            plate_1_96 = p.ref("plate_1_96", None, "96-flat", discard=True)
            plate_2_96 = p.ref("plate_2_96", None, "96-flat", discard=True)
            plate_1_384 = p.ref("plate_1_384", None, "384-flat", discard=True)
            plate_2_384 = p.ref("plate_2_384", None, "384-flat", discard=True)
    
            # A full-plate transfer between two 96 or 384-well plates
            p.stamp(plate_1_96, plate_2_96, "10:microliter")
            p.stamp(plate_1_384, plate_2_384, "10:microliter")
    
            # Defining shapes for selective stamping:
            row_rectangle = dict(rows=1, columns=12)
            two_column_rectangle = dict(rows=8, columns=2)
    
            # A transfer from the G row to the H row of another 96-well plate
            p.stamp(plate_1_96.well("G1"), plate_2_96.well("H1"),
            "10:microliter", row_rectangle)
    
            # A 2-column transfer from columns 1,2 of a 96-well plate to
            #columns 2,4 of a 384-well plate
            p.stamp(plate_1_96.well("A1"), plate_1_384.wells_from("A2", 2,
            columnwise=True), "10:microliter", two_column_rectangle)
    
            # A 2-row transfer from rows 1,2 of a 384-well plate to rows 2,3
            #of a 96-well plate
            p.stamp(plate_1_384.wells(["A1", "A2", "B1", "B2"]), plate_1_96.
            wells(["B1", "B1", "C1", "C1"]), "10:microliter",
            shape=row_rectangle)
    
        Parameters
        ----------
        source_origin : Container, Well, WellGroup, List of Wells
            Top-left well or wells where the rows/columns will be defined with
            respect to the source transfer.
            If a container is specified, stamp will be applied to all
            quadrants of the container.
        dest_origin : Container, Well, WellGroup, List of Wells
            Top-left well or wells where the rows/columns will be defined with
            respect to the destination transfer.
            If a container is specified, stamp will be applied to all
            quadrants of the container
        volume : str, Unit, list
            Volume(s) of liquid to move from source plate to destination
            plate. Volume can be specified as a single string or Unit, or can
            be given as a list of volumes.  The length of a list of volumes
            must match the number of destination wells given unless the same
            volume is to be transferred to each destination well.
        shape : dictionary, list, optional
            The shape(s) parameter is optional and will default to a rectangle
            corresponding to a full 96-well plate (8 rows by 12 columns).
            The rows and columns will be defined wrt the specified origin.
            The length of a list of shapes must match the number of
            destination wells given unless the same shape is to be used for
            each destination well. If the length of shape is greater than 1,
            one_tip=False.
    
            Example
    
            .. code-block:: python
    
                rectangle = {}
                rectangle["rows"] = 8
                rectangle["columns"] = 12
    
        mix_after : bool, optional
            Specify whether to mix the liquid in destination wells after
            liquid is transferred.
        mix_before : bool, optional
            Specify whether to mix the liquid in source wells before
            liquid is transferred.
        mix_vol : str, Unit, optional
            Volume to aspirate and dispense in order to mix liquid in wells
            before and/or after it is transfered.
        repetitions : int, optional
            Number of times to aspirate and dispense in order to mix
            liquid in wells before and/or after it is transfered.
        flowrate : str, Unit, optional
            Speed at which to mix liquid in well before and/or after each
            transfer step in units of "microliter/second".
        dispense_speed : str, Unit, optional
            Speed at which to dispense liquid into the destination well.  May
            not be specified if dispense_target is also specified.
        aspirate_source : fn, optional
            Can't be specified if aspirate_speed is also specified.
        dispense_target : fn, optional
            Same but opposite of  aspirate_source.
        pre_buffer : str, Unit, optional
            Volume of air aspirated before aspirating liquid.
        disposal_vol : str, Unit, optional
            Volume of extra liquid to aspirate that will be dispensed into
            trash afterwards.
        transit_vol : str, Unit, optional
            Volume of air aspirated after aspirating liquid to reduce presence
            of bubbles at pipette tip.
        blowout_buffer : bool, optional
            If true the operation will dispense the pre_buffer along with the
            dispense volume. Cannot be true if disposal_vol is specified.
        one_source : bool, optional
            Specify whether liquid is to be transferred to destination origins
            from a group of origins all containing the same substance. Volume
            of all wells in the shape must be equal to or greater than the
            volume in the origin well. Specifying origins with overlapping
            shapes can produce undesireable effects.
        one_tip : bool, optional
            Specify whether all transfer steps will use the same tip or not.
            If multiple different shapes are used, one_tip cannot be true.
        new_group : bool, optional
    
            Example
    
            .. code-block:: python
    
                p.stamp(plate_1_96.well("A1"), plate_2_96.well("A1"),
                "10:microliter")
                p.stamp(plate_1_96.well("A1"), plate_2_96.well("A1"),
                "10:microliter")
    
            Autoprotocol Output:
    
            .. code-block:: json
    
                "instructions": [
                    {
                      "groups": [
                        {
                          "transfer": [
                            {
                              "volume": "10.0:microliter",
                              "to": "plate_2_96/0",
                              "from": "plate_1_96/0"
                            }
                          ],
                          "shape": {
                            "rows": 8,
                            "columns": 12
                          },
                          "tip_layout": 96
                        }
                      ],
                      "op": "stamp"
                    },
                    {
                      "groups": [
                        {
                          "transfer": [
                            {
                              "volume": "10.0:microliter",
                              "to": "plate_2_96/0",
                              "from": "plate_1_96/0"
                            }
                          ],
                          "shape": {
                            "rows": 8,
                            "columns": 12
                          },
                          "tip_layout": 96
                        }
                      ],
                      "op": "stamp"
                    }
                  ]
    
        """
        
        
        self._complete_mix_kwargs(mix_kwargs,source_origin,dest_origin,volume)
        
        source_wells, dest_wells = convert_stamp_shape_to_wells(source_origin, 
                                                               dest_origin, 
                                                               shape=shape, 
                                                               one_source=one_source)
        
        
        
        
        self._assert_safe_to_pipette_from(source_wells)
        
        super(CustomProtocol,self).stamp(source_origin=source_origin, dest_origin=dest_origin, 
                                         volume=volume, 
                                         shape=shape, 
                                         aspirate_speed=aspirate_speed, dispense_speed=dispense_speed, aspirate_source=aspirate_source,
                                         dispense_target=dispense_target, pre_buffer=pre_buffer, disposal_vol=disposal_vol,
                                         transit_vol=transit_vol, blowout_buffer=blowout_buffer, one_source=one_source,
                                         one_tip=one_tip, new_group=new_group,
                                         **mix_kwargs)
        
        self._assert_valid_transfer(source_wells, dest_wells)
        

        
    def consolidate(self, sources, dest, volumes, allow_carryover=True,
                        aspirate_speed=None, dispense_speed=None, aspirate_source=None,
                        dispense_target=None, pre_buffer=None, transit_vol=None,
                        blowout_buffer=None, tip_type=None, new_group=False,
                        **mix_kwargs):
            """
            Aspirates from each source well, in order, the volume specified, then
            dispenses the sum volume into the target well. Be aware that the same
            tip will be used to aspirate from all the source wells, so if you want
            to avoid contaminating any of them you should use a separate transfer
            group. Consolidate is limited by the maximum volume of the disposable
            tip. If the total volume you want to dispense into the target well
            exceeds the volume that will fit in one tip, you must either specify
            `allow_carryover` to allow the tip to carry on pipetting from the
            source wells after it has touched the target well, or break up your
            operation into multiple groups with separate tips.
        
            Parameters
            ----------
            sources : Well, WellGroup
                Well or wells to transfer liquid from.
            dest : Well
                Well to which to transfer consolidated liquid.
            volumes : str, Unit, list
                The volume(s) of liquid to be transferred from source well(s) to
                destination well.  Volume can be specified as a single string or
                Unit, or can be given as a list of volumes.  The length of a list
                of volumes must match the number of source wells given.
            mix_after : bool, optional
                Specify whether to mix the liquid in the destination well after
                liquid is transferred.
            mix_vol : str, Unit, optional
                Volume to aspirate and dispense in order to mix liquid in a wells
                before and/or after each transfer step.
            repetitions : int, optional
                Number of times to aspirate and dispense in order to mix
                liquid in well before and/or after each transfer step.
            flowrate : str, Unit, optional
                Speed at which to mix liquid in well before and/or after each
                transfer step.
            aspirate speed : str, Unit, optional
                Speed at which to aspirate liquid from source well.  May not be
                specified if aspirate_source is also specified. By default this
                is the maximum aspiration speed, with the start speed being half
                of the speed specified.
            dispense_speed : str, Unit, optional
                Speed at which to dispense liquid into the destination well. May
                not be specified if dispense_target is also specified.
            aspirate_source : fn, optional
                Options for aspirating liquid. Cannot be specified if
                aspirate_speed is also specified.
            dispense_target : fn, optional
                Options for dispensing liquid. Cannot be specified if
                dispense_speed is also specified.
            pre_buffer : str, Unit, optional
                Volume of air aspirated before aspirating liquid.
            transit_vol : str, Unit, optional
                Volume of air aspirated after aspirating liquid to reduce
                presence of bubbles at pipette tip.
            blowout_buffer : bool, optional
                If true the operation will dispense the pre_buffer along with the
                dispense volume cannot be true if disposal_vol is specified.
        
            Raises
            ------
            TypeError
                If supplying more than one destination well for consolidation.
            ValueError
                If a volume list is supplied and the length does not match the
                number of source wells.
            """    
            
            self._complete_mix_kwargs(mix_kwargs,sources,dest,volumes)
            
            #check that we have enough volume in source for this operation
            
            #check that we have enough space in desintation for this operation
            
            self._assert_safe_to_pipette_from(sources)
            
            if not isinstance(volumes, list) and volumes>MAX_PIPETTE_TIP_VOLUME:
                
                remaining_volume = volumes
                
                while remaining_volume>ul(0):
                    volume_to_consolidate = min(remaining_volume,MAX_PIPETTE_TIP_VOLUME)
                    
                    super(CustomProtocol,self).consolidate(sources, dest, volume_to_consolidate, 
                                                           allow_carryover=allow_carryover,
                                                           aspirate_speed=aspirate_speed, dispense_speed=dispense_speed, 
                                                           aspirate_source=aspirate_source,
                                                           dispense_target=dispense_target, 
                                                           pre_buffer=pre_buffer, transit_vol=transit_vol,
                                                           blowout_buffer=blowout_buffer, 
                                                           tip_type=tip_type, new_group=new_group,
                                                           **mix_kwargs)    
                    
                    remaining_volume-=volume_to_consolidate
                
                
            else:
                
            
                super(CustomProtocol,self).consolidate(sources, dest, volumes, 
                                                       allow_carryover=allow_carryover,
                                                       aspirate_speed=aspirate_speed, dispense_speed=dispense_speed, 
                                                       aspirate_source=aspirate_source,
                                                       dispense_target=dispense_target, 
                                                       pre_buffer=pre_buffer, transit_vol=transit_vol,
                                                       blowout_buffer=blowout_buffer, 
                                                       tip_type=tip_type, new_group=new_group,
                                                       **mix_kwargs)
                
            self._assert_valid_transfer(sources, dest)
                
            
                
    def image_plate(self, ref, mode, dataref):
        """
        Cover and Capture an image of the specified container.
    
        Example Usage:
    
        .. code-block:: python
    
            p = Protocol()
    
            agar_plate = p.ref("agar_plate", None, "1-flat", discard=True)
            bact = p.ref("bacteria", None, "micro-1.5", discard=True)
    
            p.spread(bact.well(0), agar_plate.well(0), "55:microliter")
            p.incubate(agar_plate, "warm_37", "18:hour")
            p.image_plate(agar_plate, mode="top", dataref="my_plate_image_1")
    
    
        Autoprotocol Output:
    
        .. code-block:: json
    
            {
              "refs": {
                "bacteria": {
                  "new": "micro-1.5",
                  "discard": true
                },
                "agar_plate": {
                  "new": "1-flat",
                  "discard": true
                }
              },
              "instructions": [
                {
                  "volume": "55.0:microliter",
                  "to": "agar_plate/0",
                  "from": "bacteria/0",
                  "op": "spread"
                },
                {
                  "where": "warm_37",
                  "object": "agar_plate",
                  "co2_percent": 0,
                  "duration": "18:hour",
                  "shaking": false,
                  "op": "incubate"
                },
                {
                  "dataref": "my_plate_image_1",
                  "object": "agar_plate",
                  "mode": "top",
                  "op": "image_plate"
                }
              ]
            }
    
    
        Parameters
        ----------
        ref : str, Container
            Container to take image of
        mode : str
            Imaging mode (currently supported: "top")
        dataref : str
            Name of data reference of resulting image
    
        """
        
        self.cover(ref)
        super(CustomProtocol,self).image_plate(ref, mode, dataref)
        
    def assert_valid_state(self):
        for n, ref in self.refs.items():
            for well in ref.container.all_wells():
                assert_non_negative_well(well)        
        
        
              
        
    def _get_final_protocol(self):

        #if everything is provisioned, then this protocol is the final protocol
        if all([volume==ul(0) for volume in self.unprovisioned_stock_reagent_volumes.values()]):
            return self
        
        #we have unprovisioned reagent
        
        p = CustomProtocol(parent_protocol=self)
        
        
        for reagent,volume_needed in self.unprovisioned_stock_reagent_volumes.items():
            reagent_info = self.our_inventory[reagent]
            #convert Reagent.culture_medium to culture_medium
            stock_plate = p.ref(reagent.name, 
                                cont_type=reagent_info['cont_type'],
                                discard=True)
            reagent_info['refill_function'](p, stock_plate,volume=volume_needed)
            
        #append our instructions
        p.instructions+=self.instructions
        setattr(p, "time_constraints", (getattr(
            p, "time_constraints", []) + \
                                    getattr(self, "time_constraints", [])))
        
        #build outs for each
        super(CustomProtocol,self).as_dict()
        super(CustomProtocol,p).as_dict()
        
        other_outs = getattr(p, "outs", {})
        my_outs = getattr(self, "outs", {})
        
        other_outs.update(my_outs)        

        if other_outs:
            setattr(p, "outs", other_outs)    
        
        #append our refs (except the reagent plates we just provisioned)
        for ref_name,ref in self.refs.items():
            if ref_name not in [reagent.name for reagent in self.unprovisioned_stock_reagent_volumes.keys()]:
                p.refs[ref_name] = ref
             
        return p
    
    
    #We have modified this to look for containers in the parent protocol in case
    #our instrunctions came from it
    def _ref_for_container(self, container):
        for k in self.refs:
            v = self.refs[k]
            if v.container is container:
                return k    
            
        if self.parent_protocol:
            return self.parent_protocol._ref_for_container(container)
    
    def as_dict(self, finalize=True):
        """
        Return the entire protocol as a dictionary.
    
        Example Usage:
    
        .. code-block:: python
    
            from autoprotocol.protocol import Protocol
            import json
    
            p = Protocol()
            sample_ref_2 = p.ref("sample_plate_2",
                                  id="ct1cxae33lkj",
                                  cont_type="96-pcr",
                                  storage="ambient")
            p.seal(sample_ref_2)
            p.incubate(sample_ref_2, "warm_37", "20:minute")
    
            print json.dumps(p.as_dict(), indent=2)
    
        Autoprotocol Output:
    
        .. code-block:: json
    
            {
              "refs": {
                "sample_plate_2": {
                  "id": "ct1cxae33lkj",
                  "store": {
                    "where": "ambient"
                  }
                }
              },
              "instructions": [
                {
                  "object": "sample_plate_2",
                  "op": "seal"
                },
                {
                  "duration": "20:minute",
                  "where": "warm_37",
                  "object": "sample_plate_2",
                  "shaking": false,
                  "op": "incubate"
                }
              ]
            }
    
        Returns
        -------
        dict
            dict with keys "refs" and "instructions" and optionally "time_constraints" and "outs",
            each of which contain the "refified" contents of their corresponding Protocol attribute.
    
        """
        
        if finalize:
            p = self._get_final_protocol()
            return p.as_dict(finalize=False)
        
        self.assert_valid_state()
                
        return super(CustomProtocol,self).as_dict()
    
    
    def spin(self, ref, acceleration, duration, flow_direction=None, spin_direction=None):
        """
        Apply acceleration to a container.
    
        Example Usage:
    
        .. code-block:: python
    
            p = Protocol()
            sample_plate = p.ref("sample_plate",
                                 None,
                                 "96-flat",
                                 storage="warm_37")
    
            p.spin(sample_plate, "1000:g", "20:minute", flow_direction="outward)
    
        Autoprotocol Output:
    
        .. code-block:: json
    
            "instructions": [
                {
                  "acceleration": "1000:g",
                  "duration": "20:minute",
                  "flow_direction": "outward",
                  "spin_direction": [
                    "cw",
                    "ccw"
                  ]
                  "object": "sample_plate",
                  "op": "spin"
                }
            ]
    
        Parameters
        ----------
        ref : Container
            The container to be centrifuged.
        acceleration: str
            Acceleration to be applied to the plate, in units of `g` or
            `meter/second^2`.
        duration: str, Unit
            Length of time that acceleration should be applied.
        flow_direction: str
            Specifies the direction contents will tend toward with respect to
            the container. Valid directions are "inward" and "outward", default
            value is "inward".
        spin_direction: list of strings
            A list of "cw" (clockwise), "cww" (counterclockwise). For each
            element in the list, the container will be spun in the stated
            direction for the set "acceleration" and "duration". Default values
            are derived from the "flow_direction" parameter. If
            "flow_direction" is "outward", then "spin_direction" defaults to
            ["cw", "ccw"]. If "flow_direction" is "inward", then
            "spin_direction" defaults to ["cw"].
    
        Raises
        ------
        TypeError:
            If ref to spin is not of type Container.
        TypeError:
            If spin_direction or flow_direction are not properly formatted.
        ValueError:
            If spin_direction or flow_direction do not have appropriate values.
    
        """
        
        super(CustomProtocol,self).spin(ref, acceleration, duration, flow_direction, spin_direction)
        
        #fix for bug https://github.com/autoprotocol/autoprotocol-python/issues/140
        
        if flow_direction=='outward':
            for well in ref.all_wells():
                well.volume = get_well_dead_volume(well)
                
            #sanity check to make sure that we uncovered before spinning
            assert not ref.cover
            
    
    #we modified this to default to mirroring because humans think of time constraints as mirrored by default
    #without mirroring, the "to" event can happen anytime before the "from" event 
    #(its only considering the time after "from")
    def add_time_constraint(self, from_dict, to_dict, time_between, mirror=True):
        """
        Add time constraints from from_dict to to_dict. Time constraints
        ensure that the time from the from_dict to the to_dict does not exceed
        the time_between given. Care should be taken when applying time
        constraints as constraints may make some protocols impossible to
        schedule or run.
    
        Though autoprotocol orders instructions in a list, instructions do
        not need to be run in the order they are listed and instead depend on
        the preceding dependencies. Time constraints should be added with such
        limitations in mind.
    
        Constraints are directional; use mirror=True if the time constraint
        should be added in both directions.
    
        Example Usage:
    
        .. code-block:: python
    
            plate_1 = protocol.ref("plate_1", id=None, cont_type="96-flat",
                                   discard=True)
            plate_2 = protocol.ref("plate_2", id=None, cont_type="96-flat",
                                   discard=True)
    
            protocol.cover(plate_1)
            time_point_1 = protocol.get_instruction_index()
    
            protocol.cover(plate_2)
            time_point_2 = protocol.get_instruction_index()
    
            protocol.add_time_constraint(
                {"mark": plate_1, "state": "start"},
                {"mark": time_point_1, "state": "end"},
                "1:minute")
            protocol.add_time_constraint(
                {"mark": time_point_2, "state": "start"},
                {"mark": time_point_1, "state": "start"},
                "1:minute", True)
    
    
        Autoprotocol Output:
    
        .. code-block:: json
    
            {
              "refs": {
                "plate_1": {
                  "new": "96-flat",
                  "discard": true
                },
                "plate_2": {
                  "new": "96-flat",
                  "discard": true
                }
              },
              "time_constraints": [
                {
                  "to": {
                    "instruction_end": 0
                  },
                  "less_than": "1.0:minute",
                  "from": {
                    "ref_start": "plate_1"
                  }
                },
                {
                  "to": {
                    "instruction_start": 0
                  },
                  "less_than": "1.0:minute",
                  "from": {
                    "instruction_start": 1
                  }
                },
                {
                  "to": {
                    "instruction_start": 1
                  },
                  "less_than": "1.0:minute",
                  "from": {
                    "instruction_start": 0
                  }
                }
              ],
              "instructions": [
                {
                  "lid": "standard",
                  "object": "plate_1",
                  "op": "cover"
                },
                {
                  "lid": "standard",
                  "object": "plate_2",
                  "op": "cover"
                }
              ]
            }
    
    
        Parameters
        ----------
        from_dict: dict
            Dictionary defining the initial time constraint condition.
            Composed of keys: "mark" and "state"
    
            mark: int or Container
                instruction index of container
            state: "start" or "end"
                specifies either the start or end of the "mark" point
    
        to_dict: dict
            Dictionary defining the end time constraint condition.
            Specified in the same format as from_dict
        time_between : str, Unit
            max time between from_dict and to_dict
        mirror: bool, optional
            choice to mirror the from and to positions when time constraints
            should be added in both directions
    
        Raises
        ------
        ValueError
            If an instruction mark is less than 0
        TypeError
            If mark is not container or integer
        TypeError
            If state not in ['start', 'end']
        KeyError
            If to_dict or from_dict does not contain 'mark'
        KeyError
            If to_dict or from_dict does not contain 'state'
        ValueError
            If time is less than '0:second'
        RuntimeError
            If from_dict and to_dict are equal
        RuntimeError
            If from_dict["marker"] and to_dict["marker"] are equal and
            from_dict["state"] = "end"
    
        """    
        
        super(CustomProtocol,self).add_time_constraint(from_dict, to_dict, 
                                                       time_between, mirror)
        
    
    def dispense(self, ref, reagent, columns, speed_percentage=None, is_resource_id=False):
        """
        Dispense specified reagent to specified columns.
    
        Example Usage:
    
        .. code-block:: python
    
            p = Protocol()
            sample_plate = p.ref("sample_plate",
                                 None,
                                 "96-flat",
                                 storage="warm_37")
    
            p.dispense(sample_plate,
                       "water",
                       [{"column": 0, "volume": "10:microliter"},
                        {"column": 1, "volume": "20:microliter"},
                        {"column": 2, "volume": "30:microliter"},
                        {"column": 3, "volume": "40:microliter"},
                        {"column": 4, "volume": "50:microliter"},
                        {"column": 5, "volume": "60:microliter"},
                        {"column": 6, "volume": "70:microliter"},
                        {"column": 7, "volume": "80:microliter"},
                        {"column": 8, "volume": "90:microliter"},
                        {"column": 9, "volume": "100:microliter"},
                        {"column": 10, "volume": "110:microliter"},
                        {"column": 11, "volume": "120:microliter"}
                       ])
    
    
        Autoprotocol Output:
    
        .. code-block:: json
    
            "instructions": [
                {
                  "reagent": "water",
                  "object": "sample_plate",
                  "columns": [
                    {
                      "column": 0,
                      "volume": "10:microliter"
                    },
                    {
                      "column": 1,
                      "volume": "20:microliter"
                    },
                    {
                      "column": 2,
                      "volume": "30:microliter"
                    },
                    {
                      "column": 3,
                      "volume": "40:microliter"
                    },
                    {
                      "column": 4,
                      "volume": "50:microliter"
                    },
                    {
                      "column": 5,
                      "volume": "60:microliter"
                    },
                    {
                      "column": 6,
                      "volume": "70:microliter"
                    },
                    {
                      "column": 7,
                      "volume": "80:microliter"
                    },
                    {
                      "column": 8,
                      "volume": "90:microliter"
                    },
                    {
                      "column": 9,
                      "volume": "100:microliter"
                    },
                    {
                      "column": 10,
                      "volume": "110:microliter"
                    },
                    {
                      "column": 11,
                      "volume": "120:microliter"
                    }
                  ],
                  "op": "dispense"
                }
              ]
    
        Parameters
        ----------
        ref : Container
            Container for reagent to be dispensed to.
        reagent : str
            Reagent to be dispensed to columns in container.
        columns : list
            Columns to be dispensed to, in the form of a list of dicts specifying
            the column number and the volume to be dispensed to that column.
            Columns are expressed as integers indexed from 0.
            [{"column": <column num>, "volume": <volume>}, ...]
        speed_percentage : int, optional
            Integer between 1 and 100 that represents the percentage of the
            maximum speed at which liquid is dispensed from the reagent
            dispenser.
        is_resource_id : bool, optional
            If true, interprets reagent as a resource ID
    
        """
        
        super(CustomProtocol,self).dispense(ref, reagent, columns, 
                                            speed_percentage=speed_percentage, is_resource_id=is_resource_id)  
        
        
        dest_wells = set()
        for col_info in columns:
            dest_wells.update(get_column_wells(ref, col_info['column']))
        
        self._assert_valid_add_volume(dest_wells)
          
        
    #--------------------------------------------
    #--------------- Bacterial Methods ----------
    #--------------------------------------------
    
    def add_antibiotic(p,wellsorcontainer,antibiotic, broth_volume=None,
                       mix_after=None,**mix_kwargs):
        """
    
        Ensures that all wells provided have the appropriate concentration of the given antibiotic.
        ignores empty wells
    
        """
        assert isinstance(antibiotic, Antibiotic)
        assert isinstance(p,Protocol)
        wells = convert_to_wellgroup(wellsorcontainer)
    
    
        
    
        if broth_volume and antibiotic.broth:
            
            #user didn't set a mix_after argument and reagents are pre-mixed
            if mix_after==None:
                mix_after=False
            
            p.provision_by_name(antibiotic.broth,wells,broth_volume, mix_after=mix_after, **mix_kwargs)
            return
    
    
        #we are mixing antibiotics here, so default to mix after
        if mix_after==None:
            mix_after=True    
    
        if broth_volume:
            p.provision_by_name(Reagent.lb_miller,wells,broth_volume)
    
        volumes = []
    
        for well in wells:
            #reagent_concentration is now in microliter/microgram
            reagent_concentration = 1.0/antibiotic.reagent_concentration.to('microgram/microliter')
            effective_concentration = antibiotic.effective_concentration.to('microgram/microliter')
            antibiotic_volume_to_add = ceil_volume(effective_concentration * well.volume.to('microliter') * reagent_concentration,1)
    
            volumes.append(antibiotic_volume_to_add)
    
        p.provision_by_name(antibiotic.reagent,wells,volumes, mix_after=mix_after, **mix_kwargs)    
    
    def measure_bacterial_density(self,wells, ):
            
        wells = convert_to_wellgroup(wells)
        
        source_container = wells[0].container
        
        num_rows = source_container.container_type.row_count()
        
        assert isinstance(source_container,Container)
        
        transfer_to_absorbance_plate = False
        if source_container.container_type.shortname != '96-flat':
            transfer_to_absorbance_plate = True
    
        if transfer_to_absorbance_plate:    
    
            if not self.absorbance_plate:
                self.absorbance_plate = self.ref('absorbance_plate', cont_type="96-flat", discard=True)        
            
            if not source_container.container_type.shortname.startswith('96-'):
                raise Exception("OD measurement only works for 96-well source, please update")
            
            
            if self.next_absorbance_plate_index>=12:
                raise Exception('code only knows how to run 12 absorbance calls, please update to increase this')
            
            column_index = wells[0].index
            
            if set(get_column_wells(source_container,column_index))!=set(wells):
                raise Exception('Entire columns allowed only if transfer required for OD measurement')
            
            self.transfer_column(source_container, column_index, 
                                self.absorbance_plate, 
                                self.next_absorbance_plate_index,
                                ul(100))
            
            wells = get_column_wells(self.absorbance_plate,
                                     self.next_absorbance_plate_index)
            
            self.next_absorbance_plate_index+=1
            
        
        self.absorbance(wells[0].container, wells.indices(),
                        wavelength="600:nanometer",
                        dataref='cell_density_600nm_absorbance_%s'%self.absorbance_measurement_count, flashes=25)    
        
        self.absorbance_measurement_count+=1
    
    
    def create_agar_plate(self, name,container_type_name,antibiotic=None,discard=True,storage=None):
        
        plates = {'6-flat':
                  {Antibiotic.kan: "ki17rs7j799zc2",
                      Antibiotic.amp: "ki17sbb845ssx9",
                      Antibiotic.spc: "ki17sbb9r7jf98",
                      Antibiotic.cam: "ki17urn3gg8tmj",
                      None: "ki17reefwqq3sq"},
                  '1-flat':{Antibiotic.kan: "ki17t8j7kkzc4g",
                      Antibiotic.amp: "ki17t8jcebshtr",
                      Antibiotic.spc: "ki17t8jaa96pw3",
                      Antibiotic.cam: "ki17urn592xejq",
                      None: "ki17t8jejbea4z"}     
                  }
        
        try:
            kit_id = plates[container_type_name][antibiotic]
        except KeyError,e:
            raise Exception('we can\'t find a plate for %s with antibiotic %s'%(container_type_name, antibiotic.name))
        
        return self.ref_kit_container(name, container_type_name, kit_id, 
                                     discard, storage)
        
    
    def ref_kit_container(self, name, container_type_name, kit_id, discard=True, storage=None):
        kit_item = Container(None, self.container_type(container_type_name))
        if storage:
            self.refs[name] = Ref(name, {"reserve": kit_id, "store": {"where": storage}}, kit_item)
        else:
            self.refs[name] = Ref(name, {"reserve": kit_id, "discard": discard}, kit_item)
        return kit_item    
    
    def miniprep(self, sources, dests):
        """
        
        Run a miniprep to isolate plasmids and genomic DNA.
        
        Removes x volume
        
        To well will have 45uL of solution
        
        """
        
  
        sources = convert_to_wellgroup(sources)
        dests = convert_to_wellgroup(dests)
        
    
        #source wells must have sufficient volume (1ml) although 1.8mL is recomended
        #it might be ideal to automatically move part of the sample over to a 96-deep and incubate if the volume is too low 
        #[we would still want to run an oD600 check before minipreping ]
  
        for source_well in sources:
            if get_volume(source_well,aspiratable=True) < ml(1):
                raise Exception('At least 1ml of aspirtable volume is needed'+\
                                'to mini prep but %s in well %s'%(get_volume(source_well,aspiratable=True),
                                                                 source_well))
    
    
        if len(sources) != len(dests):
            raise RuntimeError("Must have the same number of source and destination wells for miniprep") 

        groups = []
        for source, dest in zip(sources,dests):
            groups.append({
                "from": source,
                "to": dest
            })
            
            
        self.instructions.append(MiniPrep(groups))
        
        for i in range(0,len(sources)): 
            dests[i].name = '%s_miniprep'%sources[i].name   
            
            #minipreping uses all available volume
            sources[i].volume = get_well_dead_volume(sources[i])
        
        
    def pcr(self, template_well, primer1_well, primer2_well,
            annealing_temp_c, greater_than_65_percent_gc_primer,
            product_length,
            product_name='pcr product',
            negative_control=True,
            discard_pcr_plate=False,
            run_gel=True,
            measure_concentration=True):
        """
        
        #use the following calculator to calc Tm for primers - http://tmcalculator.neb.com/#!/ (w/ 500nM primer)
        
        Assumes primer well concentration is 10uM
        
        Bands in the negative control indicate that there was contamination in the primers'
        
        #primer Tm should be between 55-60C and within 2C of each other
        
        """
        
        if product_length>20*1000:
            raise Exception('Q5 polymerase is not appropriate for PCR products more than 20kb')
        
        
        #annealing_temp = min(primer1_tm_c, primer2_tm_c) + 1
        
        #if max(primer1_tm_c, primer2_tm_c)>72:
            #raise Exception('primer melting temps are far too high, aim for between 55-60')
        
        
        #determine if we need GC enhancer
        
        #dilute primers to 10uM
        
       
        
        # Temporary tubes for use, then discarded (you can't set storage if you are going to discard)
        mastermix_well = self.ref("mastermix", cont_type="micro-1.5", discard=True).well(0)
        pcr_plate =      self.ref("%s_pcr_plate"%product_name, cont_type="96-pcr", storage=Temperature.cold_4 if not discard_pcr_plate else None,
                                  discard = discard_pcr_plate)
        
        
        experimental_pcr_well = pcr_plate.wells(["A1"])[0]
        experimental_pcr_well.name = product_name     
        experimental_pcr_well.properties['dna_length'] = str(product_length)
        
        pcr_wells = [experimental_pcr_well]
        
        if negative_control:
            negative_control_well = pcr_plate.wells(["A2"])[0]
            negative_control_well.name = 'negative_control_no_template'
            pcr_wells.append(negative_control_well)

        # Initialize all existing inventory
        all_inventory_wells = [template_well, primer1_well, primer2_well]
        for well in all_inventory_wells:
            init_inventory_well(well)
    
        # -----------------------------------------------------
        # Q5 PCR protocol
        # www.neb.com/protocols/2013/12/13/pcr-using-q5-high-fidelity-dna-polymerase-m0491
        #
        # 25ul reaction (we will run it 2 times, totally 50uL of solution)
        # -------------                      3rxn
        # Q5 reaction buffer      5    ul --> 15uL
        # Q5 polymerase           0.25 ul --> 0.75uL
        # 10mM dNTP               0.5  ul --> 1.5uL
        #if >65% GC
        # GC enhancer             5    ul --> 15uL
        # 10uM forward primer     1.25 ul --> 3.75
        # 10uM reverse primer     1.25 ul --> 3.75
        # you only need 12k copies of a target, ultimately we need moles of target but this is a rough estimate for the DNA we use
        # 1ng Template (1 ng/ul) 1 ul --> 1uL (1ng/ul concentration) 
        # water                   
        # -------------------------------
        # Sum                     14.25 ul --> 42.75
        # water                   10.75 uL --> 32.25
        #if <65% GC, add 5uL to water
        #
        #
        
        water_ul = ul(32.25 + (0 if greater_than_65_percent_gc_primer else 15))
        
        # Mastermix tube will have 96ul of stuff, leaving space for 3x1ul aliquots of template
        self.provision_by_name(Reagent.water, mastermix_well, water_ul)
        self.provision_by_name(Reagent.q5_buffer_5x, mastermix_well, ul(15))
        #provisioning is to nearest .1
        self.provision_by_name(Reagent.q5_hf_polymerase, mastermix_well, ul(.8))
        self.provision_by_name(Reagent.dntp_10mM, mastermix_well, ul(1))
        if greater_than_65_percent_gc_primer:
            self.provision_by_name(Reagent.q5_high_gc_enhancer, mastermix_well, ul(15))
        self.transfer(primer1_well, mastermix_well, ul(3.75), mix_before=True, mix_after=True)
        self.transfer(primer2_well, mastermix_well, ul(3.75), mix_before=True, mix_after=True)
        
        # Transfer mastermix to pcr_plate without template
        self.transfer(mastermix_well, pcr_wells, ul(24))
        
        # Finally add template
        self.transfer(template_well,  experimental_pcr_well, ul(1), mix_after=True)
        
        if negative_control:
            self.provision_by_name(Reagent.water, negative_control_well, ul(1))
        
        # ---------------------------------------------------------
        # Thermocycle with Q5 and hot start
        # min Tm + 1 annealing temperature is recommended by NEB protocol
        # self.seal is enforced by transcriptic
        #
        
        seconds_per_kb = 30
        
        if product_length>=6*1000:
            seconds_per_kb = 50
        
        extension_time = max(10, product_length*seconds_per_kb/1000)
        
        
        cycles = [{"cycles":  1, "steps": [{"temperature": "98:celsius", "duration": "30:second"}]}] + \
            touchdown_pcr(74.8, annealing_temp_c, [8, 25, extension_time], stepsize=1) +\
            [{"cycles": 20, "steps": [{"temperature": "98:celsius", "duration": "8:second"},
                                      {"temperature": "%s:celsius"%annealing_temp_c, "duration": "25:second"},
                                      {"temperature": "72:celsius", "duration": "{:d}:second".format(extension_time)}]},
             {"cycles":  1, "steps": [{"temperature": "72:celsius", "duration": "2:minute"}]}]
        self.seal(pcr_plate)
        self.thermocycle(pcr_plate, cycles, volume=ul(25))
        self.unseal(pcr_plate)

        if run_gel:
            
            if product_length <= 2*1000:
                ladder = "ladder1"
                duration = "10:minute"
            else:
                ladder = "ladder2"   
                duration = "15:minute"
            
            self.mix(pcr_wells)
            
            #2% for 10min can handle 2kb - 100bp
            self.gel_separate(pcr_wells,
                           ul(10), "agarose(10,2%)", ladder, duration, 'pcr_gel')
        
        
        #return the experimental pcr well
        return experimental_pcr_well
    
    
    def simple_gel_purify(self, source_wells, dna_lengths, dest_wells=None, dataref='gel_purify'):
        """
        
        A simple version of gel purify.
        
        min(all,25uL) DNA inside source wells will be purified and moved to dest_wells.
        
        If dest_wells is blank, a pcr plate will be created and a well assigned as a dest well for each source well.
        
        This protocol is intended mostly for PCR cleanup where you are extracting a single band 
        
        """
        
        source_wells = convert_to_wellgroup(source_wells)
        if dest_wells:
            dest_wells = convert_to_wellgroup(dest_wells)
        else:
            dest_wells = []
            
            dest_plate = self.ref("purified_dna_plate", None,
                                  "96-pcr", storage="cold_4")      
            next_well_index = 0
            
            for i, source_well in enumerate(source_wells):
                if source_well.name:
                    source_name = source_well.name
                elif source_well.container.name:
                    source_name = source_well.container.name
                else:
                    source_name = 'gene_product_%s'%i
                    
                dest_well = dest_plate.well(next_well_index)
                
                dest_well.name = source_name
                
                set_property(dest_well, 'dna_length', dna_lengths[i])
                
                next_well_index+=1
                
                dest_wells.append(dest_well)
            
        
        if isinstance(dna_lengths, int):
            dna_lengths = [dna_lengths]*len(source_wells)
        
    
        assert isinstance(dna_lengths,list)
        
        min_well_volume = min([get_volume(well, aspiratable=True) for well in source_wells])
        
        volume = min(ul(25),min_well_volume)
        
        if volume<ul(5):
            raise Exception('must have at least 5uL in source available for gel extract')
        
        if max(dna_lengths) <= 2*1000:
            ladder = "ladder1"
        else:
            ladder = "ladder2"
            
        if max(dna_lengths) > 10*1000:
            raise Exception('current gels don\'t support more than 10kb')
        
        if min(dna_lengths) <50 :
            raise Exception('current gels don\'t support less than 50bp')        
        
        
        #each item in this list represents extractions from a single lane
        extracts = [make_gel_extract_params(
            w,
            make_band_param(
                "water",
                "25:microliter",
                #we take a wide window since bands can smear
                min(10*1000,dna_lengths[i]+50),
                max(50,dna_lengths[i]-450),
                dest_wells[i]))
                    for i, w in enumerate(source_wells)]
        
        self.mix(source_wells)
        
        self.gel_purify(extracts, volume,
                        "size_select(8,2%)", 
                        ladder,
                        dataref=dataref)      
        
        return dest_wells
        
        
    def linearize_dna(self, vector_source_well,
                      negative_control=True,
                      gel_verify=True):        
        """
        
        Cuts a vector with hindiii and sali.  
        @TODO: update this method to take enzymes and dynamically generate protocol from http://nebcloner.neb.com/#!/redigest
        
        """
        pcr_plate  = self.ref("linearized_%s"%vector_source_well.name,  cont_type="96-pcr", storage=Temperature.cold_4)
        experiment_well = pcr_plate.well(0)        
        experiment_well.name = "linearized_%s"%vector_source_well.name
        
        pcr_wells = [experiment_well]
        
        if negative_control:
            negative_control_well = pcr_plate.well(1)
            negative_control_well.name = 'uncut_%s'%vector_source_well.name
            pcr_wells.append(negative_control_well)
        
        
        # -------------------------------------------------------------
        # Restriction enzyme cutting vector
        # 2 experiments (1 without re)
        # 50ul total reaction volume for cutting 1ug of DNA:
        # yuL water
        # 5ul CutSmart 10x
        # xul vector (1ug of DNA)
        # 1ul of each enzyme 
        # 1ul fastap (dephosphorylation)
        
        cutsmart_volume = ul(5)
        enzyme_volume = ul(1)
        antarctic_phos_volume = ul(1) # 5 units
        antarctic_phos_buffer_volume =  ul(5)
        vector_volume = convert_mass_to_volume(ug(1), vector_source_well)
        water_volume = ul(50) - (2*enzyme_volume + vector_volume + cutsmart_volume + antarctic_phos_volume +\
                                 antarctic_phos_buffer_volume)
        
        self.provision_by_name(Reagent.water, pcr_wells, water_volume)
        self.mix(vector_source_well)
        self.transfer(vector_source_well, pcr_wells, vector_volume, mix_after=True,
                      one_source=True, one_tip=True)
        self.provision_by_name(Reagent.cutsmart_buffer_10x, pcr_wells, cutsmart_volume)
        if negative_control:
            self.provision_by_name(Reagent.water, negative_control_well, 2*enzyme_volume+ \
                                   antarctic_phos_buffer_volume + antarctic_phos_volume)
        
        self.provision_by_name(Reagent.hindiii, experiment_well, enzyme_volume)
        self.provision_by_name(Reagent.sali, experiment_well, enzyme_volume)
        
        self.provision_by_name(Reagent.antarctic_phosphatase, experiment_well, antarctic_phos_volume)
        self.provision_by_name(Reagent.antarctic_phosphatase_buffer_10x, experiment_well, antarctic_phos_buffer_volume)
        
        
        assert all([well.volume == ul(50) for well in pcr_wells])
        
        self.mix(pcr_wells)
        
        #enzyme only takes 15min but we want to leave time for dephosphorylation
        self.incubate(pcr_plate, "warm_37", "30:minute", shaking=False)
        
        #inactivate RE's (see this site for thermal temps)
        #https://www.neb.com/tools-and-resources/usage-guidelines/heat-inactivation
        self.thermocycle(pcr_plate, [{"cycles":  1, "steps": [{"temperature": "80:celsius", "duration": "20:minute"}]}], volume=ul(50))
        
        if gel_verify:
            # --------------------------------------------------------------
            # Gel electrophoresis, to ensure the cutting worked
            #
            self.mix(pcr_wells)
            self.gel_separate(pcr_wells, ul(10), "agarose(10,2%)", 'ladder2', '15:minute', 'verify_linearization_gel')            
                
        return experiment_well