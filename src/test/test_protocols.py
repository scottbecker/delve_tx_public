from __future__ import print_function
import unittest
import json
import importlib
import os
import tempfile
from transcriptic_tools.harness import run

class TestProtocols(unittest.TestCase):
    longMessage = True
    
    def setUp(self):

        self.temp_file = tempfile.NamedTemporaryFile(delete=False)
        self.temp_file_path = self.temp_file.name
        self.temp_file.close()
        
    def tearDown(self):
        os.unlink(self.temp_file_path)

def make_test_function(protocol_name,inputs,module_path):
    def test(self):
    
        module = importlib.import_module(module_path)
        
        
        f = open(self.temp_file_path,'wb')
        f.write(json.dumps(inputs,indent=4))       
        f.close()
       
        
        run(module.main,protocol_name, config=self.temp_file_path)
    
    return test


contents = open('manifest.json','rb').read().decode('utf8')
manifest = json.loads(contents)

for protocol in manifest['protocols']:
    protocol_name = protocol['name']
    input_sets = [protocol['preview']]
    if 'additional_previews' in protocol:
        input_sets+= protocol['additional_previews']
    
    #pull out the module path
    module_path = protocol['command_string'].split(' ')[2]
    
    for i in range(0,len(input_sets)):
        inputs = input_sets[i]
        test_func = make_test_function(protocol_name,inputs,module_path)
        setattr(TestProtocols, 'test_{0}_{1}'.format(protocol_name,i), test_func)
    
    
