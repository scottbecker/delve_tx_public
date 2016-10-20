import json
import os

manifest_contents = open('../src/manifest.json','rb').read().decode('utf8')
manifest = json.loads(manifest_contents)

public_protocols = json.loads(open('public_protocols.json','rb').read().decode('utf8'))


for protocol in manifest['protocols']:
   if protocol['name'] not in public_protocols:
      continue
   
   protocol_path = protocol['command_string'].split(' ')[2].replace('.','/').replace('/protocol','')
   os.system('rsync -av --delete ../src/%s ../../dtx_public/src/protocols/ --exclude *.pyc'%protocol_path)

   
