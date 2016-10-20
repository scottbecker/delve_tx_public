from __future__ import print_function
import json
contents = open('../src/manifest.json','rb').read().decode('utf8')
manifest = json.loads(contents)

for protocol in manifest['protocols']:
    if protocol.get('preview'):
        protocol['preview'] = {}
      
f = open('/tmp/manifest.json','wb')
f.write(json.dumps(manifest,indent=True))
