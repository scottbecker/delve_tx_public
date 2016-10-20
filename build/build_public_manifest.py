import json

manifest_contents = open('../src/manifest.json','rb').read().decode('utf8')
manifest = json.loads(manifest_contents)

public_protocols = json.loads(open('public_protocols.json','rb').read().decode('utf8'))

public_manifest = {
   'format': 'python',
   'license': 'Proprietary',
   'protocols': []
}

for protocol in manifest['protocols']:
   if protocol['name'] not in public_protocols:
      continue
   public_manifest['protocols'].append(protocol)

f = open('../../dtx_public/src/manifest.json','wb')
f.write(json.dumps(public_manifest,indent=4))