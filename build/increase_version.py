import json
contents = open('../src/manifest.json','rb').read().decode('utf8')
manifest = json.loads(contents)

for protocol in manifest['protocols']:
   major_version, minor_version, patch = protocol['version'].split('.')
   patch=int(patch)+1
   new_version = '%s.%s.%s'%(major_version,minor_version,patch)
   protocol['version'] = new_version

f = open('../src/manifest.json','wb')
f.write(json.dumps(manifest,indent=4))