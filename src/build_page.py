import json,re
D='/tmp/claude-0/-home-claude/ee499b9d-9071-5623-bb08-56ecb36a99eb/scratchpad/endo/'
t=open(D+'template.html').read()
bb=[l for l in open('/mnt/project-files/namd-injections/namd-injections.html') if 'class="brandbar"' in l][0].strip()
data=json.load(open(D+'data.json'))
nq=sum(1 for l in open(D+'qlog.tsv') if l.strip())+sum(1 for l in open(D+'querylog.tsv') if l.strip())
nmeta=sum(1 for l in open(D+'meta_store.jsonl'))
t=t.replace('/*BRANDBAR*/',bb).replace('/*DATA*/',json.dumps(data,separators=(',',':'))).replace('/*QCOUNT*/',str(nq)).replace('/*META*/',f'{nmeta:,}')
open('/mnt/project-files/endophthalmitis/endophthalmitis.html','w').write(t)
print(len(t),nq,nmeta)
