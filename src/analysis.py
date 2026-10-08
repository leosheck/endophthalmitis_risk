# Builds data.json for the endophthalmitis page from the extracted records.
import json, math, re
from scipy.stats import beta, norm
D='/tmp/claude-0/-home-claude/ee499b9d-9071-5623-bb08-56ecb36a99eb/scratchpad/endo/'
ivi={r['key']:r for r in json.load(open(D+'ivi_raw.json'))}
cat={r['key']:r for r in json.load(open(D+'cat.json'))}

# ---- references -------------------------------------------------------
REF={}
_AF=json.load(open(D+'authfix.json'))
def addref(key, cite, pmid=None, doi=None, src='PubMed'):
    if cite.startswith('(authors not captured') and pmid in _AF:
        _au=_AF[pmid];cite=', '.join(_au[:3])+(' et al.' if len(_au)>3 else '.')+' '+cite[cite.index(')')+1:].strip()
    elif cite.startswith('(CATT Research Group)'):
        cite='CATT Research Group.'+cite[len('(CATT Research Group)'):]
    REF[key]=dict(key=key,cite=cite,pmid=pmid,doi=doi,src=src)
for k,r in ivi.items(): addref(k, r['citation'], r.get('pmid'), r.get('doi'))
for k,r in cat.items(): addref(k, r['citation'], r.get('pmid'), r.get('doi'))
# items held in Bear's library (wiki / Zotero) that are cited with their numbers
LIB={'Rohowetz2026','Israilevich2024','Tao2025PFS_NMA','Zhang2025PFS_MA','Uzzan2024','Ton2023','DilutePI2019','Shorstein2021KP','Friling2022NCR_ISBCS','Lacy2021IRIS_ISBCS','Malwankar2021ISBCS','Zafar2021Medicare','Gower2017Cochrane','Kessel2015','Huang2016','ESCRS2007','Barry2006ESCRSprelim','Tan2012Singapore','Garat2009Barcelona','RomeroAroca2012Reus'}
addref('Garat2005', 'Garat M, Moser CL, Alonso-Tarres C et al. Intracameral cefazolin to prevent endophthalmitis in cataract surgery: 3-year retrospective study. J Cataract Refract Surg. 2005.', None, '10.1016/j.jcrs.2005.08.053', 'Library')
for k in LIB:
    if k in REF: REF[k]['src']='PubMed and Bear\'s library'
used=[]
def cite(k):
    assert k in REF, k
    if k not in used: used.append(k)
    return k

# ---- statistics -------------------------------------------------------
def exact(x,n):
    lo = 0.0 if x==0 else beta.ppf(0.025,x,n-x+1)
    hi = 1.0 if x==n else beta.ppf(0.975,x+1,n-x)
    return x/n, lo, hi
def rr_counts(a,n1,b,n2):
    if a==0 or b==0:
        a+=.5;b+=.5;n1+=1;n2+=1
    rr=(a/n1)/(b/n2); se=math.sqrt(1/a-1/n1+1/b-1/n2)
    return rr, rr*math.exp(-1.96*se), rr*math.exp(1.96*se), se
def dl(ests):  # ests: list of (logest, se)
    k=len(ests)
    if k==0: return None
    w=[1/s**2 for _,s in ests]; sw=sum(w)
    fe=sum(wi*y for wi,(y,_) in zip(w,ests))/sw
    Q=sum(wi*(y-fe)**2 for wi,(y,_) in zip(w,ests))
    df=k-1
    C=sw-sum(wi**2 for wi in w)/sw
    tau2=max(0,(Q-df)/C) if df>0 else 0
    wr=[1/(s**2+tau2) for _,s in ests]; swr=sum(wr)
    mu=sum(wi*y for wi,(y,_) in zip(wr,ests))/swr; se=math.sqrt(1/swr)
    i2=max(0,(Q-df)/Q)*100 if Q>0 and df>0 else 0
    pi=None
    if k>=3:
        from scipy.stats import t
        tq=t.ppf(.975,k-2); h=tq*math.sqrt(tau2+se**2); pi=(mu-h,mu+h)
    return dict(mu=mu,lo=mu-1.96*se,hi=mu+1.96*se,k=k,i2=round(i2),tau2=tau2,pi=pi)
def logit(p): return math.log(p/(1-p))
def ilogit(x): return 1/(1+math.exp(-x))
def pool_rates(rows):
    e=[]
    for r in rows:
        x,n=r['cases'],r['n']
        if x==0: x,n=.5,n+1
        p=x/n; e.append((logit(p), math.sqrt(1/x+1/(n-x))))
    d=dl(e)
    out=dict(est=ilogit(d['mu']),lo=ilogit(d['lo']),hi=ilogit(d['hi']),k=d['k'],i2=d['i2'],
             cases=sum(r['cases'] for r in rows), n=sum(r['n'] for r in rows))
    if d['pi']: out['pilo'],out['pihi']=ilogit(d['pi'][0]),ilogit(d['pi'][1])
    return out

def parse_ci(s):
    m=re.findall(r'(\d+(?:\.\d+)?)',s); return float(m[0]),float(m[1])

# row builders ------------------------------------------------------------
def C(key,label,a,n1,b,n2,note=None,kind='study',pool=True,derived=False,group=None):
    rr,lo,hi,se=rr_counts(a,n1,b,n2)
    return dict(ref=cite(key),label=label,kind=kind,measure='RR',est=rr,lo=lo,hi=hi,se=se,
                a=f'{a:,}/{n1:,}',b=f'{b:,}/{n2:,}',ra=a/n1,rb=b/n2,pool=pool,note=note,derived=derived,group=group,
                counts=(a,n1,b,n2))
def R(key,label,est,lo,hi,measure,invert=False,note=None,kind='study',pool=True,group=None):
    if invert: est,lo,hi=1/est,1/hi,1/lo
    se=(math.log(hi)-math.log(lo))/3.92
    return dict(ref=cite(key),label=label,kind=kind,measure=measure+(' (inverted)' if invert else ''),est=est,lo=lo,hi=hi,se=se,
                a=None,b=None,pool=pool,note=note,group=group)
def M(key,label,est,lo,hi,measure,invert=False,note=None):
    return R(key,label,est,lo,hi,measure,invert,note,kind='meta',pool=False)

def factor(id,title,expo,ref,rows,note,pool=True,subpools=None,scale='ratio'):
    f=dict(id=id,title=title,expo=expo,refgroup=ref,rows=rows,note=note,scale=scale)
    if pool:
        pr=[r for r in rows if r['kind']=='study' and r['pool']]
        if len(pr)>=2:
            d=dl([(math.log(r['est']),r['se']) for r in pr])
            f['pooled']=dict(est=math.exp(d['mu']),lo=math.exp(d['lo']),hi=math.exp(d['hi']),k=d['k'],i2=d['i2'],
                             pilo=math.exp(d['pi'][0]) if d['pi'] else None,pihi=math.exp(d['pi'][1]) if d['pi'] else None)
        if subpools:
            f['subpools']=[]
            for g in subpools:
                pr=[r for r in rows if r['kind']=='study' and r['pool'] and r.get('group')==g]
                if len(pr)>=2:
                    d=dl([(math.log(r['est']),r['se']) for r in pr])
                    f['subpools'].append(dict(group=g,est=math.exp(d['mu']),lo=math.exp(d['lo']),hi=math.exp(d['hi']),k=d['k'],i2=d['i2']))
    for r in rows: r.pop('se',None); r.pop('counts',None)
    return f

# ===================== INJECTIONS: overall ===============================
EXCL={ # key: reason (overlap or not a per-injection series)
 'Rayess2016':'overlaps Wills Eye series [Israilevich 2024]','Rayess2015Diagnosis':'overlaps Wills Eye series','Pancholy2021':'overlaps Wills Eye series',
 'Pancholy2022DEX':'overlaps Wills Eye series','Patel2019WillsCulture':'overlaps Wills Eye series','Patel2020PhysMask':'overlaps Wills Eye series',
 'Patel2022Season':'overlaps Wills Eye series','Storey2014Abx':'overlaps Wills Eye series',
 'Gregori2015':'Bascom Palmer, replaced by 2005-2025 total','Moshfeghi2011':'Bascom Palmer, replaced by 2005-2025 total',
 'Yannuzzi2018':'Bascom Palmer, replaced by 2005-2025 total','ReyesCapo2021':'Bascom Palmer, replaced by 2005-2025 total','Rohowetz2026':'replaced by 2005-2025 total',
 'Baudin2018':'French national data, replaced by Baudin 2022 national database study','Baudin2022Abx':'French national data, replaced by Baudin 2022 national database study',
 'VanderBeek2015Bev':'same claims database as Vander Beek 2015','Chen2011Houston':'same practice as Li 2016',
 'Karimi2022Bev':'same centre as Karimi 2023','Tanaka2019':'same centre as Tanaka 2022','Morioka2020Abx':'same network as Morioka 2025',
 'IRIS2026FirstInjection':'first injection only','VictoriaDEX2026':'steroid implant only','Stem2017DEX':'steroid implant only',
}
SKIPDESIGN=re.compile(r'meta|pooled|systematic|outbreak|review|pharmacovigilance|literature|in vitro',re.I)
def grp(design):
    d=design.lower()
    if re.search(r'claims|registry|national|database|emr|medicare|cdm',d): return 'National, claims or registry'
    if re.search(r'multi|practices|centres|centers|hospitals|sites|departments|2 centres|two-centre|2 large|2 swiss|3 european|2 tertiary',d): return 'Multi-centre'
    return 'Single centre'
ser=[]; excluded=[]
for k,r in ivi.items():
    o=r.get('overall') or {}
    n,x=o.get('injections'),o.get('cases')
    if k in EXCL: excluded.append((k,EXCL[k])); continue
    if SKIPDESIGN.search(r['design'] or ''): continue
    if not n or x is None or n<10000: continue
    ser.append(dict(key=k,n=n,cases=x,design=r['design'],country=r.get('country'),period=r.get('period'),group=grp(r['design'])))
ser.append(dict(key='Rohowetz2026',n=374232,cases=44,design='consecutive series, all agents 2005-2025 (Bascom Palmer)',country='USA',period='2005-2025',group='Single centre',note='totals summed from per-drug counts'))
for s in ser:
    p,lo,hi=exact(s['cases'],s['n']); s.update(rate=p,lo=lo,hi=hi)
ser.sort(key=lambda s:-s['n'])
def ay(k):
    c=REF[k]['cite']; ys=re.findall(r'(?:19|20)\d{2}',c); return re.split(r'[ ,.]',c)[0]+(' '+ys[-1] if ys else '')
for s in ser: s['ref']=cite(s['key']); s['label']=ay(s['key'])+(', 2005-2025' if s['key']=='Rohowetz2026' else '')
ivi_pool=pool_rates(ser)
ivi_pool_grp={g:pool_rates([s for s in ser if s['group']==g]) for g in ['National, claims or registry','Multi-centre','Single centre']}
ivi_metas=[]
for k,lab in [('McCannel2011','US literature 2005-2010'),('Fileta2014','43 articles'),('Merani2015','20 large series'),('Reibaldi2018','60 articles'),('Menchini2018','13 studies')]:
    o=ivi[k]['overall']; p,lo,hi=exact(o['cases'],o['injections'])
    ivi_metas.append(dict(ref=cite(k),label=lab,rate=p,lo=lo,hi=hi,cases=o['cases'],n=o['injections']))
# cumulative per-eye observations
cumobs=[dict(ref=cite('Israilevich2024'),pts=[[11,0.20],[35,0.38],[66,0.46],[126,0.50]],label='Wills Eye cohort'),
        dict(ref=cite('Zhang2026IRIS'),pts=[[10,0.31],[60,0.58]],label='IRIS Registry')]
cumRef=cite('Bjerager2024Denmark')

# ===================== INJECTION FACTORS ===============================
IF=[]
IF.append(factor('pfs','Prefilled syringe','Prefilled syringe (manufacturer)','Drug drawn up from a vial',[
 C('Louis2023','Ranibizumab',17,112809,13,34790),
 C('Louis2023','Aflibercept 2 mg',7,29899,41,68076),
 C('Storey2019PFS','Ranibizumab',12,78407,43,165347),
 C('Feng2023','Ranibizumab',12,77925,21,93073),
 C('Dhoot2021Vestrum','Ranibizumab',96,427763,233,884061),
 C('Morioka2025','Aflibercept 2 mg',2,43684,10,108355),
 C('Rohowetz2026','Ranibizumab',6,31291,2,11689),
 C('Rohowetz2026','Aflibercept 2 mg',15,81460,9,87294),
 C('Rohowetz2026','Faricimab',1,3226,0,12366),
 C('David2026','Aflibercept 2 mg, 2023 vs 2022',15,35589,16,31991,note='Counts derived from reported rates (1 in 2,372 and 1 in 1,991).',derived=True),
 R('Baudin2018','Ranibizumab, adjusted',1.63,1.15,2.30,'IRR',invert=True,note='Reported as vial versus prefilled; inverted.'),
 M('Zhang2025PFS_MA','Meta-analysis, 11 cohorts',0.53,0.44,0.64,'RR'),
 M('Tao2025PFS_NMA','Network meta-analysis, 20 studies',0.45,0.40,0.49,'OR'),
], 'Comparisons within the same practice or database. Most are before-and-after or concurrent observational data, so era and practice policy changes (for example masks) run alongside the change in preparation.'))
IF.append(factor('pfs_cp','Prefilled syringe, culture-positive cases only','Prefilled syringe','Drug drawn up from a vial',[
 C('Storey2019PFS','Ranibizumab, culture-positive',2,78407,22,165347),
 M('Zhang2025PFS_MA','Meta-analysis, 3 studies',0.23,0.07,0.72,'RR'),
 M('Tao2025PFS_NMA','Network meta-analysis',0.15,0.06,0.41,'OR'),
], 'Restricting to culture-proven infection removes most sterile inflammation from the count.', pool=False))
IF.append(factor('comp','Pharmacy-compounded syringe','Pharmacy-compounded syringe','Drug drawn up from a vial',[
 C('Blom2020','Mixed agents',5,68150,6,44776),
 M('Zhang2025PFS_MA','Meta-analysis, 6 studies',0.52,0.38,0.72,'RR'),
 M('Tao2025PFS_NMA','Network meta-analysis',0.69,0.64,0.74,'OR'),
], 'Compounded bevacizumab was the source of the 2011 to 2013 US outbreaks (streptococcal and fungal), which are not part of these rates.', pool=False))
IF.append(factor('steroid','Steroid injection','Intravitreal steroid (triamcinolone or dexamethasone implant)','Anti-VEGF injection',[
 C('Dhoot2021Vestrum','Dexamethasone implant',53,49464,1395,4192245),
 C('Dhoot2021Vestrum','Triamcinolone',110,75038,1395,4192245,pool=False,note='Shares the anti-VEGF reference group with the row above; not pooled.'),
 C('Zabaneh2026','Dexamethasone implant',15,11954,106,318618),
 C('Pancholy2022DEX','Dexamethasone implant vs ranibizumab 0.5 mg',5,4973,43,163974),
 C('Louis2023','Dexamethasone implant or triamcinolone',8,4906,93,302443),
 C('Chen2025Taiwan','Triamcinolone',5,939,2,17426),
 C('Cheung2012','Triamcinolone vs ranibizumab',2,935,3,9453),
 R('VanderBeek2015','Steroids, adjusted',6.92,3.54,13.52,'OR'),
 R('Baudin2018','Steroids, adjusted',3.21,2.33,4.44,'IRR'),
 R('Emerson2021Medicare','Steroids, adjusted',3.91,2.75,5.56,'OR'),
 M('Tao2025PFS_NMA','Triamcinolone, network meta-analysis',3.44,3.03,3.89,'OR'),
 M('Tao2025PFS_NMA','Dexamethasone implant, network meta-analysis',2.71,2.38,3.10,'OR'),
], 'Implants are given through a larger-bore applicator. Some cases after triamcinolone are sterile reactions to the drug that are hard to tell apart from infection.'))
IF.append(factor('afl','Aflibercept','Aflibercept 2 mg','Bevacizumab or ranibizumab',[
 C('Dhoot2021Vestrum','vs bevacizumab',687,1412699,379,1467722),
 C('Rayess2016','vs bevacizumab',14,40356,60,153812),
 C('Louis2023','vs compounded bevacizumab',48,97975,15,56869),
 C('Feng2023','vs bevacizumab',45,122947,10,32045),
 C('Rohowetz2026','vs bevacizumab, 2005-2025',24,163071,11,142691),
 R('Kiss2018Claims','vs ranibizumab, adjusted',2.19,1.68,2.85,'OR'),
 R('Seong2022Korea','vs ranibizumab, adjusted',1.8,1.16,2.79,'IRR'),
 R('Bavinger2019','vs bevacizumab, adjusted',1.34,0.99,1.81,'OR'),
 R('Mishra2023CHX','vs bevacizumab',3.48,2.09,7.24,'OR'),
 R('Emerson2021Medicare','vs other agents, adjusted',1.47,1.19,1.80,'OR'),
], 'In most of these datasets aflibercept 2 mg was drawn from a vial while bevacizumab came as a pharmacy-filled syringe and ranibizumab was partly prefilled, so preparation and drug are mixed. Several studies also report more sterile inflammation with aflibercept, which some case definitions count.'))
IF.append(factor('bev','Bevacizumab','Bevacizumab (compounded)','Ranibizumab',[
 C('Rayess2016','',60,153812,109,309722),
 C('Kiss2018Claims','',268,481572,94,201013),
 C('Fintak2008','',3,12585,3,14320),
 C('Rohowetz2026','2005-2025',11,142691,8,48088),
 C('Forooghian2017','both compounded',22,54101,2,5614),
 C('Inoue2011','vs ranibizumab or pegaptanib',5,1209,0,4027),
 R('VanderBeek2015Bev','adjusted',0.66,0.39,1.09,'OR'),
 R('Bavinger2019','adjusted',1.25,0.97,1.61,'OR',invert=True,note='Reported as ranibizumab versus bevacizumab; inverted.'),
 C('Sigford2015','pooled literature',103,176124,79,269379,kind='meta',pool=False),
], 'Excludes the outbreak clusters from single compounding pharmacies.'))
IF.append(factor('abx','Topical antibiotic drops','Topical antibiotic (before and/or after injection)','No antibiotic',[
 C('Storey2014Abx','After injection',28,57654,11,34900),
 C('Bhavsar2016DRCR','DRCR.net trials',6,11565,3,17208),
 C('Cheung2012','5 days after',5,8259,2,5266),
 C('Torres-Costa2020','After injection, earlier period',5,14828,8,18687),
 C('Ramel2014','After injection',3,10144,3,1306),
 C('Park2013Korea','3 days before only',0,8649,3,8683),
 C('Comcali2025','After injection',7,13121,6,8218),
 C('Bhavsar2009DRCR','DRCR.net, several days after',3,2201,0,1276),
 C('Meyer2007','3 days after',2,860,0,984),
 M('Menchini2018','Meta-analysis, 13 studies',1.33,0.75,2.38,'OR'),
 M('BenoistdAzy2016','Meta-analysis, 9 studies',0.804,0.384,1.682,'OR'),
 M('Bande2017','Meta-analysis, 8 studies',1.7,1.08,2.66,'RR'),
 M('Liu2025AbxMA','Meta-analysis, 18 studies',1.85,0.72,4.76,'OR'),
 M('AleksanderIvanov2026','Meta-analysis, anti-VEGF subgroup',1.56,1.03,2.36,'OR'),
], 'Most comparisons are before-and-after changes of practice policy.'))
IF.append(factor('chx','Chlorhexidine instead of povidone-iodine','Chlorhexidine','Povidone-iodine 5 to 10%',[
 C('Bates2025','Aqueous chlorhexidine 0.05%',5,35284,59,267190),
 C('Ton2023','Aqueous chlorhexidine 0.05%',9,34893,13,32802),
 C('Mishra2023CHX','Chlorhexidine 0.05% in alcohol vs 10% iodine',20,31135,29,139817),
 M('CHX_MA_2025b','Meta-analysis, 3 comparative studies',1.27,0.50,3.22,'RR'),
 M('AqCHX_MA_2026','Meta-analysis, aqueous chlorhexidine, 8 studies',0.80,0.54,1.19,'OR'),
], 'Chlorhexidine is usually reserved for patients who report iodine sensitivity, so the groups differ.'))
IF.append(factor('nopi','No povidone-iodine','No antiseptic, or reduced povidone-iodine','Standard 5% povidone-iodine',[
 C('DilutePI2019','No povidone-iodine',1,319,12,29281),
 C('Modjtahedi2020IodineSens','No povidone-iodine',1,92,4,21185),
 C('Modjtahedi2020IodineSens','Reduced povidone-iodine (fewer drops)',6,769,4,21185,pool=False,note='Shares the reference group with the row above; not pooled.'),
 C('Bhavsar2016DRCR','No povidone-iodine',2,13,9,28773),
 C('DilutePI2019','Dilute povidone-iodine 0.625 to 2.5%',1,5460,12,29281,pool=False),
 C('Ajlan2025','2.5% applied twice',0,3826,3,3534,pool=False),
], 'Pooled estimate uses the no-antiseptic rows only. Very small exposed groups.'))
IF.append(factor('mask','Universal face masks','Period with universal masking (staff and patients)','Period before masking',[
 C('Patel2021Masks','12 US centres',45,211454,85,294514),
 C('Louis2023','',28,88218,73,219131),
 C('Liddle2024','',17,58270,29,63114),
 C('Haliyur2024','',10,38962,20,72194),
 C('Naguib2021','',7,22418,41,111679),
 C('Hebert2021','',7,14638,5,13863),
 C('Mirghorbani2022','Patients masked, mask not taped',7,19650,7,34277),
 C('Karimi2023','',4,10717,9,28085),
 C('Tanaka2022','',1,14725,1,31173),
 M('Tao2024Mask_NMA','Network meta-analysis',0.83,0.67,1.02,'OR'),
], 'Masking began with COVID-19 in 2020, when clinic volumes and patient mix also changed.'))
IF.append(factor('talk','Injector mask or no-talking policy','Injector face mask or no-talking policy','Usual care',[
 C('Starr2022Mayo','Injector mask vs no mask',41,98726,20,66098),
 C('Patel2020PhysMask','Injector mask vs no-talking (no mask)',9,30162,168,453460,pool=False,note='Different reference group (no-talking policy); not pooled.'),
 C('Garg2015','No-talking vs usual care, culture-positive',8,82658,9,47155,pool=False),
 M('Tao2024Mask_NMA','Injector mask, network meta-analysis',0.72,0.53,0.99,'OR'),
 M('Tao2024Mask_NMA','No-talking, network meta-analysis',0.56,0.39,0.82,'OR'),
], 'Oral streptococci cause a larger share of infections after injection than after surgery, which is the rationale for these measures.', pool=False))
IF.append(factor('setting','Clinic room instead of operating theatre','Office or clinic room','Operating theatre',[
 C('Abell2012','Single surgeon',4,3376,0,8873),
 C('Tabandeh2014','Two physicians',3,8647,2,3063),
 M('Li2021OfficeVsOR_MA','Meta-analysis, 31 studies',3.06,0.07,139.75,'OR'),
 M('Li2021OfficeVsOR_MA','Meta-analysis, culture-positive only',21.52,2.39,193.55,'OR'),
], 'Few direct comparisons. Large theatre-based European series report low rates (for example 10 in 134,701 [ref]), but they also differ in drug preparation and protocol.'))
IF.append(factor('bilat','Both eyes on the same day','Injection given in a same-day bilateral session','Injection given alone',[
 C('PakistanBilatOR2022','',1,15338,3,14920),
 C('BilatSpain2017','',0,1612,1,6560),
], 'A systematic review of 17 studies (69,239 bilateral sessions, 138,478 injections) found 33 cases and no case affecting both eyes [ref].'))
IF.append(factor('injector','Who gives the injection','Trainee or non-certified injector','Retina specialist',[
 C('Day2020Resident','Residents vs attendings',3,3877,2,3057),
 C('Falavarjani2015','Resident vs attending',1,4921,0,3116),
 R('Emerson2021Medicare','Not board-certified vs certified, adjusted',0.72,0.523,0.996,'OR',invert=True,note='Reported as certified versus not; inverted.'),
], 'Nurse-delivered programmes report rates in the usual range without a concurrent doctor comparison: see the rate table below the chart.'))
IF.append(factor('patient','Patient factors','Factor present','Factor absent',[
 R('Kulkarni2026','Diabetes',5.05,1.82,16.27,'OR',pool=False),
 C('Rayess2015Diagnosis','Diabetic eye disease vs neovascular AMD',12,24541,119,296017,pool=False),
 R('Patel2021Immuno','Systemic immunosuppression',9.86,4.0,24.3,'OR',pool=False),
 R('Baudin2018','Male vs female',0.78,0.63,0.96,'IRR',pool=False),
 R('Bavinger2019','Male vs female',1.3,1.05,1.61,'OR',invert=True,pool=False,note='Reported as female versus male; inverted.'),
], 'Different factors on one axis; not pooled.', pool=False))
IF.append(factor('anaes','Anaesthetic and antiseptic order','Technique named in the row','Comparison named in the row',[
 C('Uner2023','Lidocaine gel before iodine vs iodine first',59,65307,13,37601,pool=False),
 C('Lad2012','Lidocaine gel vs no gel',4,4682,4,4120,pool=False),
 R('Stem2019','2% lidocaine jelly, adjusted',11.28,3.39,37.46,'OR',pool=False),
 R('Stem2019','0.5% tetracaine gel, adjusted',3.95,1.15,13.50,'OR',pool=False),
 C('Tustin2014','Subconjunctival lidocaine vs other anaesthesia',0,6853,8,8189,pool=False),
 C('Levinson2017','Iodine reapplied after speculum vs not',2,11952,20,16155,pool=False),
], 'Gel applied before the iodine may shield the surface from it. Lid speculum use itself has no comparison with counts; one series reported 0.10% with a bladed speculum and 0.066% without (P = 0.27) [ref].', pool=False))
for f in IF:
    f['note']=f['note'].replace('[ref]','')
# refs used in notes
nurse=[]
for k,lab in [('Michelotti2014','Nurse injectors, UK'),('Precivia2021','Nurse injectors with guide device'),('Samalia2016NZ','Nurse injectors, New Zealand'),('NurseOzurdex2020','Nurse injectors, dexamethasone implant')]:
    o=ivi[k]['overall']; p,lo,hi=exact(o['cases'],o['injections']) if k!='Precivia2021' else exact(3,20421)
    n= o['injections'] if k!='Precivia2021' else 20421; x=o['cases'] if k!='Precivia2021' else 3
    nurse.append(dict(ref=cite(k),label=lab,cases=x,n=n,rate=p,lo=lo,hi=hi))
noteRefs={'setting':cite('Freiberg2017'),'bilat':cite('Okonkwo2024BilatSR'),'anaes':cite('Shah2011')}
bil_ub=1-0.05**(1/69239)

# pooled arm rates for PFS vs vial (for the cumulative calculator)
pf=[(17,112809,13,34790),(7,29899,41,68076),(12,78407,43,165347),(12,77925,21,93073),(96,427763,233,884061),(2,43684,10,108355),(6,31291,2,11689),(15,81460,9,87294),(1,3226,0,12366)]
pfs_rate=pool_rates([dict(cases=a,n=n) for a,n,_,_ in pf]); vial_rate=pool_rates([dict(cases=b,n=n) for _,_,b,n in pf])

# ===================== CATARACT: overall ===============================
CS=[ # key, IC status, label override, n, cases
 ('Pershing2020IRIS','Mixed or not recorded'),('Keay2012Medicare','No routine intracameral'),('Daien2016France','Mixed or not recorded'),
 ('WangJ2025Aier','No routine intracameral'),('Friling2022NCR_ISBCS','Routine intracameral'),('KimSH2019Korea','Mixed or not recorded'),
 ('Nowak2019Poland','Mixed or not recorded'),('Bro2026NCR','Routine intracameral'),('Mortemousque2026France','Mixed or not recorded'),
 ('Freeman2010Quebec','No routine intracameral'),('West2005Medicare','No routine intracameral'),('Song2025AravindAir','Routine intracameral'),
 ('Hatch2009Ontario','No routine intracameral'),('Llovet2024Baviera','Routine intracameral'),('Herrinton2016KP','Mixed or not recorded'),
 ('Shorstein2021KP','Routine intracameral'),('Yao2013China','Mixed or not recorded'),('Hage2024QuinzeVingts','Mixed or not recorded'),
 ('Wai2018Malaysia','Mixed or not recorded'),('Moussa2020WestMidlands','Mixed or not recorded'),('Ting2024Moorfields','Mixed or not recorded'),
 ('Day2015NOD','Mixed or not recorded'),('LiuN2025China','Mixed or not recorded'),('Li2004EPSWA','No routine intracameral'),
 ('Oshika2007Japan','No routine intracameral'),('Inoue2017Japan','Mixed or not recorded'),('MaX2020Cefuroxime','Mixed or not recorded'),
 ('Zhu2017China','Mixed or not recorded'),('Tan2012Singapore','Mixed or not recorded'),('Lalitha2005Aravind','No routine intracameral'),
 ('Ravindran2009Pondicherry','No routine intracameral'),('Malmin2020Norway','Mixed or not recorded'),('Rathi2020EPS','Routine intracameral'),
 ('WongTY2004SNEC','No routine intracameral'),('Khanna2015LVPEI','No routine intracameral'),('Rock2014Tubingen','Mixed or not recorded'),
 ('Ness2011Freiburg','No routine intracameral'),('Segretin2022BuenosAires','Mixed or not recorded'),('Haapala2005Finland','No routine intracameral'),
 ('Nentwich2015Munich','No routine intracameral'),('Schmitz1999GermanSurvey','Mixed or not recorded'),('Wykoff2010BPEI','Mixed or not recorded'),
 ('Norregaard1997Denmark','No routine intracameral'),('Rudnisky2014Alberta','Mixed or not recorded'),('Ianchulev2016Office','Mixed or not recorded'),
 ('Montan2002StErik','Routine intracameral'),('Bohigian1999PVI','No routine intracameral'),('Lloyd2009Toronto','No routine intracameral'),
 ('Moser2019Barcelona','Routine intracameral'),('TuniPicado2018','Routine intracameral'),('Weston2015InjectorIOL','Mixed or not recorded'),
 ('Miller2005BPEI','Mixed or not recorded'),('Jensen2008FQ','No routine intracameral'),('Moshirfar2006FQ','No routine intracameral'),
]
cser=[]
for k,ic in CS:
    r=cat[k]; o=r['overall']; n,x=o.get('operations_n'),o.get('cases')
    if not n or x is None or n<10000: continue
    cser.append(dict(key=k,n=n,cases=x,ic=ic,design=r['design'],country=r['country'],period=r['period']))
# Aravind split by intracameral use
cser.append(dict(key='Haripriya2019Aravind',n=993009,cases=692,ic='No routine intracameral',design='registry, before intracameral moxifloxacin',country='India',period='2011-2018',label='Haripriya 2019 no IC'))
cser.append(dict(key='Haripriya2019Aravind',n=1069634,cases=185,ic='Routine intracameral',design='registry, intracameral moxifloxacin',country='India',period='2011-2018',label='Haripriya 2019 IC'))
for s in cser:
    p,lo,hi=exact(s['cases'],s['n']); s.update(rate=p,lo=lo,hi=hi)
cser.sort(key=lambda s:-s['n'])
for s in cser:
    s['ref']=cite(s['key'])
    if 'label' not in s:
        a=REF[s['key']]['cite'].split(',')[0].split(' ')[0]; s['label']=a+' '+(cat[s['key']].get('year') or '')
cat_pool=pool_rates(cser)
cat_pool_grp={g:pool_rates([s for s in cser if s['ic']==g]) for g in ['Routine intracameral','No routine intracameral','Mixed or not recorded']}
since2010=[s for s in cser if re.search(r'(20(1\d|2\d))',s['period'] or '') and int(re.findall(r'(\d{4})',s['period'])[0])>=2010]
cat_pool_recent=pool_rates(since2010)
cat_metas=[]
for k,lab,rate in [('Taban2005','215 studies, 1963-2003',None),('Shi2022','39 studies, phacoemulsification',None),('Kato2022','51 articles, network meta-analysis',None),('Cao2013','42 studies',None)]:
    r=cat[k]; o=r['overall']; cat_metas.append(dict(ref=cite(k),label=lab,rate=float(str(o['rate_pct']).split()[0])/100,lo=None,hi=None,cases=o.get('cases'),n=o.get('operations_n')))
cat_metas[1]['lo'],cat_metas[1]['hi']=0.00083,0.00101
cite('Zafar2021Medicare')

# ===================== CATARACT FACTORS ===============================
CF=[]
CF.append(factor('ic','Intracameral antibiotic','Intracameral antibiotic at the end of surgery','No intracameral antibiotic',[
 C('Barry2006ESCRSprelim','Cefuroxime 1 mg, ESCRS randomised trial',5,6836,23,6862,group='Randomised trials'),
 C('Sharma2026RCT','Moxifloxacin 0.5 mg, randomised trial',6,30000,16,30000,group='Randomised trials'),
 C('Melega2019RCT','Moxifloxacin 0.15 mg, randomised trial',1,1818,7,1822,group='Randomised trials'),
 C('Haripriya2019Aravind','Moxifloxacin, Aravind',185,1069634,692,993009,group='Moxifloxacin'),
 C('Shenoy2021Rural','Moxifloxacin',92,112967,179,101815,group='Moxifloxacin'),
 C('Bhatta2021Nepal','Moxifloxacin',8,31340,116,80643,group='Moxifloxacin'),
 C('DaveVP2022EPS3','Moxifloxacin',15,34318,21,32649,group='Moxifloxacin'),
 C('Matsuura2013Japan','Moxifloxacin',3,18794,8,15958,group='Moxifloxacin'),
 C('Vieira2017Brazil','Moxifloxacin',1,3680,8,3515,group='Moxifloxacin'),
 C('Galvis2014Colombia','Moxifloxacin',0,1618,1,1056,group='Moxifloxacin'),
 C('Friling2018StErik','Cefuroxime',33,75144,62,34390,group='Cefuroxime'),
 C('Katz2015Israel','Cefuroxime',10,29431,22,26663,group='Cefuroxime'),
 C('Ng2016HongKong','Cefuroxime',0,23096,8,7332,group='Cefuroxime'),
 C('RodriguezCaravaca2013Alcorcon','Cefuroxime',5,12868,39,6595,group='Cefuroxime'),
 C('deSousaCasavechia2023','Cefuroxime',3,16977,6,6207,group='Cefuroxime'),
 C('ConciBrazil2023','Cefuroxime',4,13495,7,3407,group='Cefuroxime'),
 C('Beselga2013Portugal','Cefuroxime',0,13390,6,2299,group='Cefuroxime'),
 C('SharmaS2015Bhubaneswar','Cefuroxime',8,7366,12,7756,group='Cefuroxime'),
 C('Barreau2012Limoges','Cefuroxime',1,2289,35,2826,group='Cefuroxime'),
 C('Tan2012Singapore','Cefazolin',2,20638,19,29539,group='Cefazolin'),
 C('RomeroAroca2012Reus','Cefazolin',7,13305,76,11696,group='Cefazolin'),
 C('Garat2009Barcelona','Cefazolin',6,12649,25,5930,group='Cefazolin'),
 C('Au2016Westmead','Vancomycin',6,12266,11,2539,group='Vancomycin'),
 C('Rush2015ASC','Vancomycin',0,9386,11,11333,group='Vancomycin'),
 C('Li2019ColeEye','Surgeons preferring intracameral vs topical',11,16325,29,16201,group='Mixed agents'),
 R('Herrinton2016KP','Kaiser Permanente, adjusted',0.58,0.38,0.91,'OR',group='Adjusted cohorts'),
 R('Daien2016France','France national, cefuroxime, adjusted',0.61,0.56,0.68,'OR',group='Adjusted cohorts'),
 R('CreuzotGarcher2016France','France hospital database, adjusted',0.53,0.50,0.57,'OR',group='Adjusted cohorts'),
 M('Kessel2015','Meta-analysis, 1 trial and 17 cohorts',0.12,0.08,0.18,'RR'),
 M('Huang2016','Meta-analysis, vancomycin or moxifloxacin',0.20,0.10,0.42,'RR'),
 M('Bowen2018','Meta-analysis, 17 studies',0.20,0.13,0.32,'OR'),
 M('Kato2022','Network meta-analysis, any intracameral',0.19,0.12,0.30,'OR'),
 M('Ramos2026','Bayesian meta-analysis, 54 studies',0.208,0.161,0.268,'OR'),
 M('AbuZaid2026','Meta-analysis of randomised trials, moxifloxacin',0.22,0.07,0.77,'RR'),
], 'Apart from the three randomised trials, these are before-and-after comparisons at the same unit, in which other changes over time also occurred. Several pre-intracameral baselines (Spain, France) were 0.4 to 1.2%, well above contemporary rates, which inflates the ratio.',
 subpools=['Randomised trials','Moxifloxacin','Cefuroxime','Cefazolin','Vancomycin']))
CF.append(factor('icagent','Moxifloxacin instead of cefuroxime','Intracameral moxifloxacin','Intracameral cefuroxime',[
 C('Shorstein2021KP','Kaiser Permanente',21,106854,13,97801),
 C('Friling2024NCR','Swedish register',25,122340,126,514916),
 C('Rathi2020EPS','15 rural Indian centres',12,24534,3,17932,derived=True,note='Counts derived from reported rates (0.049% and 0.017%).'),
], 'Kaiser used 0.1 mL of 0.1% moxifloxacin; the authors suggest the dose may have been low for deeper anterior chambers.'))
CF.append(factor('pcr','Posterior capsule rupture','Posterior capsule rupture or vitreous loss','Uncomplicated surgery',[
 C('Haripriya2019Aravind','Aravind, no intracameral antibiotic',63,14505,629,978504),
 C('Haripriya2019Aravind','Aravind, with intracameral moxifloxacin',25,13847,160,1055787,pool=False,note='Same registry as the row above; not pooled.'),
 R('Day2015NOD','UK National Ophthalmology Database',7.94,3.35,18.83,'OR'),
 R('Herrinton2016KP','Kaiser Permanente, adjusted',3.68,1.89,7.20,'OR'),
 R('CreuzotGarcher2016France','France, adjusted',5.24,4.11,6.68,'OR'),
 R('Segers2022EUREQUO_PCR','European registry, adjusted',4.4,2.48,7.81,'OR'),
 R('Hatch2009Ontario','Ontario, anterior vitrectomy, adjusted',9.56,6.43,14.2,'OR'),
 R('KimSH2019Korea','Korea, anterior vitrectomy',8.71,6.71,11.32,'OR'),
 R('Yao2013China','China, vitreous communication',9.96,4.54,21.84,'OR'),
 R('Wai2018Malaysia','Malaysia registry',3.773,1.915,7.432,'OR'),
 R('Lalitha2005Aravind','Aravind 2002-2003',6.57,1.5,28.5,'RR'),
 R('WongTY2004SNEC','Singapore',8.0,3.1,20.7,'RR'),
 R('Khanna2015LVPEI','LV Prasad',6.98,2.22,21.98,'OR'),
 R('ESCRS2007','ESCRS trial, any surgical complication',4.95,1.68,14.6,'OR'),
 M('Cao2013','Meta-analysis',6.34,4.23,9.52,'OR'),
 M('Zhang2025','Meta-analysis',9.18,3.31,25.43,'OR'),
 M('ChenKY2025','Meta-analysis',3.434,1.789,6.591,'OR'),
], 'Communication with the vitreous lets organisms reach a medium that supports their growth.'))
CF.append(factor('pcric','Intracameral antibiotic when the capsule ruptures','Intracameral antibiotic','None, in eyes with posterior capsule rupture',[
 C('Haripriya2019Aravind','Moxifloxacin, Aravind',25,13847,63,14505),
 R('Herrinton2016KP','Kaiser Permanente, adjusted',0.49,0.08,3.09,'OR'),
 M('Arteche2026','Meta-analysis, 9 studies',0.41,0.26,0.65,'OR'),
], 'France reported 0.37% with cefuroxime and 0.51% without in eyes with capsular rupture (P = .001); group sizes were not given in the abstract.'))
CF.append(factor('incision','Clear corneal incision','Clear corneal incision','Scleral tunnel or limbal incision',[
 C('Nagaki2003RCT','Randomised trial, Japan',11,3831,2,3901),
 C('Oshika2007Japan','Japan national survey',9,20894,38,78170),
 C('Miller2005BPEI','Bascom Palmer',6,11462,1,4458),
 R('ESCRS2007','ESCRS trial cohort',5.88,1.34,25.9,'OR'),
 R('Schmitz1999GermanSurvey','German national survey',0.35,0.24,0.51,'OR',invert=True,note='Reported as scleral versus clear corneal; inverted.'),
 M('Taban2005','Meta-analysis 2005',2.55,1.75,3.71,'RR'),
 M('Cao2013','Meta-analysis, prospective studies',5.65,3.75,8.52,'OR'),
 M('Cao2013','Meta-analysis, retrospective studies',2.44,1.46,4.09,'OR'),
], 'Most data predate routine intracameral antibiotics and modern square, multiplanar wound construction. The Swedish register found 0.053% versus 0.036% (P = .14).'))
CF.append(factor('topical','Topical antibiotic drops','Topical antibiotic (named in the row)','Comparison named in the row',[
 R('Gower2017Cochrane','Levofloxacin drops vs placebo (ESCRS)',0.72,0.32,1.61,'RR',pool=False),
 R('Herrinton2016KP','Drops added to intracameral vs intracameral alone',1.63,0.48,5.47,'OR',pool=False),
 C('Jensen2008FQ','Fourth- vs third-generation fluoroquinolone',7,12566,33,16710,pool=False),
 R('Herrinton2016KP','No recorded antibiotic vs topical only',1.95,1.22,3.11,'OR',pool=False),
 M('Huang2016','Meta-analysis, topical vs none',0.65,0.43,0.99,'RR'),
], 'Swedish and French national data report no additional benefit of drops on top of an intracameral antibiotic.', pool=False))
CF.append(factor('isbcs','Both eyes on the same day','Immediate sequential bilateral surgery','Delayed or unilateral surgery',[
 C('Friling2022NCR_ISBCS','Swedish register',14,92238,408,1364934),
 C('Herrinton2017KP_ISBCS','Kaiser Permanente',1,10494,2,38736),
 R('Lacy2021IRIS_ISBCS','IRIS Registry, adjusted',1.08,0.87,1.31,'OR'),
], 'Same-day bilateral surgery is offered to selected patients, so lower-risk eyes are over-represented. Medicare reported 1.74 versus 1.01 per 1,000 (P = .15).'))
CF.append(factor('male','Male sex','Male','Female',[
 R('Keay2012Medicare','Medicare',1.23,1.15,1.31,'RR'),
 R('CreuzotGarcher2016France','France',1.48,1.40,1.56,'OR'),
 R('KimSH2019Korea','Korea',1.42,1.25,1.61,'OR'),
 R('Hatch2009Ontario','Ontario',1.4,1.19,1.64,'OR'),
 R('Freeman2010Quebec','Quebec',1.44,1.24,1.66,'OR'),
 R('Zafar2021Medicare','Medicare 2011-2019',0.89,0.86,0.92,'OR',invert=True,note='Reported as female versus male; inverted.'),
 R('Tan2012Singapore','Singapore',2.96,1.15,7.65,'OR'),
 R('Wai2018Malaysia','Malaysia',2.121,1.464,3.015,'OR'),
 R('ESCRS2007','ESCRS, proven cases',2.7,1.07,6.8,'OR'),
 M('Cao2013','Meta-analysis',1.41,1.22,1.63,'OR'),
 M('Zhang2025','Meta-analysis',0.81,0.75,0.87,'OR',invert=True,note='Reported as female versus male; inverted.'),
], 'Keay 2012 and Zafar 2021 are separate Medicare periods.'))
CF.append(factor('diabetes','Diabetes','Diabetes','No diabetes',[
 R('KimSH2019Korea','Korea',1.59,1.31,1.93,'OR'),
 R('Hou2023TaiwanDM','Taiwan, adjusted',1.09,1.03,1.16,'OR'),
 R('Herrinton2016KP','Diabetic retinopathy, Kaiser, adjusted',1.31,0.87,1.98,'OR'),
 M('ChenKY2025','Meta-analysis',1.174,1.109,1.242,'OR'),
 M('Zhang2025','Meta-analysis (I² 95%)',4.9,2.41,9.95,'OR'),
], 'The two meta-analyses disagree by a factor of four; the larger figure comes from highly heterogeneous small studies.'))
CF.append(factor('other','Other patient and surgical factors','Factor named in the row','Comparison named in the row',[
 R('Keay2012Medicare','Age 85 or over vs 65-74',1.53,1.38,1.69,'RR',pool=False,group='Patient'),
 R('Herrinton2016KP','Age 80 or over vs under 70',1.46,1.05,2.05,'OR',pool=False,group='Patient'),
 R('Ghoraba2025IRIS','Uveitis',7.0,5.88,8.33,'OR',pool=False,group='Patient'),
 R('Ghoraba2025IRIS','Systemic immunosuppressive drugs',2.37,1.91,2.95,'OR',pool=False,group='Patient'),
 R('Keay2012Medicare','Black vs White',1.17,1.03,1.33,'RR',pool=False,group='Patient'),
 R('Zafar2021Medicare','Cataract combined with retinal surgery',2.6,2.15,3.16,'OR',pool=False,group='Surgery'),
 R('CreuzotGarcher2016France','Any combined procedure',1.77,1.53,2.05,'OR',pool=False,group='Surgery'),
 R('Zafar2021Medicare','Previous glaucoma surgery',1.4,1.18,1.65,'OR',pool=False,group='Surgery'),
 M('Cao2013','Extracapsular or intracapsular vs phaco, meta-analysis',2.19,1.40,3.42,'OR'),
 R('Lalitha2005Aravind','Extracapsular vs phaco, culture-positive',4.9,1.2,19.3,'RR',pool=False,group='Surgery'),
 R('ESCRS2007','Silicone vs acrylic lens',3.13,1.47,6.67,'OR',pool=False,group='Surgery'),
 R('KimSH2019Korea','Sutured vs sutureless wound',2.82,2.13,3.74,'OR',pool=False,group='Surgery'),
 C('Thoms2007Suture','Unsutured vs sutured clear corneal wound',5,379,0,436,pool=False,group='Surgery'),
 R('Wai2018Malaysia','Longer surgery',2.114,1.473,3.032,'OR',pool=False,group='Surgery'),
 R('GarciaArumi2007','Surgery over 45 minutes',7.2,1.7,29.7,'OR',pool=False,group='Surgery'),
 R('Segers2022EUREQUO_Anaesth','Topical vs regional anaesthesia',0.6,0.44,0.82,'OR',invert=True,pool=False,group='Surgery',note='Reported as regional versus topical; inverted.'),
 R('Keay2012Medicare','Surgeon 1-50 cases a year vs over 1,000',3.8,3.13,4.61,'RR',pool=False,group='Surgeon and setting'),
 R('Campbell2019Ontario','Late-career vs mid-career surgeon',1.41,1.01,1.98,'OR',pool=False,group='Surgeon and setting'),
 C('Kunkler2022Resident','Resident vs attending surgeon',6,6447,16,26058,pool=False,group='Surgeon and setting'),
 C('Song2025AravindAir','Laminar flow theatre vs split air conditioning',15,116981,22,114626,pool=False,group='Surgeon and setting'),
 C('Speaker1991PVI','5% povidone-iodine vs silver protein',2,3489,11,4594,pool=False,group='Antisepsis'),
 C('Nentwich2015Munich','1% iodine irrigation vs 1% iodine drop',6,14676,15,6696,pool=False,group='Antisepsis'),
], 'Single estimates; not pooled. Each row states its own comparison.', pool=False))

out=dict(
 ivi=dict(series=ser,pooled=ivi_pool,pooled_groups=ivi_pool_grp,metas=ivi_metas,cumobs=cumobs,factors=IF,nurse=nurse,
          pfs_rate=pfs_rate,vial_rate=vial_rate,cumRef=cumRef,bilat_ub=bil_ub,excluded=excluded,noteRefs=noteRefs),
 cat=dict(series=cser,pooled=cat_pool,pooled_groups=cat_pool_grp,pooled_recent=cat_pool_recent,metas=cat_metas,factors=CF),
)
# ic arm rates from comparative studies
ic_rows=[r for r in CF[0]['rows'] if r['a'] and r['kind']=='study']
def arm(s): x,n=s.split('/'); return dict(cases=int(x.replace(',','')),n=int(n.replace(',','')))
out['cat']['ic_arm']=pool_rates([arm(r['a']) for r in ic_rows]); out['cat']['noic_arm']=pool_rates([arm(r['b']) for r in ic_rows])
for k,_ in excluded: cite(k)
refs=[dict(n=i+1,**REF[k]) for i,k in enumerate(used)]
num={k:i+1 for i,k in enumerate(used)}
def renum(o):
    if isinstance(o,dict):
        for kk,v in o.items():
            if kk=='ref' and isinstance(v,str): o[kk]=num[v]
            else: renum(v)
    elif isinstance(o,list):
        for v in o: renum(v)
renum(out)
out['ivi']['noteRefs']={k:num[v] if isinstance(v,str) else v for k,v in noteRefs.items()}
out['ivi']['cumRef']=num[cumRef]
out['refs']=refs
json.dump(out,open(D+'data.json','w'),indent=0,default=float)
print('refs',len(refs),'ivi series',len(ser),'cat series',len(cser))
print('IVI pooled',ivi_pool); print({g:(round(v['est']*100,4),round(v['lo']*100,4),round(v['hi']*100,4),v['k']) for g,v in ivi_pool_grp.items()})
print('PFS',pfs_rate,'vial',vial_rate)
print('CAT pooled',cat_pool); print({g:(round(v['est']*100,4),round(v['lo']*100,4),round(v['hi']*100,4),v['k']) for g,v in cat_pool_grp.items()})
print('cat recent',cat_pool_recent); print('ic arm',out['cat']['ic_arm'],'noic',out['cat']['noic_arm'])
for f in IF+CF:
    p=f.get('pooled'); print(f['id'], p and (round(p['est'],2),round(p['lo'],2),round(p['hi'],2),p['k'],p['i2']), [ (s['group'],round(s['est'],2),round(s['lo'],2),round(s['hi'],2),s['k']) for s in f.get('subpools',[])])
print('excluded',excluded)
