"""Independent raw-export checks of selection tables and era KM estimates."""
import csv,math,sys,statistics
from pathlib import Path
from collections import Counter
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path('.')
out=Path(sys.argv[2]) if len(sys.argv)>2 else root/'04_results/revision/selection_and_calendar'
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig')))
data=[];raw=[]
with (root/'01_data_extraction/pacc_pdac_seer17_2004_2023_raw.txt').open(encoding='utf-8-sig',newline='') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  code=x['Histologic Type ICD-O-3']
  if code not in ('8140','8500','8550'):continue
  year=int(x['Year of diagnosis']);sc=x['RX Summ--Surg Prim Site (1998-2022)'] if year<=2022 else x['RX Summ--Surg Prim Site 2023 (2023+)']
  stage=next((s for s in ('Localized','Regional','Distant') if x['Combined Summary Stage with Expanded Regional Codes (2004+)'].startswith(s)),'Unknown')
  y=dict(code=code,histology='pACC' if code=='8550' else 'PDAC',year=year,age=int(x['Age recode with single ages and 85+'].split(' ')[0].replace('+','')),sex=x['Sex'],stage=stage,era='2004-2009' if year<=2009 else '2010-2015' if year<=2015 else '2016-2017' if year<=2017 else '2018-2023',site='Head' if x['Primary Site - labeled'].startswith('C25.0') else 'Body/tail' if x['Primary Site - labeled'].startswith(('C25.1','C25.2')) else 'Other/unspecified',surgery='No/unknown' if sc in ('NA','Blank(s)','00','A000','99','A990') else 'Yes',chemotherapy='Yes' if x['Chemotherapy recode (yes, no/unk)']=='Yes' else 'No/unknown',confirmation=x['Diagnostic Confirmation'])
  y['confirmed']=y['confirmation'] in ('Positive histology','Pos hist AND immunophenotyping AND/OR pos genetic studies')
  y['eligible']=x['Type of Reporting Source'] not in ('Autopsy only','Death certificate only') and x['Survival Days']!='NA' and x['Vital status recode (study cutoff used)']!='NA'
  raw.append(y)
  if y['confirmed'] and y['eligible']:data.append((code,year,max(float(x['Survival Days']),.5)/365.25,int(x['Vital status recode (study cutoff used)']=='Dead')))
for q in read(out/'selection_flow.csv'):
 z=[x for x in raw if x['histology']==q['histology']]
 vals=[len(z),sum(x['confirmed'] for x in z),sum(not x['confirmed'] for x in z),sum(x['confirmed'] and not x['eligible'] for x in z),sum(x['confirmed'] and x['eligible'] for x in z),sum(not x['confirmed'] and x['eligible'] for x in z)]
 assert vals==[int(q[k]) for k in ('target','confirmation_pass','confirmation_excluded','later_excluded','final','excluded_confirmation_only_otherwise_eligible')]
checks=[];cache={}
for q in read(out/'included_excluded_characteristics.csv'):
 key=(q['universe'],q['histology'])
 if key not in cache:
  z=[x for x in raw if x['histology']==key[1] and (key[0]=='confirmation_step' or x['eligible'])];cache[key]=([x for x in z if x['confirmed']],[x for x in z if not x['confirmed']])
 a,b=cache[key];v=q['variable'];assert len(a)==int(q['included_n']) and len(b)==int(q['excluded_n'])
 if v in ('age','year'):
  av=[x[v] for x in a];bv=[x[v] for x in b];smd=abs(statistics.mean(av)-statistics.mean(bv))/math.sqrt((statistics.variance(av)+statistics.variance(bv))/2)
  assert statistics.median(av)==float(q['included'].split()[0]) and statistics.median(bv)==float(q['excluded'].split()[0])
 else:
  na=sum(x[v]==q['level'] for x in a);nb=sum(x[v]==q['level'] for x in b);assert na==int(q['included'].split()[0]) and nb==int(q['excluded'].split()[0]);pa=na/len(a);pb=nb/len(b);den=math.sqrt((pa*(1-pa)+pb*(1-pb))/2);smd=abs(pa-pb)/den if den else (0 if pa==pb else math.nan)
 if q['SMD']!='NA':assert abs(smd-float(q['SMD']))<1e-10
 checks.append(dict(check='selection characteristic',key=str(key)+(v,q['level']).__str__(),passed=True))
maxerr=0
for q in read(out/'era_conditional_OS.csv'):
 lm=float(q['landmark']);hz=float(q['horizon']);lo=int(q['diagnosis_start']);hi=int(q['diagnosis_end']);g=q['group']
 z0=[(t,e) for c,y,t,e in data if lo<=y<=hi and (c=='8550' if g=='pACC' else c in ('8140','8500') if g=='PDAC pooled' else c=='8500')];z=[(t-lm,e) for t,e in z0 if t>lm]
 assert [len(z0),len(z),sum(t>=hz for t,e in z),sum(e and t<=hz for t,e in z),sum(not e and t<hz for t,e in z)]==[int(q[k]) for k in ('n_diagnosis','n_landmark','n_end','deaths','censored')]
 supported=len(z)>=2 and max((t for t,e in z),default=-1)>=hz;assert supported==(q['supported']=='TRUE')
 if supported:
  times=Counter(t for t,e in z);deaths=Counter(t for t,e in z if e);risk=len(z);s=1.;gw=0.
  for t in sorted(times):
   if t>hz:break
   ev=deaths[t]
   if ev:
    s*=1-ev/risk
    if risk>ev:gw+=ev/(risk*(risk-ev))
   risk-=times[t]
  low=up=s
  if 0<s<1:
   cen=math.log(-math.log(s));wid=1.959963984540054*math.sqrt(gw)/abs(math.log(s));low=math.exp(-math.exp(cen+wid));up=math.exp(-math.exp(cen-wid))
  for val,k in zip((s,low,up),('estimate','lower','upper')):
   if q[k]!='NA':err=abs(val-float(q[k]));maxerr=max(maxerr,err);assert err<1e-10
 checks.append(dict(check='era KM and risk sets',key=str((g,q['era'],lm,hz)),passed=True))
with (out/'independent_remaining_validation.csv').open('w',newline='',encoding='utf8') as f:
 w=csv.DictWriter(f,fieldnames=checks[0]);w.writeheader();w.writerows(checks)
print('All',len(checks),'selection/era checks passed. Maximum KM/CI error:',maxerr)
