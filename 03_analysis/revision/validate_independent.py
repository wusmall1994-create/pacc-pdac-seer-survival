"""Validate all new OS estimates directly from raw export using Python stdlib.
No R objects or survival libraries are used. Outputs only aggregate checks.
"""
import csv, math, sys
from collections import Counter
from pathlib import Path

project = Path(sys.argv[1]) if len(sys.argv)>1 else Path('.')
out = Path(sys.argv[2]) if len(sys.argv)>2 else project/'04_results/revision/conditional_OS'
data=[]
with (project/'01_data_extraction/pacc_pdac_seer17_2004_2023_raw.txt').open(encoding='utf-8-sig',newline='') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        code=r['Histologic Type ICD-O-3']
        if code not in ('8140','8500','8550'): continue
        if r['Diagnostic Confirmation'] not in ('Positive histology','Pos hist AND immunophenotyping AND/OR pos genetic studies'): continue
        if r['Type of Reporting Source'] in ('Autopsy only','Death certificate only'): continue
        if r['Survival Days']=='NA' or r['Vital status recode (study cutoff used)']=='NA': continue
        st=r['Combined Summary Stage with Expanded Regional Codes (2004+)']
        st=next((s for s in ('Localized','Regional','Distant') if st.startswith(s)),'Unknown')
        data.append((code,int(r['Year of diagnosis']),st,max(float(r['Survival Days']),.5)/365.25,int(r['Vital status recode (study cutoff used)']=='Dead')))
assert len(data)==116516
rows=list(csv.DictReader((out/'conditional_OS_all_analyses.csv').open()))
checks=[]; cache={}
for q in rows:
    key=(q['comparator'],int(q['diagnosis_end']),q['histology'],q['stage'])
    if key not in cache:
        comp,cutoff,h,st=key
        cache[key]=[(t,e) for code,yr,s,t,e in data if yr<=cutoff and (st=='All' or s==st) and (('pACC' if code=='8550' else 'PDAC')==h) and (comp=='8140+8500' or code in ('8500','8550'))]
    source=cache[key]; lm=float(q['landmark']); hz=float(q['horizon'])
    z=[(t-lm,e) for t,e in source if t>lm]
    supported=len(z)>=2 and max((t for t,e in z),default=-1)>=hz
    assert supported==(q['supported']=='TRUE')
    assert len(source)==int(q['n_diagnosis']) and len(z)==int(q['n_landmark'])
    assert sum(t>=hz for t,e in z)==int(q['n_at_horizon'])
    assert sum(e==1 and t<=hz for t,e in z)==int(q['deaths_within_horizon'])
    assert sum(e==0 and t<hz for t,e in z)==int(q['censored_before_horizon'])
    errors=[0.,0.,0.]
    if supported:
        all_times=Counter(t for t,e in z); deaths=Counter(t for t,e in z if e)
        risk=len(z); surv=1.; greenwood=0.
        for t in sorted(all_times):
            if t>hz: break
            ev=deaths[t]
            if ev:
                surv*=1-ev/risk
                if risk>ev: greenwood+=ev/(risk*(risk-ev))
            risk-=all_times[t]
        if 0<surv<1:
            se=math.sqrt(greenwood)/abs(math.log(surv))
            center=math.log(-math.log(surv)); width=1.959963984540054*se
            lower=math.exp(-math.exp(center+width)); upper=math.exp(-math.exp(center-width))
        else: lower=upper=surv
        for j,(v,col) in enumerate(zip((surv,lower,upper),('estimate','lower95','upper95'))):
            if q[col]!='NA':
                errors[j]=abs(v-float(q[col])); assert errors[j]<1e-10,(key,lm,hz,col,v,q[col])
    checks.append(dict(comparator=q['comparator'],cohort=q['cohort'],diagnosis_end=q['diagnosis_end'],histology=q['histology'],stage=q['stage'],landmark=q['landmark'],horizon=q['horizon'],supported=supported,estimate_error=errors[0],lower_error=errors[1],upper_error=errors[2],counts_match=True))
with (out/'independent_python_validation.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(checks[0]));w.writeheader();w.writerows(checks)
print('Validated rows:',len(checks))
print('Supported rows:',sum(x['supported'] for x in checks))
print('Maximum estimate/CI absolute errors:',*[max(x[c] for x in checks) for c in ('estimate_error','lower_error','upper_error')])
print('All raw-cohort counts, risk sets, window deaths and censor counts agree.')
