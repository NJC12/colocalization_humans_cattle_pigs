import csv, glob, json, os, re, collections, bisect
S="/n/scratch/users/n/njc12/snakemake/farm_sims_for_publication"
O="/n/scratch/users/n/njc12/pubfreeze_sandbox"
ARMS=["causal_maf001_paired","causal_maf001_unpaired","causal_power_n200000","causal_power_n30000"]

# Selection-coefficient bins, carried over verbatim from the round-3 dashboard so
# the two are readable side by side.
SELCO_BREAKS=[-float("inf"),-9.8e-4,-4.9e-4,-2.4e-4,-1.2e-4,-6.1e-5,-3.1e-5,0.0,float("inf")]
SELCO_LABELS=["< -9.8e-4","-9.8e-4 .. -4.9e-4","-4.9e-4 .. -2.4e-4","-2.4e-4 .. -1.2e-4",
              "-1.2e-4 .. -6.1e-5","-6.1e-5 .. -3.1e-5","-3.1e-5 .. 0","neutral (0)"]
def selco_bin(v):
    i=bisect.bisect_right(SELCO_BREAKS,v)-1
    return min(max(i,0),len(SELCO_LABELS)-1)
MAF_LABELS=[f"{i*5}–{(i+1)*5}%" for i in range(10)]
def maf_bin(m):
    if m is None or m<=0: return None
    return min(int(m/0.05),9)

# MAF histogram at 0.001 resolution -> quantiles without holding every value.
NB=501
# Keyed by (category, selection bin) and deduplicated by (category, replicate).
# The GWAS variant panel is a property of the stage-1 tree sequence and the
# min_maf floor, both of which are identical across the four sampling strategies
# -- verified: the three 15-replicate arms report the same variant count to the
# unit. Keying it per arm would count each panel up to four times and let the
# arm with fewer finished replicates silently reweight the distribution.
hist=collections.defaultdict(lambda:[0]*NB)          # (cat,selcobin) -> hist
seen_panel=set()                                     # (cat,rundir)
def add_maf(k,m):
    hist[k][min(int(round(m*1000)),NB-1)]+=1
def quants(h):
    n=sum(h)
    if n==0: return None
    cum=0; want={0.0:None,0.25:None,0.5:None,0.75:None,1.0:None}; qs=sorted(want)
    out={}; qi=0
    for i,c in enumerate(h):
        if not c: continue
        lo=cum; cum+=c
        while qi<len(qs):
            target=qs[qi]*(n-1)
            if lo<=target<cum: out[qs[qi]]=i/1000.0; qi+=1
            else: break
    for q in qs: out.setdefault(q,None)
    return dict(n=n,min=out[0.0],q1=out[0.25],med=out[0.5],q3=out[0.75],max=out[1.0])

# traits, keyed by workdir
traits=collections.defaultdict(list)
with open(O+"/dash_traits.tsv") as fh:
    for r in csv.DictReader(fh,delimiter="\t"):
        traits[(r["root"], r["rep"].rstrip("/"))].append(r)

runs=[]
with open(S+"/RUNS.tsv") as fh:
    rd=csv.DictReader(fh,delimiter="\t")
    for r in rd: runs.append(r)

outcome=collections.Counter()   # (arm,cat,panel,rcp,mafbin,outcome) -> n
done=0
for r in runs:
    arm, rundir = r["arm"], r["run_dir"]
    wd=os.path.join(S,arm,rundir)
    vf=glob.glob(os.path.join(wd,"stage2","*","gwas_*_gtex_*_maf_*","*gwas_vars_*.tsv"))
    if not vf: continue
    pos2maf={}
    panel_key=(r["category"], rundir)
    novel = panel_key not in seen_panel
    seen_panel.add(panel_key)
    with open(vf[0]) as fh:
        for v in csv.DictReader(fh,delimiter="\t"):
            try:
                m=float(v["maf"]); sc=float(v["selco"]); p=int(float(v["position"]))
            except (ValueError,KeyError,TypeError): continue
            if m>0.5: m=1.0-m
            if novel: add_maf((r["category"],selco_bin(sc)),m)
            pos2maf[p]=m
    for t in traits.get((arm, rundir),()):
        mb=maf_bin(pos2maf.get(int(t["trait"][2:])))
        if mb is None: continue
        for rcp,scol,ocol in ((50,"self_rcp","n_other_50"),(90,"self_rcp","n_other_90")):
            thr=0.5 if rcp==50 else 0.9
            s=float(t[scol]); o=int(t[ocol])
            oc = "correct" if s>thr else ("wrong" if o>0 else "none")
            outcome[(arm,t["cat"],int(t["gtex_n"]),rcp,mb,oc)]+=1
    done+=1
    if done%50==0: print(f"  {done}/{len(runs)} runs", flush=True)

box=[]
for (cat,sb),h in hist.items():
    q=quants(h)
    if q: box.append(dict(cat=cat,sb=sb,**q))
stack=[dict(arm=a,cat=c,panel=p,rcp=rc,mb=mb,oc=oc,n=n) for (a,c,p,rc,mb,oc),n in outcome.items()]
json.dump(dict(selco_labels=SELCO_LABELS,maf_labels=MAF_LABELS,box=box,stack=stack),
          open(O+"/dash_maf.json","w"),separators=(",",":"))
assert stack, "outcome join produced nothing -- the (root,rep) key did not match any run"
print("runs scanned:",done," box rows:",len(box)," stack rows:",len(stack))
print("bytes:",os.path.getsize(O+"/dash_maf.json"))
