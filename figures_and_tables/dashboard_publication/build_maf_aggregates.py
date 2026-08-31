import csv, glob, json, os, collections, bisect
S="/n/scratch/users/n/njc12/snakemake/farm_sims_for_publication"
O="/n/scratch/users/n/njc12/pubfreeze_sandbox"

SELCO_BREAKS=[-float("inf"),-9.8e-4,-4.9e-4,-2.4e-4,-1.2e-4,-6.1e-5,-3.1e-5,0.0,float("inf")]
SELCO_LABELS=["< -9.8e-4","-9.8e-4 – -4.9e-4","-4.9e-4 – -2.4e-4","-2.4e-4 – -1.2e-4",
              "-1.2e-4 – -6.1e-5","-6.1e-5 – -3.1e-5","-3.1e-5 – 0","0 (neutral)"]
def selco_bin(v):
    return min(max(bisect.bisect_right(SELCO_BREAKS,v)-1,0),len(SELCO_LABELS)-1)
MAF_LABELS=[f"({i*5/100:.2f},{(i+1)*5/100:.2f}]" for i in range(10)]
def maf_bin(m):
    if m is None or m<=0: return None
    return min(int(m/0.05),9)

NB=501                                    # MAF histogram at 0.001 resolution
def stats_from_hist(h):
    """Tukey box stats without holding every value: quartiles from the histogram,
    whiskers to the most extreme datum inside 1.5 IQR, outliers counted."""
    n=sum(h)
    if n==0: return None
    def q(p):
        target=p*(n-1); cum=0
        for i,c in enumerate(h):
            if c and cum+c>target: return i/1000.0
            cum+=c
        return (NB-1)/1000.0
    q1,med,q3 = q(.25),q(.5),q(.75)
    iqr=q3-q1; lo_f,hi_f = q1-1.5*iqr, q3+1.5*iqr
    lo=hi=None; n_lo=n_hi=0; vmin=vmax=None
    for i,c in enumerate(h):
        if not c: continue
        v=i/1000.0
        if vmin is None: vmin=v
        vmax=v
        if v<lo_f: n_lo+=c
        elif v>hi_f: n_hi+=c
        else:
            if lo is None: lo=v
            hi=v
    if lo is None: lo=hi=med
    return dict(n=n,q1=q1,med=med,q3=q3,wlo=lo,whi=hi,nlo=n_lo,nhi=n_hi,vmin=vmin,vmax=vmax)

# Keyed per (arm, category): the plots are small multiples, one panel per
# category x sampling strategy, so each panel must describe its OWN replicate
# set. The three 15-replicate arms share their panels exactly (verified: same
# variant count to the unit); causal_power_n30000 pools 10 replicates and so
# legitimately differs.
hist=collections.defaultdict(lambda:[0]*NB)

traits=collections.defaultdict(list)
with open(O+"/dash_traits.tsv") as fh:
    for r in csv.DictReader(fh,delimiter="\t"):
        traits[(r["root"], r["rep"].rstrip("/"))].append(r)

runs=list(csv.DictReader(open(S+"/RUNS.tsv"),delimiter="\t"))
outcome=collections.Counter()
done=0
for r in runs:
    arm, rundir, cat = r["arm"], r["run_dir"], r["category"]
    wd=os.path.join(S,arm,rundir)
    vf=glob.glob(os.path.join(wd,"stage2","*","gwas_*_gtex_*_maf_*","*gwas_vars_*.tsv"))
    if not vf: continue
    pos2maf={}
    with open(vf[0]) as fh:
        for v in csv.DictReader(fh,delimiter="\t"):
            try:
                m=float(v["maf"]); sc=float(v["selco"]); p=int(float(v["position"]))
            except (ValueError,KeyError,TypeError): continue
            if m>0.5: m=1.0-m
            hist[(arm,cat,selco_bin(sc))][min(int(round(m*1000)),NB-1)]+=1
            pos2maf[p]=m
    for t in traits.get((arm,rundir),()):
        mb=maf_bin(pos2maf.get(int(t["trait"][2:])))
        if mb is None: continue
        s=float(t["self_rcp"])
        for rcp,ocol in ((50,"n_other_50"),(90,"n_other_90")):
            thr=0.5 if rcp==50 else 0.9
            o=int(t[ocol])
            if   s>thr:  oc="tp"      # colocalized with the correct trait
            elif o>0:    oc="fp"      # colocalized, but only with a wrong one
            elif s>0:    oc="uc"      # a signal at the right pairing, under the cutoff
            else:        oc="ufm"     # no signal at the right pairing at all
            outcome[(arm,t["cat"],int(t["gtex_n"]),rcp,mb,oc)]+=1
    done+=1
    if done%75==0: print(f"  {done}/{len(runs)}", flush=True)

box=[dict(arm=a,cat=c,sb=sb,**st) for (a,c,sb),h in hist.items() if (st:=stats_from_hist(h))]
stack=[dict(arm=a,cat=c,panel=p,rcp=rc,mb=mb,oc=oc,n=n) for (a,c,p,rc,mb,oc),n in outcome.items()]
assert stack and box, "aggregation produced nothing"
json.dump(dict(selco_labels=SELCO_LABELS,maf_labels=MAF_LABELS,box=box,stack=stack),
          open(O+"/dash_maf.json","w"),separators=(",",":"))
print("runs:",done," box:",len(box)," stack:",len(stack)," bytes:",os.path.getsize(O+"/dash_maf.json"))
