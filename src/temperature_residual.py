#!/usr/bin/env python3
"""Temperature-controlled snap-rate residual at UAP incident times.

Snapping-shrimp snap rate is dominated by temperature. So the meaningful
quantity for a UAP search is the RESIDUAL: observed snap rate minus what
in-situ temperature predicts. A residual rise coincident with a UAP event
is a snap-rate change temperature does NOT explain.

Method (CI01, Channel Islands):
  1. per-minute snap rate for every CI01 file (cached to npz);
  2. align the RBR in-situ temperature (5 s) to each minute;
  3. fit snap_rate = f(temperature) globally (binned-median lookup) from ALL
     pooled data -> the temperature-expected snap rate;
  4. residual = observed - expected;
  5. for each UAP incident night, residual after/before change (+/-20 min)
     at 06:00 UTC, vs a control-day null of residual changes.
"""
import glob, re, datetime, json, os
import numpy as np
import soundfile as sf
import scipy.io as sio
from scipy import signal
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

CB=(2000.0,20000.0); BIN=60; HALF=20
INCIDENTS={"2019-07-15":"Jul14 (Kidd)","2019-07-16":"Jul15 (Omaha)",
           "2019-07-26":"Jul25","2019-07-31":"Jul30"}


def dn2dt(dn): return datetime.datetime.fromordinal(int(dn)-366)+datetime.timedelta(days=dn%1)
def dt2dn(dt): return dt.toordinal()+366+(dt.hour*3600+dt.minute*60+dt.second)/86400.0
def fstart(p):
    m=re.search(r"_(\d{6})(\d{6})\.flac$",p); return datetime.datetime.strptime(m.group(1)+m.group(2),"%y%m%d%H%M%S")


def snap_rate(path):
    cache=path+".snaprate.npz"
    if os.path.exists(cache):
        z=np.load(cache); return z["f0off"].item(), z["rate"]
    info=sf.info(path); fs=info.samplerate; n=BIN*fs
    sos=signal.butter(4,CB,btype="bandpass",fs=fs,output="sos")
    rate=[]
    with sf.SoundFile(path) as fh:
        for off in range(0,info.frames-n,n):
            fh.seek(off); x=fh.read(n,dtype="float64")
            if x.ndim>1: x=x[:,0]
            xb=signal.sosfiltfilt(sos,x); env=np.abs(signal.hilbert(xb))
            med=np.median(env); mad=np.median(np.abs(env-med))+1e-30
            rate.append(len(signal.find_peaks(env,height=med+8*mad,distance=int(0.01*fs))[0]))
    rate=np.array(rate,float)
    np.savez(cache,f0off=fstart(path).timestamp(),rate=rate)
    return fstart(path), rate


def main():
    d=sio.loadmat("data/omaha/CI01_02_temperature.mat",struct_as_record=False,squeeze_me=True)["RSK"].data
    ts=np.asarray(d.tstamp).ravel().astype(float); tp=np.asarray(d.values).ravel().astype(float)
    good=(tp>5)&(tp<30)  # drop out-of-water artifacts
    ts,tp=ts[good],tp[good]

    files=(sorted(glob.glob("data/omaha/SanctSound_CI01*.flac"))
           +sorted(glob.glob("data/omaha/inc_*CI01*.flac"))
           +sorted(glob.glob("data/omaha/ctrl_*CI01*.flac")))
    recs=[]  # (datestr, f0, minutes_arr, rate_arr, temp_arr)
    for p in files:
        f0=fstart(p) if not p.split("/")[-1].startswith(("inc_","ctrl_")) else fstart(re.sub(r".*_(\d{6})(\d{6})","\\1\\2",p) if False else p)
        f0=fstart(p)
        _,rate=snap_rate(p)
        mins=np.arange(len(rate))
        dn=dt2dn(f0)+mins/1440.0
        temp=np.interp(dn,ts,tp)
        recs.append((f0.strftime("%Y-%m-%d"),f0,mins,rate,temp))

    # global temperature lookup f(temp) from pooled data
    allt=np.concatenate([r[4] for r in recs]); allr=np.concatenate([r[3] for r in recs])
    bins=np.linspace(np.percentile(allt,1),np.percentile(allt,99),25)
    idx=np.digitize(allt,bins)
    fexp=np.array([np.median(allr[idx==i]) if np.any(idx==i) else np.nan for i in range(len(bins)+1)])
    def expected(temp):
        i=np.clip(np.digitize(temp,bins),0,len(fexp)-1); return fexp[i]
    r_global=np.corrcoef(allt,allr)[0,1]
    print(f"pooled CI01 snap-rate vs temperature: r = {r_global:+.3f} "
          f"(n={len(allt)} min over {len(recs)} files)")

    def resid_series(rate,temp): return rate-expected(temp)
    half=HALF  # minutes == bins (BIN=60s)
    # control null of residual after/before changes (control days only)
    null=[]
    ctrl_dates=set()
    for datestr,f0,mins,rate,temp in recs:
        if any(datestr==k for k in INCIDENTS): continue
        ctrl_dates.add(datestr)
        res=resid_series(rate,temp)
        for c in range(half,len(res)-half):
            b=res[c-half:c].mean(); a=res[c:c+half].mean()
            null.append(a-b)
    null=np.array(null)

    rows=[]
    for datestr,label in INCIDENTS.items():
        rec=next((r for r in recs if r[0]==datestr),None)
        if rec is None:
            rows.append({"date":datestr,"label":label,"status":"no data"}); continue
        _,f0,mins,rate,temp=rec
        res=resid_series(rate,temp)
        c=int((datetime.datetime.fromisoformat(datestr+"T06:00:00")-f0).total_seconds()/60)
        if not (half<=c<len(res)-half):
            rows.append({"date":datestr,"label":label,"status":"no window"}); continue
        raw_b,raw_a=rate[c-half:c].mean(),rate[c:c+half].mean()
        t_b,t_a=temp[c-half:c].mean(),temp[c:c+half].mean()
        res_change=res[c:c+half].mean()-res[c-half:c].mean()
        p=float(np.mean(np.array(null)>=res_change))
        rows.append({"date":datestr,"label":label,"status":"ok",
                     "raw_pct":round((raw_a/raw_b-1)*100,1),
                     "temp_change_C":round(t_a-t_b,3),
                     "temp_expl_pct":round((expected(np.array([t_a]))[0]/expected(np.array([t_b]))[0]-1)*100,1),
                     "residual_snaps_per_min":round(res_change,1),
                     "p_vs_control":round(p,3)})
    json.dump({"r_global":r_global,"n_null":len(null),"incidents":rows},
              open("results/temperature_residual.json","w"),indent=1)
    print(f"\ntemperature-controlled residual at UAP incidents (null n={len(null)}):")
    print(f"{'date':12s}{'label':14s}{'raw%':>7s}{'dTempC':>8s}{'tempExpl%':>10s}{'residual/min':>13s}{'p':>7s}")
    for r in rows:
        if r["status"]!="ok":
            print(f"{r['date']:12s}{r['label']:14s}  {r['status']}"); continue
        print(f"{r['date']:12s}{r['label']:14s}{r['raw_pct']:+7.1f}{r['temp_change_C']:+8.3f}"
              f"{r['temp_expl_pct']:+10.1f}{r['residual_snaps_per_min']:+13.1f}{r['p_vs_control']:7.3f}")
    print("\nresidual = snap-rate change NOT explained by temperature; low p = "
          "temperature-independent rise at that incident")


if __name__=="__main__":
    main()
