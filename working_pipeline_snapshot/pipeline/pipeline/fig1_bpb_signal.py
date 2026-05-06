"""
fig1_bpb_signal.py  --  BPB signal structure on 38 active SPARC galaxies
Data: data/sparc_active_sample.csv  (override: pass CSV path as argv[1])
Annotations: paper values hardcoded (Cohen d=1.44, thresholds, p-values)
Output: fig1a_signal_hist.pdf, fig1b_decision_space.pdf
"""
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

T_dR1=-0.131; T_s=0.320
C_POS="#4C72B0"; C_NEG="#DD8452"; FS=9

csv = sys.argv[1] if len(sys.argv)>1 else os.path.join(os.path.dirname(os.path.abspath(__file__)),"data","sparc_active_sample.csv")
df = pd.read_csv(csv)
pos = df[df["Areq"]>0]; neg = df[df["Areq"]<0]
dR1_pos,sinn_pos = pos["dR1"].values, pos["sinner"].values
dR1_neg,sinn_neg = neg["dR1"].values, neg["sinner"].values
ddo = df[df["name"]=="DDO161"]
d1,s1 = (float(ddo["dR1"].iloc[0]),float(ddo["sinner"].iloc[0])) if len(ddo) else (-0.090,0.310)
print(f"Loaded {len(df)} rows: {len(pos)} Areq>0, {len(neg)} Areq<0")

# Panel (a)
fig,ax=plt.subplots(figsize=(3.4,2.8))
bins=np.arange(-0.265,0.025,0.020)
ax.hist(dR1_pos,bins=bins,color=C_POS,alpha=0.80,label=rf"$A_{{\rm req}}>0$ ($n={len(pos)}$)")
ax.hist(dR1_neg,bins=bins,color=C_NEG,alpha=0.80,label=rf"$A_{{\rm req}}<0$ ($n={len(neg)}$)")
ax.axvline(T_dR1,color="k",ls="--",lw=1.4,label=rf"$T_{{\delta R_1}}={T_dR1}$")
ax.set_xlabel(r"$\delta R_1$",fontsize=FS); ax.set_ylabel("Count",fontsize=FS)
ax.text(0.97,0.95,"Cohen $d=1.44$,\n$p=1.6\\times10^{-5}$",transform=ax.transAxes,ha="right",va="top",fontsize=7)
ax.legend(fontsize=7,loc="upper left"); ax.tick_params(labelsize=7)
fig.tight_layout(); fig.savefig("fig1a_signal_hist.pdf",dpi=300,bbox_inches="tight"); plt.close(fig)
print("Saved: fig1a_signal_hist.pdf")

# Panel (b)
fig,ax=plt.subplots(figsize=(3.4,2.8))
ax.scatter(dR1_pos,sinn_pos,c=C_POS,s=18,zorder=3,label=r"$A_{\rm req}>0$")
ax.scatter(dR1_neg,sinn_neg,c=C_NEG,s=18,zorder=3,label=r"$A_{\rm req}<0$")
ax.scatter(d1,s1,marker="*",s=90,c=C_NEG,zorder=4,edgecolors="k",linewidths=0.5,label="DDO161")
ax.axvline(T_dR1,color="k",ls="--",lw=1.4); ax.axhline(T_s,color="k",ls="--",lw=1.4)
ax.text(T_dR1-0.003,0.455,r"$T_{\delta R_1}$",ha="right",va="top",fontsize=7)
ax.text(-0.245,T_s+0.005,r"$T_s$",ha="left",va="bottom",fontsize=7)
ax.set_xlabel(r"$\delta R_1$",fontsize=FS); ax.set_ylabel(r"$s_{\rm inner}$",fontsize=FS)
ax.legend(fontsize=7,loc="upper right"); ax.tick_params(labelsize=7)
fig.tight_layout(); fig.savefig("fig1b_decision_space.pdf",dpi=300,bbox_inches="tight"); plt.close(fig)
print("Saved: fig1b_decision_space.pdf")
