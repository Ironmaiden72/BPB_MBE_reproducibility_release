"""
fig2_areq_scatter.py  --  Predicted vs true Areq for 38 active galaxies
Data: data/sparc_active_sample.csv
Annotation values: r=0.657 (in-sample), LOO r=0.559, r=0.879 (oracle) -- paper values
Output: fig2a_areq_physical.pdf, fig2b_areq_oracle.pdf
"""
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

C_OK="#4C72B0"; C_ERR="#DD8452"; FS=9; LW=1.2

csv = sys.argv[1] if len(sys.argv)>1 else os.path.join(os.path.dirname(os.path.abspath(__file__)),"data","sparc_active_sample.csv")
df = pd.read_csv(csv)
At = df["Areq"].values
Ap_phys = df["Areq_pred_physical"].values
Ap_orac = df["Areq_pred_oracle"].values
ok_phys = df["sign_correct_physical"].values.astype(bool)

print(f"Loaded {len(df)} rows. Sign correct: {ok_phys.sum()}/38")

lim=(-2.5,4.5)

def panel(ax, At, Ap, ok, r_label, loor_label, mae_label, title):
    cols = [C_OK if o else C_ERR for o in ok]
    ax.scatter(At, Ap, c=cols, s=16, zorder=3)
    ax.plot(lim, lim, "k--", lw=LW, zorder=2)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel(r"$A_{\rm req}^{\rm true}$", fontsize=FS)
    ax.set_ylabel(r"$A_{\rm req}^{\rm pred}$", fontsize=FS)
    # Paper annotation values
    ax.text(0.05,0.95, f"$r={r_label}$\nMAE$={mae_label}$",
            transform=ax.transAxes, va="top", fontsize=7)
    ax.text(0.50,0.05, f"in-sample $r={r_label}$",
            transform=ax.transAxes, ha="center", fontsize=7, style="italic")
    ax.set_title(title, fontsize=8)
    ax.tick_params(labelsize=7)

fig,ax=plt.subplots(figsize=(3.4,3.0))
panel(ax, At, Ap_phys, ok_phys,
      r_label="0.657", loor_label="0.559", mae_label="0.859",
      title="Physical rule (LOO $r=0.559$)")
fig.tight_layout(); fig.savefig("fig2a_areq_physical.pdf",dpi=300,bbox_inches="tight"); plt.close(fig)
print("Saved: fig2a_areq_physical.pdf")

fig,ax=plt.subplots(figsize=(3.4,3.0))
panel(ax, At, Ap_orac, np.ones(len(At),dtype=bool),
      r_label="0.879", loor_label="0.879", mae_label="0.727",
      title="Oracle sign ($r=0.879$)")
fig.tight_layout(); fig.savefig("fig2b_areq_oracle.pdf",dpi=300,bbox_inches="tight"); plt.close(fig)
print("Saved: fig2b_areq_oracle.pdf")
