import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.family'] = 'DejaVu Sans'

# CREST cross-tabulation (Nissen 2021, Table 6), Rule 1: eligible = clinician "OK to use"
TP_app, FP_app = 23, 3          # application approvals among 500
TP_all, FP_all = 39, 461        # unscreened access: everyone approved
N = 500

x = np.linspace(0.0005, 0.995, 4000)
odds = x / (1 - x)
nb_app  = TP_app/N - (FP_app/N)*odds
nb_all  = TP_all/N - (FP_all/N)*odds

# crossovers
p_cross_all = 16/474            # app vs unscreened: 3.3755%
p_cross_none = 23/26            # app vs no access: 88.46% (= PPV)

fig, ax = plt.subplots(figsize=(8.6, 5.4), dpi=150)
ax.axvspan(5, 20, color='0.92', zorder=0)
ax.text(12.5, 0.093, 'Reasonable range\n(5\u201320%)', ha='center', va='top',
        fontsize=9.5, style='italic', color='0.25')

ax.plot(x*100, nb_app, color='black', lw=2.0, label='Application', zorder=3)
ax.plot(x*100, nb_all, color='0.45', lw=1.8, ls='--', label='Access without screening', zorder=2)
ax.axhline(0, color='0.2', lw=1.0, label='No access', zorder=1)

for px, py, lbl, dx, dy in [(p_cross_all*100, TP_app/N-(FP_app/N)*(p_cross_all/(1-p_cross_all)), '3.4%', 2.5, 0.012),
                            (p_cross_none*100, 0, '88.5%', -1.5, 0.012)]:
    ax.plot(px, py, 'o', mfc='white', mec='black', ms=6, zorder=4)
    ax.annotate(lbl, (px, py), xytext=(px+dx, py+dy), fontsize=9.5, ha='left')

ax.set_xlim(0, 100); ax.set_ylim(-0.03, 0.10)
ax.set_xlabel('Threshold probability (%)', fontsize=11)
ax.set_ylabel('Net benefit', fontsize=11)
ax.spines[['top', 'right']].set_visible(False)
ax.legend(frameon=False, fontsize=10, loc='upper right', bbox_to_anchor=(1.0, 0.97))
plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/figure4_decision_curves.png', dpi=300, bbox_inches='tight')
print("crossovers:", round(p_cross_all*100, 2), round(p_cross_none*100, 2))
print("NB app at 5% and 20%:", round(TP_app/N-(FP_app/N)*(0.05/0.95), 4), round(TP_app/N-(FP_app/N)*0.25, 4))
