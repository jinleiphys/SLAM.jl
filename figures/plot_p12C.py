#!/usr/bin/env python3
"""
Publication-quality plots for p + 12C scattering comparison (Numerov vs DBMM).
PRC style: single-column width (~3.4 inches), 10pt fonts.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import os

# PRC single column width is 3.4 inches, double column is 7.0 inches
# Font size should be ~10pt to match journal text

# Set up publication-quality style for single-column figures
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif', 'Times'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.minor.width': 0.5,
    'ytick.minor.width': 0.5,
    'xtick.major.size': 4,
    'ytick.major.size': 4,
    'xtick.minor.size': 2.5,
    'ytick.minor.size': 2.5,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'mathtext.fontset': 'stix',
    'text.usetex': False,
})

# Colors (colorblind-friendly)
col_num = '#0072B2'   # Blue for Numerov
col_dbmm = '#D55E00'  # Orange for DBMM

# Data directory
data_dir = os.path.join(os.path.dirname(__file__), 'data')
fig_dir = os.path.dirname(__file__)

# System info
sys_info = r'p + $^{12}$C, $E_{\rm lab}$ = 30 MeV'

# =============================================================================
# Load data
# =============================================================================

# S-matrix data
smatrix = np.loadtxt(os.path.join(data_dir, 'smatrix.dat'))
l_vals = smatrix[:, 0].astype(int)
S_num = smatrix[:, 1] + 1j * smatrix[:, 2]
S_slam = smatrix[:, 3] + 1j * smatrix[:, 4]

# Wave function data
l_compare = [0, 2, 5]
wf_num = {}
wf_slam = {}
for l in l_compare:
    wf_num[l] = np.loadtxt(os.path.join(data_dir, f'wf_numerov_l{l}.dat'))
    wf_slam[l] = np.loadtxt(os.path.join(data_dir, f'wf_slam_l{l}.dat'))

# =============================================================================
# Plot 1: S-matrix magnitude and phase (single column)
# =============================================================================

fig1, (ax1a, ax1b) = plt.subplots(2, 1, figsize=(3.4, 4.5))
plt.subplots_adjust(hspace=0.08, left=0.18, right=0.95, top=0.92, bottom=0.10)

# |S| vs l
ax1a.plot(l_vals, np.abs(S_num), 'o-', color=col_num, markersize=5,
          linewidth=1.2, label='Numerov', markerfacecolor=col_num)
ax1a.plot(l_vals, np.abs(S_slam), 's--', color=col_dbmm, markersize=4.5,
          linewidth=1.2, label='DBMM', markerfacecolor=col_dbmm)
ax1a.set_ylabel(r'$|S_l|$')
ax1a.set_xlim(-0.5, 10.5)
ax1a.set_ylim(0, 1.12)
ax1a.set_xticks(np.arange(0, 11, 2))
ax1a.set_yticks(np.arange(0, 1.1, 0.2))
ax1a.xaxis.set_minor_locator(MultipleLocator(1))
ax1a.yaxis.set_minor_locator(MultipleLocator(0.1))
ax1a.legend(loc='lower right', frameon=True, fancybox=False,
            edgecolor='black', framealpha=1)
ax1a.set_title(sys_info, fontsize=11)
ax1a.text(0.03, 0.88, '(a)', transform=ax1a.transAxes, fontsize=11, fontweight='bold')
ax1a.tick_params(labelbottom=False)

# arg(S) vs l
phases_num = np.angle(S_num) * 180 / np.pi
phases_slam = np.angle(S_slam) * 180 / np.pi

ax1b.plot(l_vals, phases_num, 'o-', color=col_num, markersize=5,
          linewidth=1.2, label='Numerov', markerfacecolor=col_num)
ax1b.plot(l_vals, phases_slam, 's--', color=col_dbmm, markersize=4.5,
          linewidth=1.2, label='DBMM', markerfacecolor=col_dbmm)
ax1b.set_xlabel(r'$l$')
ax1b.set_ylabel(r'$\arg(S_l)$ [deg]')
ax1b.set_xlim(-0.5, 10.5)
ax1b.set_xticks(np.arange(0, 11, 2))
ax1b.xaxis.set_minor_locator(MultipleLocator(1))
ax1b.yaxis.set_minor_locator(AutoMinorLocator())
ax1b.text(0.03, 0.88, '(b)', transform=ax1b.transAxes, fontsize=11, fontweight='bold')

fig1.savefig(os.path.join(fig_dir, 'p12C_Smatrix.pdf'), dpi=300, bbox_inches='tight')
fig1.savefig(os.path.join(fig_dir, 'p12C_Smatrix.eps'), dpi=300, bbox_inches='tight')
plt.close(fig1)
print('Saved: p12C_Smatrix.pdf and .eps')

# =============================================================================
# Plot 2: Argand diagram (single column)
# =============================================================================

fig2, ax2 = plt.subplots(1, 1, figsize=(3.4, 3.0))
plt.subplots_adjust(left=0.16, right=0.95, top=0.90, bottom=0.14)

# Unit circle
theta = np.linspace(0, 2*np.pi, 200)
ax2.plot(np.cos(theta), np.sin(theta), '--', color='gray', linewidth=1.0, zorder=1, label=r'$|S|=1$')

# S-matrix trajectories
ax2.plot(S_num.real, S_num.imag, 'o-', color=col_num, markersize=5,
         linewidth=1.2, label='Numerov', markerfacecolor=col_num, zorder=3)
ax2.plot(S_slam.real, S_slam.imag, 's--', color=col_dbmm, markersize=4.5,
         linewidth=1.2, label='DBMM', markerfacecolor=col_dbmm, zorder=2)

# Use blue for all l labels
l_label_color = '#0072B2'

# Add l labels - optimized positions to avoid overlap
label_pos = {
    0: (0.06, -0.10),
    1: (-0.15, 0.00),
    2: (-0.12, 0.07),
    3: (0.06, 0.07),
    4: (0.06, 0.06),
    5: (0.06, 0.05),
    6: (0.05, 0.05),
    8: (0.05, 0.04),
    10: (0.05, -0.06),
}
for i, l in enumerate(l_vals):
    if l in label_pos:
        ox, oy = label_pos[l]
        ax2.annotate(str(l), (S_slam[i].real + ox, S_slam[i].imag + oy),
                     fontsize=9, ha='left', va='center', fontweight='bold',
                     color=l_label_color)

ax2.set_xlabel(r'Re$(S_l)$')
ax2.set_ylabel(r'Im$(S_l)$')
ax2.set_xlim(-0.7, 1.15)
ax2.set_ylim(-0.6, 0.45)
ax2.set_aspect('equal', adjustable='box')
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())
# Position legend in center-right area to avoid blocking the unit circle
ax2.legend(loc='center right', bbox_to_anchor=(0.75, 0.25), frameon=True, fancybox=False,
           edgecolor='black', framealpha=1)
ax2.set_title(sys_info, fontsize=11)

fig2.savefig(os.path.join(fig_dir, 'p12C_argand.pdf'), dpi=300, bbox_inches='tight')
fig2.savefig(os.path.join(fig_dir, 'p12C_argand.eps'), dpi=300, bbox_inches='tight')
plt.close(fig2)
print('Saved: p12C_argand.pdf and .eps')

# =============================================================================
# Plot 3: Wave functions (single column width for 2x2 panels)
# =============================================================================

l_compare_reduced = [0, 5]
fig3, axes3 = plt.subplots(2, 2, figsize=(3.4, 3.2))
plt.subplots_adjust(hspace=0.12, wspace=0.35, left=0.15, right=0.97, top=0.90, bottom=0.12)

# Add centered title at top
fig3.suptitle(sys_info, fontsize=11, fontweight='bold')

for i, l in enumerate(l_compare_reduced):
    # Get data
    r_num = wf_num[l][:, 0]
    psi_num_re = wf_num[l][:, 1]
    psi_num_im = wf_num[l][:, 2]

    r_slam = wf_slam[l][:, 0]
    psi_slam_re = wf_slam[l][:, 1]
    psi_slam_im = wf_slam[l][:, 2]

    panel_labels = ['(a)', '(b)', '(c)', '(d)']

    # Real part
    ax_re = axes3[i, 0]
    ax_re.plot(r_num, psi_num_re, '-', color=col_num, linewidth=0.8,
               label='Numerov' if i == 0 else '')
    ax_re.scatter(r_slam, psi_slam_re, s=8, color=col_dbmm, marker='o',
                  label='DBMM' if i == 0 else '', zorder=5, edgecolors='none')
    ax_re.axhline(0, color='gray', linewidth=0.4, linestyle='--', zorder=1)
    ax_re.set_xlim(0, 15)
    ax_re.set_ylabel(rf'Re[$\psi_{l}(r)$]', fontsize=9)
    ax_re.tick_params(axis='both', labelsize=8)
    ax_re.xaxis.set_minor_locator(AutoMinorLocator())
    ax_re.yaxis.set_minor_locator(AutoMinorLocator())
    ax_re.text(0.05, 0.88, panel_labels[i*2], transform=ax_re.transAxes,
               fontsize=9, fontweight='bold')
    if i == 1:
        ax_re.set_xlabel(r'$r$ [fm]', fontsize=9)
    else:
        ax_re.tick_params(labelbottom=False)

    # Imaginary part
    ax_im = axes3[i, 1]
    ax_im.plot(r_num, psi_num_im, '-', color=col_num, linewidth=0.8,
               label='Numerov' if i == 1 else '')
    ax_im.scatter(r_slam, psi_slam_im, s=8, color=col_dbmm, marker='o',
                  zorder=5, edgecolors='none', label='DBMM' if i == 1 else '')
    ax_im.axhline(0, color='gray', linewidth=0.4, linestyle='--', zorder=1)
    ax_im.set_xlim(0, 15)
    ax_im.set_ylabel(rf'Im[$\psi_{l}(r)$]', fontsize=9)
    ax_im.tick_params(axis='both', labelsize=8)
    ax_im.xaxis.set_minor_locator(AutoMinorLocator())
    ax_im.yaxis.set_minor_locator(AutoMinorLocator())
    ax_im.text(0.05, 0.88, panel_labels[i*2+1], transform=ax_im.transAxes,
               fontsize=9, fontweight='bold')
    if i == 1:
        ax_im.legend(loc='lower left', frameon=True, fancybox=False,
                     edgecolor='black', framealpha=1, fontsize=7)
        ax_im.set_xlabel(r'$r$ [fm]', fontsize=9)
    else:
        ax_im.tick_params(labelbottom=False)

fig3.savefig(os.path.join(fig_dir, 'p12C_wavefunction.pdf'), dpi=300, bbox_inches='tight')
fig3.savefig(os.path.join(fig_dir, 'p12C_wavefunction.eps'), dpi=300, bbox_inches='tight')
plt.close(fig3)
print('Saved: p12C_wavefunction.pdf and .eps')

print('\nAll figures generated successfully!')
