#!/usr/bin/env python3
"""
Publication-quality plots for p + 12C scattering comparison (Numerov vs SLAM).
PRC style: boxed frames, proper fonts, no grids.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import os

# Set up publication-quality style - VERY LARGE fonts for publication
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif', 'Times'],
    'font.size': 22,
    'axes.labelsize': 24,
    'axes.titlesize': 24,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'legend.fontsize': 20,
    'axes.linewidth': 1.8,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'xtick.minor.width': 1.0,
    'ytick.minor.width': 1.0,
    'xtick.major.size': 8,
    'ytick.major.size': 8,
    'xtick.minor.size': 5,
    'ytick.minor.size': 5,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'mathtext.fontset': 'stix',
    'text.usetex': False,
})

# Colors (colorblind-friendly)
col_num = '#0072B2'   # Blue for Numerov
col_slam = '#D55E00'  # Orange for SLAM
col_green = '#009E73'
col_gold = '#E69F00'

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
# Plot 1: S-matrix magnitude and phase
# =============================================================================

fig1, (ax1a, ax1b) = plt.subplots(2, 1, figsize=(8, 10))
plt.subplots_adjust(hspace=0.10, left=0.18, right=0.95, top=0.94, bottom=0.08)

# |S| vs l
ax1a.plot(l_vals, np.abs(S_num), 'o-', color=col_num, markersize=12,
          linewidth=2.5, label='Numerov', markerfacecolor=col_num)
ax1a.plot(l_vals, np.abs(S_slam), 's--', color=col_slam, markersize=11,
          linewidth=2.5, label='SLAM', markerfacecolor=col_slam)
ax1a.set_ylabel(r'$|S_l|$')
ax1a.set_xlim(-0.5, 10.5)
ax1a.set_ylim(0, 1.12)
ax1a.set_xticks(np.arange(0, 11, 2))
ax1a.set_yticks(np.arange(0, 1.1, 0.2))
ax1a.xaxis.set_minor_locator(MultipleLocator(1))
ax1a.yaxis.set_minor_locator(MultipleLocator(0.1))
ax1a.legend(loc='lower right', frameon=True, fancybox=False,
            edgecolor='black', framealpha=1)
ax1a.set_title(sys_info, fontsize=26)
ax1a.text(0.03, 0.86, '(a)', transform=ax1a.transAxes, fontsize=26, fontweight='bold')
# Remove x tick labels for top panel (shared axis)
ax1a.tick_params(labelbottom=False)

# arg(S) vs l
phases_num = np.angle(S_num) * 180 / np.pi
phases_slam = np.angle(S_slam) * 180 / np.pi

ax1b.plot(l_vals, phases_num, 'o-', color=col_num, markersize=12,
          linewidth=2.5, label='Numerov', markerfacecolor=col_num)
ax1b.plot(l_vals, phases_slam, 's--', color=col_slam, markersize=11,
          linewidth=2.5, label='SLAM', markerfacecolor=col_slam)
ax1b.set_xlabel(r'$l$')
ax1b.set_ylabel(r'$\arg(S_l)$ [deg]')
ax1b.set_xlim(-0.5, 10.5)
ax1b.set_xticks(np.arange(0, 11, 2))
ax1b.xaxis.set_minor_locator(MultipleLocator(1))
ax1b.yaxis.set_minor_locator(AutoMinorLocator())
# No legend in (b) - already shown in (a)
ax1b.text(0.03, 0.86, '(b)', transform=ax1b.transAxes, fontsize=26, fontweight='bold')

fig1.savefig(os.path.join(fig_dir, 'p12C_Smatrix.pdf'), dpi=300, bbox_inches='tight')
fig1.savefig(os.path.join(fig_dir, 'p12C_Smatrix.eps'), dpi=300, bbox_inches='tight')
plt.close(fig1)
print('Saved: p12C_Smatrix.pdf and .eps')

# =============================================================================
# Plot 2: Wave functions
# =============================================================================

# Use MUCH larger fonts specifically for this multi-panel figure
fig2, axes2 = plt.subplots(3, 2, figsize=(20, 20))
plt.subplots_adjust(hspace=0.10, wspace=0.35, left=0.10, right=0.97, top=0.95, bottom=0.06)

# Double the font sizes for this figure
title_fs = 48
label_fs = 44
tick_fs = 38
legend_fs = 38

for i, l in enumerate(l_compare):
    # Get data
    r_num = wf_num[l][:, 0]
    psi_num_re = wf_num[l][:, 1]
    psi_num_im = wf_num[l][:, 2]

    r_slam = wf_slam[l][:, 0]
    psi_slam_re = wf_slam[l][:, 1]
    psi_slam_im = wf_slam[l][:, 2]

    # Real part
    ax_re = axes2[i, 0]
    ax_re.plot(r_num, psi_num_re, '-', color=col_num, linewidth=3.5,
               label='Numerov' if i == 0 else '')
    ax_re.scatter(r_slam, psi_slam_re, s=120, color=col_slam, marker='o',
                  label='SLAM' if i == 0 else '', zorder=5, edgecolors='none')
    ax_re.axhline(0, color='gray', linewidth=0.8, linestyle='--', zorder=1)
    ax_re.set_xlim(0, 15)
    ax_re.set_ylabel(rf'Re[$\psi_{l}(r)$]', fontsize=label_fs)
    ax_re.tick_params(axis='both', labelsize=tick_fs)
    ax_re.xaxis.set_minor_locator(AutoMinorLocator())
    ax_re.yaxis.set_minor_locator(AutoMinorLocator())
    if i == 0:
        ax_re.legend(loc='upper right', frameon=True, fancybox=False,
                     edgecolor='black', framealpha=1, fontsize=legend_fs)
        ax_re.set_title(sys_info + r' $-$ Real part', fontsize=title_fs)
    # Only show x label on bottom row
    if i == 2:
        ax_re.set_xlabel(r'$r$ [fm]', fontsize=label_fs)
    else:
        ax_re.tick_params(labelbottom=False)

    # Imaginary part
    ax_im = axes2[i, 1]
    ax_im.plot(r_num, psi_num_im, '-', color=col_num, linewidth=3.5)
    ax_im.scatter(r_slam, psi_slam_im, s=120, color=col_slam, marker='o',
                  zorder=5, edgecolors='none')
    ax_im.axhline(0, color='gray', linewidth=0.8, linestyle='--', zorder=1)
    ax_im.set_xlim(0, 15)
    ax_im.set_ylabel(rf'Im[$\psi_{l}(r)$]', fontsize=label_fs)
    ax_im.tick_params(axis='both', labelsize=tick_fs)
    ax_im.xaxis.set_minor_locator(AutoMinorLocator())
    ax_im.yaxis.set_minor_locator(AutoMinorLocator())
    if i == 0:
        ax_im.set_title(sys_info + r' $-$ Imag part', fontsize=title_fs)
    # Only show x label on bottom row
    if i == 2:
        ax_im.set_xlabel(r'$r$ [fm]', fontsize=label_fs)
    else:
        ax_im.tick_params(labelbottom=False)

fig2.savefig(os.path.join(fig_dir, 'p12C_wavefunction.pdf'), dpi=300, bbox_inches='tight')
fig2.savefig(os.path.join(fig_dir, 'p12C_wavefunction.eps'), dpi=300, bbox_inches='tight')
plt.close(fig2)
print('Saved: p12C_wavefunction.pdf and .eps')


# =============================================================================
# Plot 4: Argand diagram
# =============================================================================

fig4, ax4 = plt.subplots(1, 1, figsize=(9, 8))
plt.subplots_adjust(left=0.14, right=0.95, top=0.92, bottom=0.12)

# Unit circle
theta = np.linspace(0, 2*np.pi, 200)
ax4.plot(np.cos(theta), np.sin(theta), '--', color='gray', linewidth=2, zorder=1, label=r'$|S|=1$')

# S-matrix trajectories
ax4.plot(S_num.real, S_num.imag, 'o-', color=col_num, markersize=13,
         linewidth=2.5, label='Numerov', markerfacecolor=col_num, zorder=3)
ax4.plot(S_slam.real, S_slam.imag, 's--', color=col_slam, markersize=12,
         linewidth=2.5, label='SLAM', markerfacecolor=col_slam, zorder=2)

# Add l labels - carefully positioned
label_pos = {
    0: (0.08, -0.10),
    1: (-0.18, 0.05),
    2: (-0.16, 0.08),
    3: (0.08, 0.09),
    4: (0.08, 0.09),
    5: (0.08, 0.07),
    6: (0.07, 0.07),
    8: (0.07, 0.06),
    10: (0.07, -0.07),
}
for i, l in enumerate(l_vals):
    if l in label_pos:
        ox, oy = label_pos[l]
        ax4.annotate(str(l), (S_slam[i].real + ox, S_slam[i].imag + oy),
                     fontsize=20, ha='left', va='center')

ax4.set_xlabel(r'Re$(S_l)$')
ax4.set_ylabel(r'Im$(S_l)$')
ax4.set_xlim(-0.7, 1.15)
ax4.set_ylim(-0.6, 0.45)
ax4.set_aspect('equal', adjustable='box')
ax4.xaxis.set_minor_locator(AutoMinorLocator())
ax4.yaxis.set_minor_locator(AutoMinorLocator())
# Move legend to upper left to avoid blocking data
ax4.legend(loc='upper left', frameon=True, fancybox=False,
           edgecolor='black', framealpha=1)
ax4.set_title(sys_info, fontsize=26)

fig4.savefig(os.path.join(fig_dir, 'p12C_argand.pdf'), dpi=300, bbox_inches='tight')
fig4.savefig(os.path.join(fig_dir, 'p12C_argand.eps'), dpi=300, bbox_inches='tight')
plt.close(fig4)
print('Saved: p12C_argand.pdf and .eps')

print('\nAll figures generated successfully!')
