#!/usr/bin/env /usr/local/miniforge3/envs/hpc/bin/python
"""
Bi-Chromatic Driven LMG Bath Coupled to an Anharmonic Quantum Piston
---------------------------------------------------------------------
This script implements a quantum thermodynamics engine setup consisting of:
1. Bi-chromatic driven LMG Bath (preserves S=N/2 collective spin, dim = N+1):
   H_bath(t) = -(J/N) Sz^2 + h1 * cos(Omega1 * t) * Sx + h2 * cos(Omega2 * t) * Sy

2. Non-Linear Anharmonic Piston (Transmon / Kerr Oscillator):
   H_piston = omega_p * b†b - (U/2) * b†b†b b

3. Jaynes-Cummings Energy-Exchange Coupling:
   H_int = (g / sqrt(N)) * (b† S- + b S+)

Author: Quantum Thermodynamics & Non-Equilibrium Physics Group
Repository: DMBL-Engine
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jn_zeros
from qutip import tensor, qeye, destroy, jmat, basis, mesolve, expect

def build_bichromatic_lmg_piston(N_spin, N_piston, J, omega_p, U, g):
    """
    Construct Hamiltonian operators for the Bi-chromatic LMG + Anharmonic Piston system.
    """
    S = N_spin / 2.0
    dim_spin = int(2 * S + 1)

    # Spin operators
    Sx = jmat(S, 'x')
    Sy = jmat(S, 'y')
    Sz = jmat(S, 'z')
    Sp = jmat(S, '+')
    Sm = jmat(S, '-')
    I_spin = qeye(dim_spin)

    # Anharmonic Piston operators (Bosonic cutoff N_piston)
    b = destroy(N_piston)
    I_piston = qeye(N_piston)

    # Composite Hilbert space operators (spin ⊗ piston)
    Sx_full = tensor(Sx, I_piston)
    Sy_full = tensor(Sy, I_piston)
    Sz_full = tensor(Sz, I_piston)
    b_full  = tensor(I_spin, b)
    n_piston = tensor(I_spin, b.dag() * b)

    # Static components
    H_lmg    = -(J / N_spin) * tensor(Sz * Sz, I_piston)
    H_piston = omega_p * n_piston - (U / 2.0) * tensor(I_spin, b.dag() * b.dag() * b * b)
    H_int    = (g / np.sqrt(N_spin)) * (tensor(Sp, b) + tensor(Sm, b.dag()))

    H_static = H_lmg + H_piston + H_int

    # Dynamic drive components along x and y
    H_drive_x = Sx_full
    H_drive_y = Sy_full

    # Initial state: max-Sx eigenstate ⊗ piston vacuum |0⟩
    evals, evecs = Sx.eigenstates()
    psi0 = tensor(evecs[np.argmax(evals)], basis(N_piston, 0))

    return H_static, H_drive_x, H_drive_y, n_piston, Sz_full, psi0

# Drive coefficients for bi-chromatic driving
def drive_x_coeff(t, args):
    return args['h1'] * np.cos(args['Omega1'] * t)

def drive_y_coeff(t, args):
    return args['h2'] * np.cos(args['Omega2'] * t)

def main():
    print("=" * 70)
    print("  Bi-Chromatic Driven LMG Bath + Anharmonic Transmon Piston Engine")
    print("=" * 70)

    # Simulation Parameters
    N_spin   = 16      # Collective spin size N (S = N/2)
    N_piston = 25      # Piston Hilbert space dimension
    J        = 1.0     # LMG interaction strength
    omega_p  = 1.0     # Piston fundamental frequency
    U        = 0.15    # Piston Kerr anharmonicity
    g        = 0.5     # Jaynes-Cummings exchange coupling

    # Bi-Chromatic Drive Parameters (Incommensurate frequencies for hyper-chaos)
    Omega1 = 0.7
    Omega2 = 1.11803398875  # sqrt(5)/2
    h1     = jn_zeros(0, 5)[1] * Omega1 / 2.0
    h2     = 0.8 * h1

    T_total = 300.0
    dt      = 0.05
    tlist   = np.linspace(0, T_total, int(T_total / dt) + 1)

    print(f"Spin Size (N)       : {N_spin} (Dim = {N_spin + 1})")
    print(f"Piston Cutoff       : {N_piston}")
    print(f"Anharmonicity U     : {U}")
    print(f"Drive Frequencies   : Omega1 = {Omega1:.4f}, Omega2 = {Omega2:.4f}")
    print(f"Drive Amplitudes    : h1 = {h1:.4f}, h2 = {h2:.4f}")
    print(f"Exchange Coupling g : {g:.4f}")

    # Build system operators
    H_static, H_drive_x, H_drive_y, n_piston, Sz_op, psi0 = build_bichromatic_lmg_piston(
        N_spin, N_piston, J, omega_p, U, g
    )

    H_td = [
        H_static,
        [H_drive_x, drive_x_coeff],
        [H_drive_y, drive_y_coeff]
    ]

    args = {'h1': h1, 'h2': h2, 'Omega1': Omega1, 'Omega2': Omega2}

    print("\nRunning mesolve simulation...")
    result = mesolve(H_td, psi0, tlist, [], [n_piston, Sz_op], args=args)

    n_t  = np.array(result.expect[0])
    Sz_t = np.array(result.expect[1])

    steady_half = len(tlist) // 2
    n_steady    = n_t[steady_half:]
    mean_n      = n_steady.mean()
    std_n       = n_steady.std()
    max_n       = n_t.max()

    print("\n--- Simulation Results ---")
    print(f"Max Piston Population <n> : {max_n:.4f}")
    print(f"Steady-State Mean <n>     : {mean_n:.4f}")
    print(f"Steady-State Std Dev σ_n  : {std_n:.4f}")

    if max_n > 0.8 * N_piston:
        print(f"[Warning] Max population approaches cutoff N_piston={N_piston}. Consider increasing N_piston.")

    # Save summary plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    ax1.plot(tlist, n_t, color='crimson', lw=1.0, label=r'$\langle b^\dagger b \rangle(t)$')
    ax1.axhline(mean_n, color='black', ls='--', label=f'Steady Mean = {mean_n:.3f}')
    ax1.set_ylabel('Piston Population')
    ax1.set_title(f'Bi-Chromatic Driven LMG Bath + Transmon Piston (N={N_spin}, N_piston={N_piston}, U={U}, g={g})')
    ax1.legend(loc='upper right')
    ax1.grid(alpha=0.3)

    ax2.plot(tlist, Sz_t, color='teal', lw=1.0, label=r'$\langle S_z \rangle(t)$')
    ax2.set_xlabel('Time [1/J]')
    ax2.set_ylabel('Bath Spin Sz')
    ax2.legend(loc='upper right')
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plot_file = '/home/daneel/gitrepos/DMBL-Engine/Codes/bichromatic_piston_simulation.png'
    plt.savefig(plot_file, dpi=120)
    print(f"\nPlot saved to {plot_file}")
    print("Execution complete.")

if __name__ == '__main__':
    main()
