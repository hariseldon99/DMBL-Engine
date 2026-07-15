import multiprocessing
import os

import numpy as np
from qutip import basis, destroy, jmat, mesolve, qeye, tensor
from scipy.special import jn_zeros
from tqdm.auto import tqdm


# ---------------------------------------------------------------------
# Self-contained helpers
# ---------------------------------------------------------------------
def build_ops(N, N_ph, J, omega_0, g):
    """Build operators for the LMG-cavity system."""
    S = N / 2.0
    dim_spin = int(2 * S + 1)

    Sx = jmat(S, 'x')
    Sz = jmat(S, 'z')
    I_spin = qeye(dim_spin)

    a = destroy(N_ph)
    I_ph = qeye(N_ph)

    Sx_full = tensor(Sx, I_ph)
    Sz_full = tensor(Sz, I_ph)
    n_ph_full = tensor(I_spin, a.dag() * a)
    X_ph_full = tensor(I_spin, a + a.dag())

    H_lmg = -(J / N) * tensor(Sz * Sz, I_ph)
    H_cav = omega_0 * n_ph_full
    H_int = (2.0 * g / np.sqrt(N)) * X_ph_full * Sz_full

    H_static = H_lmg + H_cav + H_int
    H_drive = Sx_full

    evals, evecs = Sx.eigenstates()
    psi0 = tensor(evecs[np.argmax(evals)], basis(N_ph, 0))

    return H_static, H_drive, n_ph_full, psi0


def drive_coeff(t, args):
    return args['h0'] * np.cos(args['Omega'] * t)


def cavity_population_std(task):
    """Run mesolve for one spin size and return std[n(t)]."""
    N_spin, N_ph, J, omega_0, g, h0, Omega, T_total, dt = task
    tlist_sweep = np.linspace(0, T_total, int(T_total / dt) + 1)

    H_static_i, H_drive_i, n_ph_op_i, psi0_i = build_ops(N_spin, N_ph, J, omega_0, g)
    res = mesolve(
        [H_static_i, [H_drive_i, drive_coeff]],
        psi0_i,
        tlist_sweep,
        [],
        [n_ph_op_i],
        args={'h0': h0, 'Omega': Omega},
    )

    n_t = np.asarray(res.expect[0])
    return N_ph, n_t.std()


# ---------------------------------------------------------------------
# Sweep settings: chaotic regime, fixed cavity truncation
# ---------------------------------------------------------------------
N_spin = 400
J = 1.0
omega_0 = 1.0
g = 0.10
T_total = 3000.0
dt = 0.05

# Keep the cavity-driving parameters in the chaotic regime.
Omega_sweep = 0.7
h0_sweep = jn_zeros(0, 5)[1] * Omega_sweep / 2

N_values = np.arange(2, 401, 32)
tasks = [
    (N_spin, N_ph, J, omega_0, g, h0_sweep, Omega_sweep, T_total, dt)
    for N_ph in N_values
]

# Use multiprocessing Pool to parallelize the sweep.
# Check PBS_NCPUS environment variable for HPC; otherwise use local CPU count
n_cpus = int(os.environ.get('PBS_NCPUS', multiprocessing.cpu_count()))
ctx = multiprocessing.get_context('fork')
with ctx.Pool(processes=min(len(tasks), n_cpus)) as pool:
    std_results = list(
        tqdm(
            pool.imap_unordered(cavity_population_std, tasks),
            total=len(tasks),
            desc='Sweeping N',
        )
    )

std_results.sort(key=lambda item: item[0])
N_out = np.array([item[0] for item in std_results])
std_n = np.array([item[1] for item in std_results])

checkpoint_file = 'cavity_population_sweep_checkpoint.npz'
np.savez_compressed(
    checkpoint_file,
    N_out=N_out,
    std_n=std_n,
    N_spin=N_spin,
    J=J,
    omega_0=omega_0,
    g=g,
    Omega_sweep=Omega_sweep,
    h0_sweep=h0_sweep,
    T_total=T_total,
    dt=dt,
)

print(f'Saved sweep checkpoint to {checkpoint_file}')
print(f'Used {n_cpus} CPU(s) for parallel sweep')
