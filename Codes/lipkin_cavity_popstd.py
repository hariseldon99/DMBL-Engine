#!/usr/bin/env python3
import warnings

# Suppress all warnings before any other modules are imported
warnings.filterwarnings('ignore')

import argparse
import numpy as np
from mpi4py import MPI
from qutip import basis, destroy, jmat, mesolve, qeye, tensor
from scipy.special import jn_zeros
from tqdm.auto import tqdm

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


# ---------------------------------------------------------------------
# Self-contained helpers
# ---------------------------------------------------------------------
def build_ops(N, N_ph, J, omega_0, g):
    """Build operators for the LMG-cavity system."""
    N_ph = int(N_ph)  # Cast to standard Python int for QuTiP
    N = int(N)        # Ensure N is also an int

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
# Argument Parsing & Task Construction
# ---------------------------------------------------------------------
if rank == 0:
    parser = argparse.ArgumentParser(
        description='MPI LMG-Cavity System Parameter Sweep'
    )
    parser.add_argument('--N_spin', type=int, default=400, help='Spin size')
    parser.add_argument('--J', type=float, default=1.0, help='J parameter')
    parser.add_argument('--omega_0', type=float, default=1.0, help='Cavity frequency omega_0')
    parser.add_argument('--g', type=float, default=1.0, help='Coupling strength g')
    parser.add_argument('--T_total', type=float, default=3000.0, help='Total time')
    parser.add_argument('--dt', type=float, default=0.05, help='Time step')
    parser.add_argument('--Omega_sweep', type=float, default=0.7, help='Drive frequency Omega')
    parser.add_argument('--N_min', type=int, default=2, help='Min N_ph value')
    parser.add_argument('--N_max', type=int, default=400, help='Max N_ph value')
    parser.add_argument('--N_step', type=int, default=32, help='Step size for N_ph sweep')
    parser.add_argument(
        '--output',
        type=str,
        default='cavity_population_sweep_checkpoint.npz',
        help='Output filename',
    )

    args = parser.parse_args()

    # Derived parameter
    h0_sweep = jn_zeros(0, 5)[1] * args.Omega_sweep / 2

    N_values = np.arange(args.N_min, args.N_max, args.N_step)
    tasks = [
        (
            args.N_spin,
            N_ph,
            args.J,
            args.omega_0,
            args.g,
            h0_sweep,
            args.Omega_sweep,
            args.T_total,
            args.dt,
        )
        for N_ph in N_values
    ]

    # Chunk tasks for scatter among available MPI ranks
    task_chunks = np.array_split(tasks, size)
    task_chunks = [list(chunk) for chunk in task_chunks]
    
    # Store settings dictionary to pass for saving
    run_params = vars(args)
    run_params['h0_sweep'] = h0_sweep
else:
    task_chunks = None
    run_params = None

# Distribute chunks to ranks via Scatter
local_tasks = comm.scatter(task_chunks, root=0)

# Process local tasks assigned to this rank
local_results = []
if rank == 0:
    for task in tqdm(local_tasks, desc='Rank 0 progress'):
        local_results.append(cavity_population_std(task))
else:
    for task in local_tasks:
        local_results.append(cavity_population_std(task))

# Collect results from all processes back to rank 0
all_results = comm.gather(local_results, root=0)

# ---------------------------------------------------------------------
# Save output (Rank 0 only)
# ---------------------------------------------------------------------
if rank == 0:
    # Flatten the gathered list of lists
    std_results = [item for sublist in all_results for item in sublist]

    std_results.sort(key=lambda item: item[0])
    N_out = np.array([item[0] for item in std_results])
    std_n = np.array([item[1] for item in std_results])

    checkpoint_file = run_params['output']
    np.savez_compressed(
        checkpoint_file,
        N_out=N_out,
        std_n=std_n,
        N_spin=run_params['N_spin'],
        J=run_params['J'],
        omega_0=run_params['omega_0'],
        g=run_params['g'],
        Omega_sweep=run_params['Omega_sweep'],
        h0_sweep=run_params['h0_sweep'],
        T_total=run_params['T_total'],
        dt=run_params['dt'],
    )

    print(f'Saved sweep checkpoint to {checkpoint_file}')
    print(f'Completed execution using {size} MPI process(es).') 