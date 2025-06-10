import numpy as np

EPSILON = 1e-5
Q_PARAM = 2.0
M_GRID, N_GRID = 20, 56
RT_PARAM = 0.2
R_PARAM = 2.0
X_MAX, Y_MAX = 1.0, 14.0
DTAU = 0.01
DN_PARAM, DM_PARAM = 0.5, 0.5
TAU_MAX = 20.0
GRASHOF_NUM = 5.0
PRANDTL_NUM = 0.71
G_MOD = 5.0
SCHMIDT_NUM = 0.3

M_PLUS_1 = M_GRID + 1
N_PLUS_1 = N_GRID + 1
DELTA_X = X_MAX / M_GRID
DELTA_Y = Y_MAX / N_GRID
DELTA_Y_SQ = DELTA_Y ** 2
DELTA_Y_SQ_PR = DELTA_Y_SQ * PRANDTL_NUM
DELTA_Y_SQ_SC = DELTA_Y_SQ * SCHMIDT_NUM

def solve_tridiagonal_system(first_idx, last_idx, a, b, c, d, x_out):
    n = len(b)
    beta = np.zeros(n)
    gamma = np.zeros(n)

    beta[first_idx] = b[first_idx]
    if beta[first_idx] == 0:
        raise ValueError("Zero pivot encountered at first_idx. The matrix might be singular or ill-conditioned.")

    gamma[first_idx] = d[first_idx] / beta[first_idx]

    for i in range(first_idx + 1, last_idx + 1):
        if beta[i - 1] == 0:
            raise ValueError(f"Zero pivot encountered at index {i - 1}. The matrix might be singular or ill-conditioned.")
        beta[i] = b[i] - a[i] * c[i - 1] / beta[i - 1]
        gamma[i] = (d[i] - a[i] * gamma[i - 1]) / beta[i]

    x_out[last_idx] = gamma[last_idx]
    for k in range(1, last_idx - first_idx + 1):
        i = last_idx - k
        x_out[i] = gamma[i] - c[i] * x_out[i + 1] / beta[i]

def initialize_simulation_arrays(m_plus_1, n_plus_1):
    size_2d = (m_plus_1, n_plus_1)
    size_1d_n = (n_plus_1,)
    size_1d_large = (1081,)

    return {
        "U": np.zeros(size_2d), "U_NEW": np.zeros(size_2d),
        "V": np.zeros(size_2d), "V_NEW": np.zeros(size_2d),
        "T": np.zeros(size_2d), "T_NEW": np.zeros(size_2d),
        "W": np.zeros(size_2d), "W_NEW": np.zeros(size_2d),
        "C_OLD": np.zeros(size_2d), "C_NEW": np.zeros(size_2d),

        "DTY_ARRAY": np.zeros(size_1d_large),
        "DTY_BT": np.zeros(size_1d_n),
        "DCY_ARRAY": np.zeros(size_1d_n),
        "DCY_BT": np.zeros(size_1d_n),
        "X_SH": np.zeros(size_1d_n),
        "X_NU": np.zeros(size_1d_n),
        "SKF": np.zeros(size_1d_n),
        "SKI": np.zeros(size_1d_n),

        "A_COEFF": np.zeros(size_1d_n),
        "B_COEFF": np.zeros(size_1d_n),
        "C_COEFF": np.zeros(size_1d_n),
        "D_COEFF": np.zeros(size_1d_n),
        "E_SOLUTION": np.zeros(size_1d_n)
    }

def run_fluid_simulation():
    print("\n--- FLUID DYNAMICS SIMULATION START ---")
    print(f"Grid Dimensions: M = {M_GRID}, N = {N_GRID}")
    print(f"Computational Domain: X_MAX = {X_MAX}, Y_MAX = {Y_MAX}")
    print(f"Grid Spacing: Delta X = {DELTA_X:.3f}, Delta Y = {DELTA_Y:.3f}")
    print(f"Time Step: DTAU = {DTAU}, Maximum Simulation Time: TAU_MAX = {TAU_MAX}")
    print(f"Key Parameters: Q = {Q_PARAM}, R = {R_PARAM}, Grashof = {GRASHOF_NUM}, "
          f"Prandtl = {PRANDTL_NUM}, Schmidt = {SCHMIDT_NUM}")
    print("-" * 40)

    data = initialize_simulation_arrays(M_PLUS_1, N_PLUS_1)

    U, U_NEW = data["U"], data["U_NEW"]
    V, V_NEW = data["V"], data["V_NEW"]
    T, T_NEW = data["T"], data["T_NEW"]
    W, W_NEW = data["W"], data["W_NEW"]
    C_OLD, C_NEW = data["C_OLD"], data["C_NEW"]
    A, B, C_diag, D, E = data["A_COEFF"], data["B_COEFF"], data["C_COEFF"], data["D_COEFF"], data["E_SOLUTION"]

    for i in range(M_PLUS_1):
        U[i, 0] = U_NEW[i, 0] = 1.0
        T[i, 0] = T_NEW[i, 0] = 1.0

    tau = 0.0
    time_step = 0
    has_converged = False

    print("\nStarting time integration loop...")
    while tau < TAU_MAX:
        tau += DTAU
        time_step += 1

        U[:], V[:], T[:], W[:] = U_NEW[:], V_NEW[:], T_NEW[:], W_NEW[:]
        C_OLD[:] = C_NEW[:]

        diff_U = np.max(np.abs(U_NEW - U))
        diff_T = np.max(np.abs(T_NEW - T))
        diff_W = np.max(np.abs(W_NEW - W))
        diff_C = np.max(np.abs(C_NEW - C_OLD))

        if diff_U < EPSILON and diff_T < EPSILON \
           and diff_W < EPSILON and diff_C < EPSILON:
            has_converged = True
            break

        if time_step % 100 == 0:
            print(f"Time Step: {time_step}, Current τ = {tau:.3f}, "
                  f"Max Diff U: {diff_U:.2e}, T: {diff_T:.2e}, "
                  f"W: {diff_W:.2e}, C: {diff_C:.2e}")

    print("\n--- SIMULATION COMPLETE ---")
    if has_converged:
        print(f"Simulation Converged in {time_step} steps at dimensionless time τ = {tau:.3f}.")
    else:
        print(f"Simulation did NOT converge within {TAU_MAX} maximum time. Final τ = {tau:.3f}.")
        print("Consider increasing TAU_MAX, adjusting DTAU, or checking model stability.")

if __name__ == "__main__":
    run_fluid_simulation()
