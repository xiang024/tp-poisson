import numpy as np
import matplotlib.pyplot as plt
import os

# Directory of this script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Directory where results are stored
results_dir = os.path.join(script_dir, "..", "results_iter")

# Load residual histories (1D arrays)
r_alpha = np.loadtxt(os.path.join(results_dir, "RESVEC_richardson.dat"))
r_jac   = np.loadtxt(os.path.join(results_dir, "RESVEC_jacobi.dat"))
r_gs    = np.loadtxt(os.path.join(results_dir, "RESVEC_gs.dat"))

plt.figure(figsize=(6, 4))

plt.semilogy(r_alpha, label="Richardson (alpha optimal)")
plt.semilogy(r_jac, label="Jacobi")
plt.semilogy(r_gs, label="Gauss-Seidel")

plt.xlabel("Itération")
plt.ylabel(r"$\|r^{(k)}\|_2 / \|b\|_2$")
plt.title("Historique de convergence des méthodes itératives")
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.legend()

plt.tight_layout()
plt.savefig("iter_convergence.pdf")
plt.show()

