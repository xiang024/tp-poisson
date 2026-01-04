import numpy as np
import matplotlib.pyplot as plt
import os

# Directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# CSV data file
data_path = os.path.join(script_dir, "direct_time.csv")

# Load CSV: skip header, delimiter is comma
data = np.loadtxt(data_path, delimiter=",", skiprows=1)

# Columns
# 0: method, 1: N, 2: time_ms
method = data[:, 0]
N = data[:, 1]
time_ms = data[:, 2]

# Split by method
data_trf = data[method == 0]
data_tri = data[method == 1]
data_sv  = data[method == 2]

plt.figure(figsize=(6, 4))

plt.loglog(data_trf[:, 1], data_trf[:, 2], 'o-', label="TRF (dgbtrf + dgbtrs)")
plt.loglog(data_tri[:, 1], data_tri[:, 2], 's-', label="TRI (LU tridiagonale)")
plt.loglog(data_sv[:, 1],  data_sv[:, 2],  '^-', label="SV (dgbsv)")

plt.xlabel("Taille du problème (N)")
plt.ylabel("Temps d'exécution (ms)")
plt.title("Performances des méthodes directes")
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.legend()

plt.tight_layout()
plt.savefig("direct_times.pdf")
plt.show()
