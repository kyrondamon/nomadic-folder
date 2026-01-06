
import numpy as np
import matplotlib.pyplot as plt
import time

# ==========================================
# CONFIGURATION
# ==========================================
SEQUENCE = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
KERNEL_WIDTH = 15     # Biological coherence length (residues)
BOND_LENGTH = 3.8     # C-alpha spacing (Angstroms)

# Kyte-Doolittle Hydrophobicity Scale
KD_SCALE = {
    'A': 1.8, 'R':-4.5, 'N':-3.5, 'D':-3.5, 'C': 2.5,
    'Q':-3.5, 'E':-3.5, 'G':-0.4, 'H':-3.2, 'I': 4.5,
    'L': 3.8, 'K':-3.9, 'M': 1.9, 'F': 2.8, 'P':-1.6,
    'S':-0.8, 'T':-0.7, 'W':-0.9, 'Y':-1.3, 'V': 4.2
}

# ==========================================
# SIGNAL PROCESSING
# ==========================================
def get_hydrophobicity_signal(seq):
    h = np.array([KD_SCALE.get(aa, 0.0) for aa in seq])
    return h - np.mean(h)

def ricker_wavelet(width):
    t = np.linspace(-width/2, width/2, width)
    sigma = width / 4.0
    norm_t = t / sigma
    return (1 - norm_t**2) * np.exp(-0.5 * norm_t**2)

def compute_nomadic_response(signal, width):
    kernel = ricker_wavelet(width)
    return np.convolve(signal, kernel, mode='same')

# ==========================================
# NOMADIC FOLDER ENGINE
# ==========================================
class NomadicFolder:
    def __init__(self, sequence):
        self.seq = sequence
        self.n = len(sequence)
        self.positions = self._initialize_chain()

        raw_h = get_hydrophobicity_signal(self.seq)
        self.response = compute_nomadic_response(raw_h, KERNEL_WIDTH)

        r_col = self.response.reshape(-1, 1)
        r_row = self.response.reshape(1, -1)
        self.affinity = r_col * r_row

    def _initialize_chain(self):
        pos = np.zeros((self.n, 3))
        for i in range(1, self.n):
            pos[i] = pos[i-1] + np.array([BOND_LENGTH, 0, 0])
        return pos

    def step(self, step_size=0.1):
        forces = np.zeros_like(self.positions)

        # A. TOPOLOGICAL ATTRACTION
        strong_contacts = np.argwhere(self.affinity > 20)
        for i, j in strong_contacts:
            if abs(i - j) < 4:
                continue
            vec = self.positions[j] - self.positions[i]
            dist = np.linalg.norm(vec)
            if dist > 5.0:
                forces[i] += (vec / dist) * self.affinity[i, j] * 0.003

        # B. CHAIN INTEGRITY
        for i in range(self.n - 1):
            vec = self.positions[i+1] - self.positions[i]
            dist = np.linalg.norm(vec)
            corr = (vec / dist) * (dist - BOND_LENGTH) * 0.5
            forces[i] += corr
            forces[i+1] -= corr

        # C. GLOBAL STERIC EXCLUSION
        for i in range(self.n):
            for j in range(i + 2, self.n):
                vec = self.positions[i] - self.positions[j]
                dist = np.linalg.norm(vec)
                if dist < 4.5 and dist > 0.01:
                    push = (vec / dist) * (4.5 - dist) * 2.0
                    forces[i] += push
                    forces[j] -= push

        self.positions += forces * step_size
        self.positions += np.random.randn(self.n, 3) * 0.05

    def get_rg(self):
        center = np.mean(self.positions, axis=0)
        return np.sqrt(np.mean(np.sum((self.positions - center)**2, axis=1)))

# ==========================================
# RUNNER
# ==========================================
if __name__ == "__main__":
    print(f"--- FOLDING UBIQUITIN ({len(SEQUENCE)} residues) ---")
    sim = NomadicFolder(SEQUENCE)
    history = []

    start = time.time()
    for i in range(600):
        sim.step()
        rg = sim.get_rg()
        history.append(rg)
        if i % 50 == 0:
            print(f"Step {i}: Rg = {rg:.2f} A")

    print(f"--- DONE in {time.time() - start:.2f}s ---")
    print(f"Final Radius of Gyration: {history[-1]:.2f} A")

    plt.plot(history, label="Trajectory", color="purple", linewidth=2)
    plt.axhline(12.0, color="green", linestyle="--", label="Native Target")
    plt.title("NomadicFolder: Ubiquitin Sprint")
    plt.xlabel("Step")
    plt.ylabel("Radius of Gyration (Å)")
    plt.legend()
    plt.show()
