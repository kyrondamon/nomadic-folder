# NomadicFolder

**Kernel-guided, sequence-specific protein compaction on consumer hardware**

NomadicFolder is a lightweight, experimental Python engine that demonstrates a fast,
sequence-specific approach to protein compaction. Instead of simulating atomic dynamics
via Molecular Dynamics (MD), it treats early-stage folding as a **signal-processing problem**:
global topological guidance is extracted directly from the amino-acid sequence using a
fixed convolutional kernel.

The goal of this project is **not** to replace MD or structure-prediction systems, but to
provide a rapid, physics-inspired **pre-folding primitive** that identifies global
constraints and compaction pathways at minimal computational cost.

---

## 🚀 Benchmark: The “Ubiquitin Sprint”

Running on consumer hardware, NomadicFolder demonstrates:

- **Target:** Ubiquitin (1UBQ, 76 residues)
- **Runtime:** ~20 seconds
- **Initial state:** Extended chain (~83 Å radius of gyration)
- **Final state:** Compact globular fold (~11.5 Å radius of gyration)
- **Behavior:** Sequence-specific collapse with non-intersecting geometry

A scrambled control sequence (same composition, randomized order) fails to compact,
stalling at large radius of gyration (>50 Å), demonstrating that **sequence order matters**.

---

## 🧠 How It Works (Mechanically)

1. **Sequence → Signal**  
   The amino-acid sequence is mapped to a one-dimensional hydrophobicity signal
   (Kyte–Doolittle scale, mean-centered).

2. **Fixed Kernel Convolution**  
   The signal is convolved with a fixed Ricker (Mexican-hat) wavelet.  
   This produces a response field highlighting sequence regions that resonate with
   a characteristic biological coherence length.

3. **Global Affinity Field**  
   The response is lifted to an interaction matrix via an outer product, creating a
   long-range, sequence-specific steering field.

4. **Geometric Relaxation**  
   A minimal geometric solver enforces:
   - Chain connectivity (bond-length constraints)
   - Global steric exclusion (hard-sphere repulsion)

   Under this guidance, the chain rapidly compacts into a globular topology.

---

## 📉 Verification: Ablation Test

To verify that the behavior is sequence-driven and not generic attraction:

- **Native ubiquitin sequence:** Converges to ~11.9 Å (folded)
- **Scrambled sequence:** Stalls at ~53.7 Å (unfolded)

All parameters are held fixed between runs.

---

## ❗ What This Is — and Is Not

**This is:**
- A fast, deterministic, physics-inspired compaction engine
- A global topology and constraint extractor
- A pre-filter or guide for downstream simulation or design

**This is not:**
- A full atomic folding engine
- A replacement for Molecular Dynamics
- A trained or data-driven model
- An AlphaFold-like predictor

---

## 📚 Origin of the Kernel

The convolutional kernel used here was originally developed in the context of the
**Nomadic Principle**, a global response framework proposed in galactic dynamics and
later explored in biological systems.

For the purposes of this repository, the kernel is treated as a **fixed, hand-specified
operator**. Users do not need to adopt or agree with the broader theoretical framework
to evaluate or use the code.

The associated paper is available on Zenodo:

> **Pre-dynamical ordering in protein folding via the Nomadic Principle**  
> Kyron Damon (2026)  
> https://doi.org/10.5281/zenodo.18160002

---

## 🛠 Requirements

- Python 3
- `numpy`
- `matplotlib`

No GPU, no training, no external datasets required.

---

## 📦 Status

This repository currently provides a reference implementation demonstrating
sequence-specific protein compaction. Future extensions may include:

- PDB export utilities
- Smooth soft-core steric potentials
- Backbone angular constraints
- Multi-domain support

---

Developed by an independent researcher.
