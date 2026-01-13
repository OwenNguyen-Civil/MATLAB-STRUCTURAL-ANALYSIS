# MATLAB-STRUCTURAL-ANALYSIS
MATLAB Structural Analysis Engine: A robust solver for arbitrary polygonal cross-sections. Features include geometric property calculation, Normal Stress (Navier), Pixel-based Shear Stress (Zhuravskii), and Von Mises failure criteria analysis. Capable of handling complex, asymmetric, and hollow shapes.
# 🏗️ PolyStress: General Structural Analysis Engine

![Banner](Images/vietnam.png)
> *Figure 1: Stress distribution analysis on a complex "Godzilla" polygon (Non-convex, 100+ vertices), demonstrating the engine's robustness against arbitrary geometries.*

---

## 📋 Overview
**PolyStress** is a MATLAB-based structural analysis engine designed to calculate geometric properties and stress distributions for **arbitrary polygonal cross-sections**.

Unlike traditional tools that are limited to standard shapes (I-beams, Rectangles), this engine utilizes **Green's Theorem** for exact geometric properties and a **Numerical Slicing Method (Pixel-based)** for shear flow analysis. This allows it to handle:
* **Complex Shapes:** Asymmetric, non-convex, or custom polygons (e.g., Map boundaries).
* **Hollow Sections:** Box girders, tubes (via cut-in/cut-out logic).
* **Mixed Loading:** Combined Axial Force ($N_z$), Biaxial Bending ($M_x, M_y$), and Shear Force ($V_y$).

## ✨ Key Features
* **📐 Geometry Engine:** - Calculates Area, Centroid ($C_x, C_y$), Moments of Inertia ($I_x, I_y$), and Product of Inertia ($I_{xy}$).
  - Automatically determines and rotates Principal Axes ($\alpha$).
* **🔴 Normal Stress Analysis:** - Computes stress based on the General Navier Formula: $\sigma = \frac{N}{A} + \frac{M_x I_y - M_y I_{xy}}{I_x I_y - I_{xy}^2}y + \dots$
  - Visualizes the rotation and translation of the **Neutral Axis (NA)**.
* **✂️ Shear Stress Analysis:** - Implements Zhuravskii's formula ($\tau = \frac{VQ}{Ib}$) using a discrete slicing algorithm to handle variable widths ($b_y$).
* **🛡️ Failure Criteria:** - Calculates **Von Mises Stress** ($\sigma_{vm} = \sqrt{\sigma^2 + 3\tau^2}$).
  - Determines **Safety Factor (SF)** based on material Yield Strength.
* **📊 Visualization:** - High-resolution Heatmaps for stress distribution.
  - Automatic detection of critical points.

---

## 📂 Project Structure
The repository is organized as follows:

```text
My-Structural-Analysis-Repo/
├── Core engine/           # 🧠 THE ENGINE (Source Code)
│   ├── calc_normal_stress.m  # Calculates Geom Props & Normal Stress
│   ├── calc_shear_stress.m   # Calculates Shear Stress (Pixel Method)
│   ├── calc_von_mises.m      # Checks Failure Criteria (Von Mises)
│   └── PolyStress_Engine.m   # (Optional) Main wrapper class
│
├── Legacy_v1/             # 🏛️ ARCHIVE (Old scripts)
│   ├── solid.m
│   ├── hollow.m
│   └── circular.m
│
├── Images/                # 📷 Gallery
│   ├── vietnam.png
│   ├── heart.png
│   └── ...
│
├── Main_Script.m          # 🚀 RUN THIS FILE (Demo)
└── README.md              # Documentation
