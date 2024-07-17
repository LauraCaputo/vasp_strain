# POSCAR Strain Application Script

This repository contains a Python script for applying homogeneous strain to VASP POSCAR files. The script prompts the user for strain percentages along the principal axes and optional shear strains, and then modifies the POSCAR file to reflect these changes.

## Features

- Applies homogeneous strain to the lattice vectors in a POSCAR file.
- Adjusts atomic positions to fit the new cell dimensions using fractional coordinates.
- Supports both direct (fractional) and Cartesian atomic positions.

## Requirements

- Python 3.x
- NumPy library

## Installation

Clone this repository to your local machine:

```bash
git clone https://github.com/LauraCaputo/vasp_strain.git
