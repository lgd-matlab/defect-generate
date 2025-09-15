#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standalone Defect Formation and Binding Energy Calculator

A self-contained interactive script for calculating defect formation energies
and binding energies from VASP OUTCAR files. No external dependencies required
beyond Python standard library.

Features:
- Interactive command-line interface
- OUTCAR convergence checking
- Energy extraction from VASP calculations
- Formation energy calculation for vacancies, interstitials, substitutions
- Binding energy calculation for defect complexes
- Support for both neutral and charged defects
- Configuration save/load functionality

Usage:
    python defect_energy_calculator.py

Author: Generated for PyDefect workflow
Date: 2025-01-04
"""

import os
import re
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, asdict
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Regular expressions for energy extraction
RE_SIGMA0 = re.compile(r"energy\(sigma->0\)\s*=\s*([+\-]?[0-9]*\.?[0-9]+(?:[eE][+\-]?\d+)?)")
RE_CONVERGENCE = re.compile(r"reached required accuracy", re.IGNORECASE)

@dataclass
class DefectConfig:
    """Configuration for a single defect calculation"""
    label: str
    defect_type: str  # vacancy, interstitial, substitution
    outcar_path: str
    element_X: Optional[str] = None
    element_Y: Optional[str] = None
    charge: int = 0
    fermi_level: float = 0.0
    reference_level: float = 0.0
    potential_alignment: float = 0.0
    correction_energy: float = 0.0
    multiplicity: int = 1

@dataclass
class ElementalReference:
    """Elemental reference for chemical potential calculation"""
    element: str
    outcar_path: str
    poscar_path: str

@dataclass
class ComplexConfig:
    """Configuration for defect complex binding energy calculation"""
    label: str
    outcar_path: str
    defect_type: str
    element_X: Optional[str] = None
    element_Y: Optional[str] = None
    component_labels: List[str] = None
    charge: int = 0
    fermi_level: float = 0.0
    reference_level: float = 0.0
    potential_alignment: float = 0.0
    correction_energy: float = 0.0

@dataclass
class CalculationConfig:
    """Complete calculation configuration"""
    host_outcar: str
    elemental_refs: List[ElementalReference]
    defects: List[DefectConfig]
    complexes: List[ComplexConfig]
    calculate_formation: bool = True
    calculate_binding: bool = False

class DefectEnergyCalculator:
    """Main calculator class for defect energies"""

    def __init__(self):
        self.results = {}

    def read_file_safely(self, filepath: str) -> str:
        """Safely read file content with error handling"""
        try:
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                return f.read()
        except Exception as e:
            raise FileNotFoundError(f"Cannot read file {filepath}: {e}")

    def check_convergence(self, outcar_path: str) -> bool:
        """Check if OUTCAR calculation converged"""
        content = self.read_file_safely(outcar_path)
        return bool(RE_CONVERGENCE.search(content))

    def extract_energy(self, outcar_path: str) -> float:
        """Extract energy(sigma->0) from OUTCAR"""
        content = self.read_file_safely(outcar_path)
        matches = list(RE_SIGMA0.finditer(content))

        if not matches:
            # Fallback search
            lines = content.splitlines()
            for line in reversed(lines):
                if "energy(sigma->0)" in line:
                    match = RE_SIGMA0.search(line)
                    if match:
                        return float(match.group(1))
            raise ValueError(f"energy(sigma->0) not found in {outcar_path}")

        return float(matches[-1].group(1))

    def count_atoms_in_poscar(self, poscar_path: str) -> int:
        """Count total atoms in POSCAR file"""
        content = self.read_file_safely(poscar_path)
        lines = [line.strip() for line in content.splitlines() if line.strip()]

        # Find the line with atom counts (should be integers)
        for i, line in enumerate(lines):
            parts = line.split()
            if all(part.isdigit() for part in parts):
                return sum(int(part) for part in parts)

        raise ValueError(f"Cannot determine atom count from POSCAR: {poscar_path}")

    def calculate_chemical_potentials(self, elemental_refs: List[ElementalReference]) -> Dict[str, float]:
        """Calculate chemical potentials from elemental references"""
        mu_values = {}

        for ref in elemental_refs:
            logger.info(f"Processing elemental reference for {ref.element}")

            # Check convergence
            if not self.check_convergence(ref.outcar_path):
                raise RuntimeError(f"Elemental OUTCAR not converged: {ref.outcar_path}")

            # Extract energy and atom count
            total_energy = self.extract_energy(ref.outcar_path)
            atom_count = self.count_atoms_in_poscar(ref.poscar_path)

            # Calculate per-atom chemical potential
            mu_values[ref.element] = total_energy / atom_count
            logger.info(f"μ[{ref.element}] = {mu_values[ref.element]:.6f} eV/atom")

        return mu_values

    def calculate_formation_energy(self, defect: DefectConfig, host_energy: float,
                                 mu_values: Dict[str, float]) -> float:
        """Calculate formation energy for a defect"""
        # Check convergence
        if not self.check_convergence(defect.outcar_path):
            raise RuntimeError(f"Defect OUTCAR not converged: {defect.outcar_path}")

        # Extract defect energy
        defect_energy = self.extract_energy(defect.outcar_path)

        # Calculate reservoir term based on defect type
        reservoir_term = 0.0

        if defect.defect_type == 'vacancy':
            if not defect.element_X or defect.element_X not in mu_values:
                raise KeyError(f"Chemical potential for {defect.element_X} not available")
            # For vacancy: n_X = -multiplicity, so reservoir = -(-multiplicity) * mu_X
            reservoir_term = defect.multiplicity * mu_values[defect.element_X]

        elif defect.defect_type == 'interstitial':
            if not defect.element_X or defect.element_X not in mu_values:
                raise KeyError(f"Chemical potential for {defect.element_X} not available")
            # For interstitial: n_X = +multiplicity, so reservoir = -multiplicity * mu_X
            reservoir_term = -defect.multiplicity * mu_values[defect.element_X]

        elif defect.defect_type == 'substitution':
            if not defect.element_X or not defect.element_Y:
                raise ValueError(f"Substitution requires both X and Y elements")
            if defect.element_X not in mu_values or defect.element_Y not in mu_values:
                raise KeyError(f"Chemical potentials for {defect.element_X} and {defect.element_Y} required")
            # For substitution: n_X = +multiplicity, n_Y = -multiplicity
            reservoir_term = (-defect.multiplicity * mu_values[defect.element_X] +
                            defect.multiplicity * mu_values[defect.element_Y])

        # Calculate charge term
        charge_term = defect.charge * (defect.fermi_level + defect.reference_level +
                                     defect.potential_alignment)

        # Total formation energy
        formation_energy = (defect_energy - host_energy + reservoir_term +
                          charge_term + defect.correction_energy)

        logger.info(f"Formation energy for {defect.label}: {formation_energy:.6f} eV")
        return formation_energy

    def calculate_binding_energy(self, complex_config: ComplexConfig,
                               component_energies: Dict[str, float],
                               host_energy: float, mu_values: Dict[str, float]) -> float:
        """Calculate binding energy for a defect complex"""
        # Calculate formation energy of the complex
        complex_defect = DefectConfig(
            label=complex_config.label,
            defect_type=complex_config.defect_type,
            outcar_path=complex_config.outcar_path,
            element_X=complex_config.element_X,
            element_Y=complex_config.element_Y,
            charge=complex_config.charge,
            fermi_level=complex_config.fermi_level,
            reference_level=complex_config.reference_level,
            potential_alignment=complex_config.potential_alignment,
            correction_energy=complex_config.correction_energy
        )

        complex_formation_energy = self.calculate_formation_energy(
            complex_defect, host_energy, mu_values)

        # Sum formation energies of components
        component_sum = sum(component_energies[label] for label in complex_config.component_labels)

        # Binding energy = sum of components - complex formation energy
        binding_energy = component_sum - complex_formation_energy

        nature = "attractive" if binding_energy > 0 else ("repulsive" if binding_energy < 0 else "neutral")
        logger.info(f"Binding energy for {complex_config.label}: {binding_energy:.6f} eV ({nature})")

        return binding_energy

    def run_calculation(self, config: CalculationConfig) -> Dict[str, Any]:
        """Run the complete defect energy calculation"""
        results = {"formation_energies": {}, "binding_energies": {}}

        # Check host convergence and extract energy
        logger.info("Processing host calculation...")
        if not self.check_convergence(config.host_outcar):
            raise RuntimeError(f"Host OUTCAR not converged: {config.host_outcar}")

        host_energy = self.extract_energy(config.host_outcar)
        logger.info(f"Host energy: {host_energy:.6f} eV")

        # Calculate chemical potentials
        logger.info("Calculating chemical potentials...")
        mu_values = self.calculate_chemical_potentials(config.elemental_refs)

        # Calculate formation energies
        if config.calculate_formation:
            logger.info("Calculating formation energies...")
            for defect in config.defects:
                formation_energy = self.calculate_formation_energy(defect, host_energy, mu_values)
                results["formation_energies"][defect.label] = formation_energy

        # Calculate binding energies
        if config.calculate_binding:
            logger.info("Calculating binding energies...")
            for complex_config in config.complexes:
                binding_energy = self.calculate_binding_energy(
                    complex_config, results["formation_energies"], host_energy, mu_values)
                results["binding_energies"][complex_config.label] = binding_energy

        return results

def validate_file_path(filepath: str, file_type: str = "file") -> bool:
    """Validate that a file path exists and is readable"""
    if not filepath or not os.path.exists(filepath):
        print(f"Error: {file_type} '{filepath}' does not exist.")
        return False
    if not os.path.isfile(filepath):
        print(f"Error: '{filepath}' is not a file.")
        return False
    return True

def get_file_path(prompt: str, file_type: str = "file") -> str:
    """Get and validate a file path from user input"""
    while True:
        filepath = input(prompt).strip()
        if validate_file_path(filepath, file_type):
            return filepath
        print("Please enter a valid file path.")

def get_yes_no(prompt: str, default: bool = False) -> bool:
    """Get yes/no input from user"""
    default_str = "Y/n" if default else "y/N"
    while True:
        response = input(f"{prompt} ({default_str}): ").strip().lower()
        if not response:
            return default
        if response in ['y', 'yes']:
            return True
        elif response in ['n', 'no']:
            return False
        print("Please enter 'y' or 'n'.")

def get_integer(prompt: str, default: int = 0, min_val: int = None, max_val: int = None) -> int:
    """Get integer input from user with validation"""
    while True:
        try:
            response = input(f"{prompt} (default: {default}): ").strip()
            if not response:
                return default
            value = int(response)
            if min_val is not None and value < min_val:
                print(f"Value must be >= {min_val}")
                continue
            if max_val is not None and value > max_val:
                print(f"Value must be <= {max_val}")
                continue
            return value
        except ValueError:
            print("Please enter a valid integer.")

def get_float(prompt: str, default: float = 0.0) -> float:
    """Get float input from user with validation"""
    while True:
        try:
            response = input(f"{prompt} (default: {default}): ").strip()
            if not response:
                return default
            return float(response)
        except ValueError:
            print("Please enter a valid number.")



def prompt_elemental_references() -> List[ElementalReference]:
    refs: List[ElementalReference] = []
    print("\nEnter elemental references for chemical potentials (μ_i = E_tot^elem / N):")
    while True:
        add = get_yes_no("Add an elemental reference?", default=(len(refs) == 0))
        if not add:
            if not refs:
                print("At least one elemental reference is required.")
                continue
            break
        element = input("  Element symbol (e.g., Zr): ").strip()
        outcar = get_file_path("  OUTCAR path: ")
        poscar = get_file_path("  POSCAR path: ")
        refs.append(ElementalReference(element=element, outcar_path=outcar, poscar_path=poscar))
    return refs


def prompt_defects() -> List[DefectConfig]:
    defects: List[DefectConfig] = []
    print("\nEnter defect entries (vacancy, interstitial, substitution):")
    while True:
        add = get_yes_no("Add a defect?", default=(len(defects) == 0))
        if not add:
            if not defects:
                print("At least one defect is required.")
                continue
            break
        label = input("  Defect label (e.g., V_Zr, Pu_i, Pu_on_Zr): ").strip()
        dtype = input("  Defect type [vacancy/interstitial/substitution]: ").strip().lower()
        if dtype not in ("vacancy", "interstitial", "substitution"):
            print("  Invalid type. Try again.")
            continue
        outcar = get_file_path("  Defect OUTCAR path: ")
        X = None
        Y = None
        if dtype in ("vacancy", "interstitial"):
            X = input("  Element X (e.g., Zr or Pu): ").strip()
        elif dtype == "substitution":
            X = input("  Element X (dopant/substituent): ").strip()
            Y = input("  Element Y (host site): ").strip()
        charge = get_integer("  Charge state q", default=0)
        ef = get_float("  Fermi level E_F (eV)", default=0.0)
        eref = get_float("  Reference level E_ref (eV, e.g., VBM)", default=0.0)
        dv = get_float("  Potential alignment ΔV (eV)", default=0.0)
        ecorr = get_float("  Correction energy E_corr (eV)", default=0.0)
        mult = get_integer("  Multiplicity (count)", default=1, min_val=1)
        defects.append(DefectConfig(label=label, defect_type=dtype, outcar_path=outcar,
                                    element_X=X, element_Y=Y, charge=charge, fermi_level=ef,
                                    reference_level=eref, potential_alignment=dv,
                                    correction_energy=ecorr, multiplicity=mult))
    return defects


def prompt_complexes(defects: List[DefectConfig]) -> List[ComplexConfig]:
    complexes: List[ComplexConfig] = []
    print("\nEnter defect complexes for binding energy calculation (optional):")
    while True:
        add = get_yes_no("Add a defect complex?", default=False)
        if not add:
            break
        label = input("  Complex label (e.g., Pu_i+V_Zr): ").strip()
        outcar = get_file_path("  Complex OUTCAR path: ")
        dtype = input("  Complex type [vacancy/interstitial/substitution]: ").strip().lower()
        if dtype not in ("vacancy", "interstitial", "substitution"):
            print("  Invalid type. Try again.")
            continue
        X = None
        Y = None
        if dtype in ("vacancy", "interstitial"):
            X = input("  Element X (e.g., Zr or Pu): ").strip()
        elif dtype == "substitution":
            X = input("  Element X (dopant/substituent): ").strip()
            Y = input("  Element Y (host site): ").strip()
        comps = input("  Component defect labels (comma-separated, must match added defects): ").strip()
        comp_labels = [c.strip() for c in comps.split(',') if c.strip()]
        for c in comp_labels:
            if c not in [d.label for d in defects]:
                print(f"  Unknown component label '{c}'. Please add the defect first.")
                comp_labels = []
                break
        if not comp_labels:
            continue
        charge = get_integer("  Charge state q", default=0)
        ef = get_float("  Fermi level E_F (eV)", default=0.0)
        eref = get_float("  Reference level E_ref (eV)", default=0.0)
        dv = get_float("  Potential alignment ΔV (eV)", default=0.0)
        ecorr = get_float("  Correction energy E_corr (eV)", default=0.0)
        complexes.append(ComplexConfig(label=label, outcar_path=outcar, defect_type=dtype,
                                       element_X=X, element_Y=Y, component_labels=comp_labels,
                                       charge=charge, fermi_level=ef, reference_level=eref,
                                       potential_alignment=dv, correction_energy=ecorr))
    return complexes


def save_config(config: CalculationConfig, path: str) -> None:
    data = {
        "host_outcar": config.host_outcar,
        "elemental_refs": [asdict(r) for r in config.elemental_refs],
        "defects": [asdict(d) for d in config.defects],
        "complexes": [asdict(c) for c in config.complexes],
        "calculate_formation": config.calculate_formation,
        "calculate_binding": config.calculate_binding,
    }
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print(f"Configuration saved to {path}")


def load_config(path: str) -> CalculationConfig:
    with open(path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    return CalculationConfig(
        host_outcar=data["host_outcar"],
        elemental_refs=[ElementalReference(**r) for r in data["elemental_refs"]],
        defects=[DefectConfig(**d) for d in data["defects"]],
        complexes=[ComplexConfig(**c) for c in data.get("complexes", [])],
        calculate_formation=data.get("calculate_formation", True),
        calculate_binding=data.get("calculate_binding", False),
    )


def interactive_main():
    print("\n=== Defect Formation & Binding Energy Calculator ===\n")
    calc = DefectEnergyCalculator()

    # Offer to load existing config
    if get_yes_no("Load an existing configuration file?", default=False):
        cfg_path = get_file_path("  Configuration JSON path: ")
        config = load_config(cfg_path)
    else:
        # Host OUTCAR
        host_outcar = get_file_path("Enter host OUTCAR path: ")
        # Elemental references
        refs = prompt_elemental_references()
        # Defects
        defects = prompt_defects()
        # Modes
        calc_formation = get_yes_no("Calculate formation energies?", default=True)
        calc_binding = get_yes_no("Calculate binding energies?", default=False)
        # Complexes (if needed)
        complexes: List[ComplexConfig] = []
        if calc_binding:
            complexes = prompt_complexes(defects)
        config = CalculationConfig(host_outcar=host_outcar,
                                   elemental_refs=refs,
                                   defects=defects,
                                   complexes=complexes,
                                   calculate_formation=calc_formation,
                                   calculate_binding=calc_binding)
        # Offer to save
        if get_yes_no("Save this configuration?", default=True):
            out_json = input("  Output config JSON path (default: defect_energy_config.json): ").strip() or "defect_energy_config.json"
            save_config(config, out_json)

    try:
        results = calc.run_calculation(config)
    except Exception as e:
        print(f"\nCalculation failed: {e}")
        sys.exit(1)

    # Summary
    print("\n=== Results Summary ===")
    if results.get("formation_energies"):
        print("Formation energies (eV):")
        for label, ef in results["formation_energies"].items():
            print(f"  {label}: {ef:.6f}")
    if results.get("binding_energies"):
        print("Binding energies (eV):")
        for label, eb in results["binding_energies"].items():
            nature = "attractive" if eb > 0 else ("repulsive" if eb < 0 else "neutral")
            print(f"  {label}: {eb:.6f} ({nature})")

if __name__ == "__main__":
    interactive_main()
