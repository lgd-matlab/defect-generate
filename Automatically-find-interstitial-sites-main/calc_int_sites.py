import numpy as np
import argparse
import sys
from pymatgen.analysis.defects.generators import  VoronoiInterstitialGenerator, ChargeInterstitialGenerator
from pymatgen.io.vasp import Poscar, Chgcar
from pymatgen.core import Element
from pymatgen.io.cif import CifWriter


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Find interstitial positions/sites in crystal structure automatically",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -e Li -n 3 -i POSCAR
  %(prog)s --element N --number 5 --input structure.POSCAR --which-site 1
  %(prog)s -e H -n 10 --charge-density CHGCAR --clustering-tol 0.8
        """
    )
    
    # Required arguments
    parser.add_argument('-e', '--element', 
                       required=True,
                       help='Interstitial element to place (e.g., N, Li, H)')
    
    parser.add_argument('-n', '--number', 
                       type=int, 
                       required=True,
                       help='Number of interstitials to insert')
    
    # Input file arguments
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-i', '--input', 
                            help='Input POSCAR file (default: POSCAR)')
    input_group.add_argument('--charge-density', 
                            help='Input CHGCAR file for charge density method')
    
    # Optional arguments
    parser.add_argument('--which-site', 
                       type=int, 
                       default=0,
                       help='Select specific interstitial site type (0=all, 1=first type, 2=second type, etc.) (default: 0)')
    
    parser.add_argument('--clustering-tol', 
                       type=float, 
                       default=0.75,
                       help='Clustering tolerance for site identification (default: 0.75)')
    
    parser.add_argument('--min-dist', 
                       type=float, 
                       default=0.5,
                       help='Minimum distance between interstitial sites (default: 0.5)')
    
    parser.add_argument('--output-prefix', 
                       default='modified_structure',
                       help='Output file prefix (default: modified_structure)')
    
    parser.add_argument('--mode-target', 
                       type=float, 
                       default=0.5,
                       help='Target value for moderate spacing mode (default: 0.5)')
    
    parser.add_argument('--verbose', 
                       action='store_true',
                       help='Enable verbose output')
    
    return parser.parse_args()


def main():
    """Main function to execute interstitial site calculation."""
    args = parse_arguments()
    
    # Set parameters from CLI arguments
    interstitial_element_to_place = args.element
    number_of_interstitials_to_insert = args.number
    which_interstitial_to_use = args.which_site
    
    if args.verbose:
        print(f"Parameters:")
        print(f"  Element: {interstitial_element_to_place}")
        print(f"  Number of interstitials: {number_of_interstitials_to_insert}")
        print(f"  Site type selection: {which_interstitial_to_use}")
        print(f"  Clustering tolerance: {args.clustering_tol}")
        print(f"  Minimum distance: {args.min_dist}")
        print()
    
    # Load structure based on input method
    if args.charge_density:
        if args.verbose:
            print(f"Loading structure from charge density file: {args.charge_density}")
        structure = Chgcar.from_file(args.charge_density)
        generator = ChargeInterstitialGenerator(
            clustering_tol=args.clustering_tol,
            min_dist=args.min_dist
        )
    else:
        input_file = args.input if args.input else 'POSCAR'
        if args.verbose:
            print(f"Loading structure from POSCAR file: {input_file}")
        structure = Poscar.from_file(input_file).structure
        generator = VoronoiInterstitialGenerator(
            clustering_tol=args.clustering_tol,
            min_dist=args.min_dist
        )
    
    # Execute the interstitial site calculation
    run_calculation(structure, generator, interstitial_element_to_place, 
                   number_of_interstitials_to_insert, which_interstitial_to_use,
                   args.output_prefix, args.mode_target, args.verbose)


def run_calculation(structure, generator, interstitial_element_to_place, 
                   number_of_interstitials_to_insert, which_interstitial_to_use,
                   output_prefix, mode_target, verbose):
    """Run the main interstitial site calculation logic."""
    
    """
    FIRST PART: Obtain all interstitial sites with fractional coordinates.
    """
    frac_coords = []
    unique_int = []
    unique_mult = []
    idx = 0
    frac_coords_dict = {}
    for interstitial in generator.generate(structure, "H"): #Element 'H' is here only to find the available sites in order to prevent error with oxidation states for some elements like noble gases
        frac_coords_dict[idx]=[]
        if verbose:
            print(f"\nUnique interstitial site at: {interstitial.site.frac_coords}")
            print(f"It has multiplicity of: {interstitial.multiplicity}")
            print(f"--------------------------------------------------\n\n")
            print(f"The following are all the equivalent positions (lattice coordinates) [fractional coordinates]\n:")
        unique_int.append(interstitial.site.frac_coords)
        unique_mult.append(interstitial.multiplicity)
        for site in interstitial.equivalent_sites:
            if verbose:
                print(f"Fractional: {site.frac_coords}")
            frac_coords.append(site.frac_coords)
            frac_coords_dict[idx].append(site.frac_coords)
        idx = idx + 1
    
    if verbose:
        print(f"\nThere are total of {len(unique_int)} unique interstitial sites at\n:{unique_int}, having multiplicity of {unique_mult})\n.If you wish to place interstitials only at certain type of interstitial positions, change the 'which_interstitial_to_use")

    if which_interstitial_to_use == 0:
            frac_coords_use = frac_coords
    else:
            frac_coords_use = frac_coords_dict[which_interstitial_to_use-1]

    # Calculate different spacing modes
    selected_points_farthest = select_spaced_points(frac_coords_use, n_points=number_of_interstitials_to_insert, mode='farthest')
    selected_points_nearest = select_spaced_points(frac_coords_use, n_points=number_of_interstitials_to_insert, mode='nearest')
    selected_points_moderate = select_spaced_points(frac_coords_use, n_points=number_of_interstitials_to_insert, mode='moderate', target_value=mode_target)

    # Save structures with different spacing modes
    save_structure_with_interstitials(structure, selected_points_farthest, interstitial_element_to_place, f"{output_prefix}_farthest", "FARTHEST", verbose)
    save_structure_with_interstitials(structure, selected_points_nearest, interstitial_element_to_place, f"{output_prefix}_nearest", "NEAREST", verbose)
    save_structure_with_interstitials(structure, selected_points_moderate, interstitial_element_to_place, f"{output_prefix}_moderate", "MODERATE", verbose)


def save_structure_with_interstitials(original_structure, selected_points, element, file_prefix, mode_name, verbose):
    """Save structure with interstitials to CIF and POSCAR files."""
    # Create a copy of the original structure
    structure_to_save = original_structure.copy()
    
    if verbose:
        print(f"Saving structure with {mode_name} distances between interstitials into POSCAR and CIF files...")
    
    for point in selected_points:
        structure_to_save.append(
            species=Element(element),
            coords=point,
            coords_are_cartesian=False  # Specify that the coordinates are fractional
         )

    cif_writer = CifWriter(structure_to_save)
    cif_writer.write_file(f"{file_prefix}.cif")  # Output CIF file
    poscar = Poscar(structure_to_save)
    poscar.write_file(f"{file_prefix}.POSCAR")  # Output POSCAR file


if __name__ == "__main__":
    main()

# Legacy support: If no CLI arguments provided, use original hardcoded values
if len(sys.argv) == 1:
    # USER INPUTS (Legacy mode)
    #-------------------------------------------
    interstitial_element_to_place = "N"
    number_of_interstitials_to_insert = 5
    which_interstitial_to_use = 0 # The value '0' will consider all found available interstitial positions for calculating
    #the distances. If you want to place interstitials on only the certain type of interstitial site, e.g., this method will
    #identify two different types of interstitial positions and you want to insert interstitial only on the first type,
    #then change the value for this variable to the number 1. If you want to place only on the second interstitial site type,
    #change it to value 2.
    #-------------------------------------------
    
    structure = Poscar.from_file("POSCAR").structure
    generator = VoronoiInterstitialGenerator(
        clustering_tol=0.75,
        min_dist=0.5,
    )

    """
    If you want to find the interstitial sites based on the charge density method instead of Voronoi method,
    uncomment the following two variables.This approach requires the charge density CHGCAR file from VASP as input
    (instead of the structural POSCAR file).
    Ensure you generate the CHGCAR file first by performing a VASP calculation on your initial structure

    generator = ChargeInterstitialGenerator( clustering_tol=0.75,
        min_dist=0.5)
    structure =  Chgcar.from_file("CHGCAR")
    """

    # Run the legacy calculation
    run_calculation(structure, generator, interstitial_element_to_place, 
                   number_of_interstitials_to_insert, which_interstitial_to_use,
                   'modified_structure', 0.5, False)
else:
    # CLI mode - main() already executed and exited
    pass


def wrap_coordinates(frac_coords):
    """
    Wrap fractional coordinates into the range [0, 1).
    """
    frac_coords = np.array(frac_coords)  # Ensure input is a NumPy array
    return frac_coords % 1


def compute_periodic_distance_matrix(frac_coords):
    """
    Compute a periodic distance matrix for a set of fractional coordinates.

    Args:
        frac_coords (ndarray): Fractional coordinates of shape (N, 3).

    Returns:
        ndarray: Distance matrix of shape (N, N).
    """
    n = len(frac_coords)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            delta = frac_coords[i] - frac_coords[j]
            delta = delta - np.round(delta)  # Apply periodic boundary conditions
            dist_matrix[i, j] = dist_matrix[j, i] = np.linalg.norm(delta)
    return dist_matrix


def select_spaced_points(frac_coords, n_points=5, mode="farthest", target_value=0.5):
    """
    Select n_points that are maximally spaced apart under periodic boundary conditions.

    Args:
        frac_coords (list of list of float): Fractional coordinates of points.
        n_points (int): Number of points to select.

    Returns:
        list of list of float: Selected fractional coordinates.
    """
    frac_coords = wrap_coordinates(frac_coords)  # Wrap fractional coordinates
    dist_matrix = compute_periodic_distance_matrix(frac_coords)
    # Greedy selection of points
    selected_indices = [0]  # Start with the first point
    for _ in range(1, n_points):
        remaining_indices = [i for i in range(len(frac_coords)) if i not in selected_indices]

        if mode == "farthest":
            # Select the point that is farthest from already selected points
            next_index = max(
                remaining_indices,
                key=lambda i: min(dist_matrix[i, j] for j in selected_indices)
            )
        elif mode == "nearest":
            # Select the point that is nearest to already selected points
            next_index = min(
                remaining_indices,
                key=lambda i: min(dist_matrix[i, j] for j in selected_indices)
            )
        elif mode == "moderate":
            # Select the point with average distance closest to the target_value
            next_index = min(
                remaining_indices,
                key=lambda i: abs(sum(dist_matrix[i, j] for j in selected_indices) / len(selected_indices) - target_value)
            )
        else:
            raise ValueError("Invalid mode. Choose from 'nearest', 'farthest', or 'moderate'.")

        selected_indices.append(next_index)

    return frac_coords[selected_indices].tolist()
