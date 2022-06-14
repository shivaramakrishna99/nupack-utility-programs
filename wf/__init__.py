"""
Carry out all NUPACK's Utilities functions
"""

from enum import Enum
from pathlib import Path
from latch import small_task, workflow
from latch.types import LatchFile
from natsort import as_ascii

from nupack import *  # Import NUPACK
import matplotlib.pyplot as plt

class Material(Enum):
    dna = "DNA"
    rna = "RNA"
    rna95 = "RNA95"

class Ensemble(Enum):
    stacking = "stacking"
    nostacking = "nostacking"


# Define Model() object and parameters   

@small_task
def utilities(
    strand1: str = "CCC",
    strand2: str = "GGG",
    structure: str = "(((+)))",
    partition_fn: bool = False,
    structure_energy_fn: bool = False,
    structure_probability_fn: bool = False,
    pairs_fn: bool = False,
    mfe_proxy_fn: bool = False,
    ensemble_size_fn: bool = False,
    ensemble_defect_fn: bool = False,
    material: str = Material.rna,
    ensemble: str = Ensemble.stacking,
    temperature: float = 37.0,
    sodium: float = 1.0,
    magnesium: float = 0.0,
    out: str = "nupack-utils-results",
) -> LatchFile:

    nt_model = Model(material=material, ensemble=ensemble, celsius=temperature, sodium=sodium, magnesium=magnesium)

    # Sequence-based Utility Functions
    noinput = f"NOT CALCULATED"

    if partition_fn == True:
        partition_fn_val = pfunc(strands=[strand1, strand2], model=nt_model)
        partition_fn_res = f"Partition function:\n{partition_fn_val}"
    else:
        partition_fn_res = f"Partition function: {noinput}"

    if pairs_fn == True:
        pairs_fn_val = pairs(strands=[strand1, strand2], model=nt_model)
        pairs_fn_res = f"Equilibrium base pairing probabilities:\n\n{pairs_fn_val}"
    else:
        pairs_fn_res = f"Equilibrium base pairing probabilities: {noinput}"

    if mfe_proxy_fn == True:
        mfe_structures = mfe(strands=[strand1, strand2], model=nt_model)
        mfe_proxy_fn_res = f"""
MFE Proxy Structure:
    Free energy of MFE proxy structure: {mfe_structures[0].energy}
    MFE proxy structure in dot-parens-plus notation: {mfe_structures[0].structure}
    MFE proxy structure as structure matrix:\n\n{mfe_structures[0].structure.matrix()}
"""
    else:
        mfe_proxy_fn_res = f"MFE Proxy Structure Results: {noinput}"
    
    if ensemble_size_fn == True:
        ensemble_size_fn_val = ensemble_size(strands=[strand1, strand2], model=nt_model)
        if ensemble == "stacking":
            ensemble_size_fn_res = f"Number of structures: {ensemble_size_fn_val}"
        elif ensemble == "nostacking":
            ensemble_size_fn_res = f"Number of stacking states: {ensemble_size_fn_val}"
    else:
        ensemble_size_fn_res = f"Ensemble size results: {noinput}"

    # Structure-based Utility Functions

    if structure_energy_fn == True:
        structure_energy_fn_val = structure_energy(strands=[strand1, strand2], structure=structure, model=nt_model)
        structure_energy_fn_res = f"Structure free energy: {structure_energy_fn_val}"
    else:
        structure_energy_fn_res = f"Structure free energy: {noinput}"
    
    if structure_probability_fn == True:
        structure_probability_fn_val = structure_probability(strands=[strand1, strand2], structure=structure, model=nt_model)
        structure_probability_fn_res = f"Equilibrium structure probability: {structure_probability_fn_val}"
    else:
        structure_probability_fn_res = f"Equilibrium structure probability: {noinput}"
    
    if ensemble_defect_fn == True:
        ensemble_defect_fn_val = defect(strands=[strand1, strand2], structure=structure, model=nt_model)
        ensemble_defect_fn_res = f"Complex ensemble defect: {ensemble_defect_fn_val}"
    else:
        ensemble_defect_fn_res = f"Complex ensemble defect: {noinput}"



    outFile = f"/{out}.txt"

    content = f"""
----------SEQUENCE ONLY UTILITY RESULTS----------
    
{partition_fn_res}

{pairs_fn_res}

{mfe_proxy_fn_res}

{ensemble_size_fn_res}

----------STRUCTURE UTILITY RESULTS----------

{structure_energy_fn_res}

{structure_probability_fn_res}

{ensemble_defect_fn_res}

----------INPUT SPECIFICATIONS----------

Strand 1: {strand1}
Strand 2: {strand2}
Structure: {structure}

Nucleic acid parameter set: {material}
Ensemble Type: {ensemble}
Temperature: {temperature} °C
Na+: {sodium} M
Mg++: {magnesium} M
"""
    with open(outFile, "w") as f:
        f.write(content)
    
    return LatchFile(outFile, f"latch:///{outFile}")

@workflow
def utilitiesNUPACK(
    strand1: str = "CCC",
    strand2: str = "GGG",
    structure: str = "(((+)))",
    partition_fn: bool = False,
    structure_energy_fn: bool = False,
    structure_probability_fn: bool = False,
    pairs_fn: bool = False,
    mfe_proxy_fn: bool = False,
    ensemble_size_fn: bool = False,
    ensemble_defect_fn: bool = False,
    material: Material = Material.rna,
    ensemble: Ensemble = Ensemble.stacking,
    temperature: float = 37.0,
    sodium: float = 1.0,
    magnesium: float = 0.0,
    out: str = "nupack-utils-results",
) -> LatchFile:
    """Carry out NUPACK Utilities functions on a complex of sequences

    # NUPACK - Utility Algorithms and Calculations
    ---

    ## **About**
    ---
    [NUPACK](https://docs.nupack.org/#about) is a growing software suite for the analysis and design of nucleic acid structures, devices, and systems serving the needs of researchers in the fields of nucleic acid nanotechnology, molecular programming, synthetic biology, and across the life sciences more broadly.

    ## **Citations**
    ---
    ### NUPACK Analysis Algorithms
    **Complex analysis and test tube analysis**
	- M.E. Fornace, N.J. Porubsky, and N.A. Pierce (2020). A unified dynamic programming framework for the analysis of interacting nucleic acid strands: enhanced models, scalability, and speed.  [ACS Synth Biol](https://pubs.acs.org/doi/abs/10.1021/acssynbio.9b00523) , 9:2665-2678, 2020. ( [pdf](http://www.nupack.org/downloads/serve_public_file/fornace20.pdf?type=pdf) ,  [supp info](http://www.nupack.org/downloads/serve_public_file/fornace20_supp.pdf?type=pdf) )
	
    - R. M. Dirks, J. S. Bois, J. M. Schaeffer, E. Winfree, and N. A. Pierce. Thermodynamic analysis of interacting nucleic acid strands.  [SIAM Rev](http://epubs.siam.org/doi/abs/10.1137/060651100) , 49:65-88, 2007. ( [pdf](http://www.nupack.org/downloads/serve_public_file/sirev07.pdf?type=pdf) )
    **Pseudoknot analysis**
	
    - R. M. Dirks and N. A. Pierce. An algorithm for computing nucleic acid base-pairing probabilities including pseudoknots.  [J Comput Chem](http://onlinelibrary.wiley.com/doi/10.1002/jcc.10296/abstract) , 25:1295-1304, 2004. ( [pdf](http://www.nupack.org/downloads/serve_public_file/jcc04.pdf?type=pdf) )
	
    - R. M. Dirks and N. A. Pierce. A partition function algorithm for nucleic acid secondary structure including pseudoknots.  [J Comput Chem](http://onlinelibrary.wiley.com/doi/10.1002/jcc.20057/abstract) , 24:1664-1677, 2003. ( [pdf](http://www.nupack.org/downloads/serve_public_file/jcc03.pdf?type=pdf) ,  [supp info](http://www.nupack.org/downloads/serve_public_file/jcc03_supp.pdf?type=pdf) )

    1. **Workflow Repository** - https://github.com/shivaramakrishna99/nupack-loop-stack/
    
    2. **Acknowledgements** - https://docs.nupack.org/#acknowledgments

    3. **License** - https://docs.nupack.org/#license
    
    4. **Citation** - https://docs.nupack.org/#citation
    ---

    __metadata__:
        display_name: NUPACK - Utility Functions
        
        author:
            name: The NUPACK Team
            email: support@nupack.org
            github: https://github.com/beliveau-lab/NUPACK
        
        repository: https://github.com/beliveau-lab/NUPACK
        
        license:
            id: BSD-3-Clause

    Args:
        material:
            __metadata__:
                display_name: "Nucleic Acid Type"
                _tmp:
                    section_title: "Sequence and Structure Details"
                appearance:
                    comment: "Choose between DNA and RNA free energy parameter sets. Default is 'rna', based on Matthews et al., 1999"
        strand1:
            __metadata__:
                display_name: "Strand 1 (as nucleotides)"
                appearance:
                    comment: "Enter the nucleotide sequence of the first strand"
        strand2:
            __metadata__:
                display_name: "Strand 2 (as nucleotides)"
                appearance:
                    comment: "Enter the nucleotide sequence of the second strand"
        structure:
            __metadata__:
                display_name: "Structure of Complex (in regular/extended dot-bracket notation)"
                appearance:
                    comment: "Enter the dot bracket notation of the intended structure"
        partition_fn:
            __metadata__:
                display_name: "Partition Function"
                _tmp:
                    section_title: "Choose NUPACK Utility Functions - sequence input based functions"
                appearance:
                    comment: "Compute the partition function https://www.google.com"
        pairs_fn:
            __metadata__:
                display_name: "Equilibrium Base Pairing Probabilities"
                appearance:
                    comment: ""
        mfe_proxy_fn:
            __metadata__:
                display_name: "MFE Proxy Structure(s)"
                appearance:
                    comment: ""
        ensemble_size_fn:
            __metadata__:
                display_name: "Complex Ensemble Size"
                appearance:
                    comment: ""
        structure_energy_fn:
            __metadata__:
                display_name: "Structure Free Energy"
                _tmp:
                    section_title: "Choose NUPACK Utility Functions - sequence and structural input based functions"
                appearance:
                    comment: "Calculate the structure energy"
        structure_probability_fn:
            __metadata__:
                display_name: "Equilibrium Structure Probability"
                appearance:
                    comment: ""
        ensemble_defect_fn:
            __metadata__:
                display_name: "Complex Ensemble Defect"
                appearance:
                    comment: "ED"
        ensemble:
            __metadata__:
                display_name: "Ensemble Type"
                _tmp:
                    section_title: Additional Model Specification
                    hidden: true
                appearance:
                    comment: "Choose between stacking and non stacking ensemble states. Default is set to 'stacking'."
        temperature:
            __metadata__:
                display_name: "Temperature (in degree Celsius)"
                _tmp:
                    hidden: true
                appearance:
                    comment: "Temperature of system. Default is 37 °C"
        sodium:
            __metadata__:
                display_name: "Na+ concentration (in M)"
                _tmp:
                    hidden: true
                appearance:
                    comment: "The sum of the concentrations of (monovalent) sodium, potassium, and ammonium ions, is specified in units of molar. Default: 1.0, Range: [0.05,1.1]"
        magnesium:
            __metadata__:
                display_name: "Mg++ (in nM). Default is 0 nM"
                _tmp:
                    hidden: true
                appearance:
                    comment: "The concentration of (divalent) magnesium ions, is specified in units of molar. Default: 0.0, Range: [0.0,0.2]"
        out:
            __metadata__:
                display_name: "Output File Name"
                _tmp:
                    section_title: Output
                appearance:
                    comment: "Name your file containing results from chosen NUPACK Utility Programs"
    """
    return utilities(
    strand1=strand1,
    strand2=strand2,
    structure=structure,
    partition_fn=partition_fn,
    structure_energy_fn=structure_energy_fn,
    structure_probability_fn=structure_probability_fn,
    pairs_fn=pairs_fn,
    mfe_proxy_fn=mfe_proxy_fn,
    ensemble_size_fn=ensemble_size_fn,
    ensemble_defect_fn=ensemble_defect_fn,
    material=material,
    ensemble=ensemble,
    temperature=temperature,
    sodium=sodium,
    magnesium=magnesium,
    out=out
    )
