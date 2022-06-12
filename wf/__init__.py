"""
Carry out all NUPACK's Utilities functions
"""

from enum import Enum
from pathlib import Path
import subprocess

from latch import small_task, workflow
from latch.types import LatchFile, LatchDir
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
    strand1: str = "ATGC",
    strand2: str = "ATGC",
    structure: str = "....",
    pfunc: bool = False,
    structure_energy: bool = False,
    structure_probability: bool = False,
    pairs: bool = False,
    mfe: bool = False,
    ensemble_size: bool = False,
    defect: bool = False,
    material: str = Material.dna,
    ensemble: str = Ensemble.stacking,
    temperature: float = 37.0,
    sodium: float = 1.0,
    magnesium: float = 0.0,
    outputDir: str = "output"
) -> LatchDir:

    nt_model = Model(material=material, celsius=temperature, sodium=sodium, magnesium=magnesium)
    # partition_function = pfunc(strands=['CCC', 'GGG'], model=nt_model)

    content = f"""
        ----------OUTPUT----------
        {nt_model}{strand1}{strand2}{structure}
    """

    return LatchDir(f"/root/{outputDir}", f"latch:///{outputDir}")

@workflow
def utilitiesNUPACK(
    strand1: str = "ATGC",
    strand2: str = "ATGC",
    structure: str = "....",
    pfunc: bool = False,
    structure_energy: bool = False,
    structure_probability: bool = False,
    pairs: bool = False,
    mfe: bool = False,
    ensemble_size: bool = False,
    defect: bool = False,
    material: str = Material.dna,
    ensemble: str = Ensemble.stacking,
    temperature: float = 37.0,
    sodium: float = 1.0,
    magnesium: float = 0.0,
    outputDir: str = "nupack-utils-run"
) -> LatchDir:
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
                    section_title: Sequence and Structure Details
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

        pfunc:
            __metadata__:
                display_name: "Partition Function"
                _tmp:
                    section_title: Choose NUPACK Utility Functions
                appearance:
                    comment: "Compute the partition function https://www.google.com"

        structure_energy:
            __metadata__:
                display_name: "Structure Free Energy"
                appearance:
                    comment: "Calculate the structure energy"

        structure_probability:
            __metadata__:
                display_name: "Equilibrium Structure Probability"
                appearance:
                    comment: ""

        pairs:
            __metadata__:
                display_name: "Equilibrium Base Pairing Probabilities"
                appearance:
                    comment: ""

        mfe:
            __metadata__:
                display_name: "MFE Proxy Structure(s)"
                appearance:
                    comment: ""

        ensemble_size:
            __metadata__:
                display_name: "Complex Ensemble Size"
                appearance:
                    comment: ""

        defect:
            __metadata__:
                display_name: "Complex Ensemble Size"
                appearance:
                    comment: ""

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

        outputDir:
            __metadata__:
                display_name: "Output Name"
                _tmp:
                    section_title: Output Specification
                appearance:
                    comment: "Specify the name of your output directory."
    """
    return utilities(
    strand1=strand1,
    strand2=strand2,
    structure=structure,
    pfunc=pfunc,
    structure_energy=structure_energy,
    structure_probability=structure_probability,
    pairs=pairs,
    mfe=mfe,
    ensemble_size=ensemble_size,
    defect=defect,
    material=material,
    ensemble=ensemble,
    temperature=temperature,
    sodium=sodium,
    magnesium=magnesium,
    outputDir=outputDir,
    )
