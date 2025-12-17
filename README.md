
# giFBA: Greedy Interaction Flux Balance Analysis

This is a COBRApy extension for modeling community interaction.

## Installation
```
pip install "git+https://github.com/m-3-lab/giFBA.git@main#subdirectory=package"
```

## Directory Description
```
📦giFBA
 ┣ 📂Examples
 ┃ ┣ 📂Model_Diagrams
 ┃ ┣ 📂Results
 ┃ ┃ ┣ 📂cFBA
 ┃ ┃ ┣ 📂giFBA_AGORA
 ┃ ┃ ┣ 📂giFBA_Simple
 ┃ ┃ ┣ 📂cFBA_Models
 ┃ ┃ ┗ micom_results.csv
 ┃ ┣ create_simple_fba.ipynb
 ┃ ┣ gifba_comparisons.ipynb
 ┃ ┗ gifba_utilization.ipynb
 ┣ 📂package
 ┃ ┣ 📂gifba
 ┃ ┃ ┣ 📂Simple_Models
 ┃ ┃ ┣ __init__.py
 ┃ ┃ ┣ config.py
 ┃ ┃ ┣ gifba_object.py
 ┃ ┃ ┣ summary.py
 ┃ ┃ ┗ utils.py
 ┃ ┣ README.md
 ┃ ┗ pyproject.toml
 ┣ .gitignore
 ┗ README.md
```

### Description of Directory Tree
- `package/` : contains all package code for gifba package
    - `package/gifba/` : contains source code (objects, methods, validation models, etc.)
    - `package/gifba/Simple_Models` : All simple models used for validation of methods
- `Examples` : contains all validation, results, \& figure generation
    - `Examples/gifba_comparisons.ipynb` : Contains validation code for the following:
        - Generating compartmentalized FBA (cFBA) models (stored in `Examples/Results/cFBA_Models`)
        - Solving cFBA models \& finding solution space 
        - Creating \& solving MICOM communities
    - `Examples/gifba_utilization.ipynb` : Contains walkthrough code for the following:
        - All 12 Simple Modesl for validation with:
            - Initialization \& giFBA usage
            - Plotting/Comparison to individualized FBA, cFBA, MICOM, \& giFBA
        - 2 simulations of real AGORA2 communities (*not included in repository*) with:
            - Summary Display
            - Cytoscape-ready files
    - `Examples/create_simple_fba.ipynb` : Notebook to create each unique "organism" for simple model validation. Some models are exact duplicates or similar (internal rxn bounds changed), thus not shown. 
    - `Examples/Model_Diagrams/` : Contains all Validation model diagrams (corresponding to simple models in `package/gifba/Simple_Models`)
    - `Examples/Results/` : Contains the following:
        - `cFBA/` : csv files of each simulation and the solution space 
        - `giFBA_AGORA/` : Cytoscape-ready files of nodes \& edges for AGORA2 model examples
        - `giFBA_Simple/` : Saved images of the simulated growth of each simple, validation model
        - `micom_results.csv` : MICOM simulation results for each simulation (each simulation has 5 tradeoff parameter values stored)
