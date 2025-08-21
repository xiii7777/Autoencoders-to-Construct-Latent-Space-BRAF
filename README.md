# Kinase Activation Loop Conformational Modeling

This project investigates the conformational landscape of protein kinases, with a focus on the activation loop (A-loop) of B-Raf. By combining structural alignment, dimensionality reduction techniques (PCA and autoencoders), and regression modeling, the analysis identifies key geometric features that are closely associated with conformational transitions.

---

## Project Structure

The project is organized as follows:
```
├── 6UAN_chainD.pdb # Template kinase structure (B-Raf, chain D)
├── PDBs # a subset dataset for testing code execution
├── Pipeline_Alignment_MUSTANG.ipynb # Structure alignment using MUSTANG
├── Pipeline_Alignment_PYMOL.ipynb # Manual alignment using PyMOL
├── Pipeline_DataProcess.ipynb # Preprocessing of input structures and cleaning
├── Pipeline_Fitting.ipynb # Fitting or interpolation of coordinates
├── Pipeline_Autoencoder.ipynb # Training and evaluation of the autoencoder model
├── Pipeline_PCA.ipynb # PCA-based dimensionality reduction and visualization
├── Pipeline_prediction_activationloop.ipynb # Predictive modeling using DFG–APE features
├── Pipeline_Prediction_wholestructure.ipynb # Predictive modeling using whole-structure features
├── strip_pdb.py
├── structure-matching-IPR011009.tsv # InterPro-matched kinase structure list
├── template.pdb # Template structure used for fitting
├── requirements.txt # List of Python dependencies
└── README.md # Project documentation
```

---

## Dependencies

This project requires the following Python packages and external tools:

### Scientific computing & data handling
- `numpy`  
- `pandas`  
- `scipy`  
- `h5py`  
- `joblib`

### Plotting & visualization
- `matplotlib`  
- `seaborn`  
- `networkx`  
- `nglview`  
- `mpl_toolkits`

### Machine learning
- `scikit-learn`  
- `torch` (PyTorch)  
- `joblib`  
- `molearn`

### Bioinformatics & structural biology
- `biopython`  
- `MDAnalysis`  
- `mdtraj`  
- `pymol` (or `pymol-open-source`)  
- `modeller`

### System & utility
- `tqdm`  
- `requests`  
- `concurrent.futures`  
- `subprocess`  
- `multiprocessing`  
- `xml.etree.ElementTree`  
- `tempfile`  
- `html`  
- `sys`  
- `os`  
- `glob`  
- `re`  
- `pickle`  
- `urllib`  
- `collections`

---

## Usage Instructions

To reproduce the analysis pipeline, run the notebooks in the following order:

### 1. Data Preparation
- `Pipeline_DataProcess.ipynb`  
  Download and clean raw PDB files.

### 2. Structural Alignment
- `Pipeline_Alignment_PYMOL.ipynb`  
  Manual PyMOL alignment.
- `Pipeline_Alignment_MUSTANG.ipynb`  
  Automated structure alignment using MUSTANG.

### 3. Fitting and Outlier Removal
- `Pipeline_Fitting.ipynb`  
  Interpolation of coordinates and cleaning outliers.

### 4. Dimensionality Reduction
- `Pipeline_PCA.ipynb`  
  Perform PCA on fitted structures.
- `Pipeline_Autoencoder.ipynb`  
  Train a FoldingNet-based autoencoder on the activation loop.

### 5. Regression Evaluation
- `Pipeline_prediction_wholestructure.ipynb`  
  Predict PCA/latent variables using whole-structure features.
- `Pipeline_prediction_activationloop.ipynb`  
  Predict latent variables using DFG–APE features.

---

## Acknowledgement

I would like to acknowledge Sam Martino for providing the initial version of the code, which laid the foundation for this project. I am also grateful to Marco for his subsequent modifications and improvements, which greatly facilitated the later analyses.
