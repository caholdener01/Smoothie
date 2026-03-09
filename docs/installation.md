### **Recommended Install**
```bash
# 1. Create a new Conda environment (Choose favorite Python version)
conda create -n smoothie_env python=3.10 -y

# 2. Activate the environment
conda activate smoothie_env

# 3. Install the package directly from PyPI
pip install smoothie-st

# 4. (Optional) Make the environment available in Jupyter Lab
conda install ipykernel -y
python -m ipykernel install --user --name=smoothie_env --display-name="Python (smoothie_env)"
```
To use Smoothie after installation, create a new copy of one of the blank scripts from [examples](https://caholdener01.github.io/Smoothie/examples/single_dataset_pipeline/) and follow the steps in the pipeline.
Refer to the [graphical overview](https://caholdener01.github.io/Smoothie/graphical_overview/) and the supporting files from [guides](https://caholdener01.github.io/Smoothie/guides/convert_to_anndata/) to correctly and confidently use Smoothie.


### Quick Install
```bash
pip install smoothie-st
```

### Github Clone Install
```bash
git clone --depth 1 https://github.com/caholdener01/Smoothie.git
cd Smoothie
```

### Reproducibility Note
Smoothie uses the Infomap random-walk community detection algorithm to identify gene modules. Because Infomap is non-deterministic, clustering results can vary slightly across runs.

To improve reproducibility, Smoothie sets a random seed (default=0). However, exact reproducibility also depends on your Python version and installed package versions.

To reproduce the exact clustering results from the tutorials on this website, do the following installation steps:
```bash
# 1. Clone the repository
git clone https://github.com/caholdener01/Smoothie.git
cd Smoothie

# 2. Create a new Conda environment
conda create -n smoothie_tutorials_env python=3.10 -y

# 3. Activate the environment
conda activate smoothie_tutorials_env

# 4. Install pinned and flexible dependencies
pip install -r tutorials_reproducibility.txt
pip install .

# 5. (Optional) Register the environment in Jupyter Lab
conda install ipykernel -y
python -m ipykernel install --user --name=smoothie_tutorials_env --display-name="Python (smoothie_tutorials_env)"
```

### Previous Version Install
To install the previous Smoothie version, you can clone an old tagged version:
```bash
git clone -b v0.2.4 --single-branch https://github.com/caholdener01/Smoothie.git
```
Please note that Smoothie v0.2.3 had a bug in the second-order correlation gene stability analysis. Smoothie v0.2.4 fixes this bug, leaving the rest of the code unchanged.