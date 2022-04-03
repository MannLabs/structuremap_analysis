# Create python environment
conda create --name structuremap_analysis python=3.8 jupyter -y
# Activate python environment
conda activate structuremap_analysis
# Install alphamap (this also installs structuremap as dependency)
pip install alphamap
# Make structuremap_analysis kernel available
conda install nb_conda_kernels -y
# Start jupyter notebook for the analysis
jupyter notebook
