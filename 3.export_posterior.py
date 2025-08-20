# Export the final h5ad file with abundance information after deconvolution by the cell2location model.
# Runs on a super-computing cluster - better on a high memory node (to load the model, or will have to reduce batch size which will increase the time).

import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
from scipy import sparse
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42



# define results directory which stores h5ad files and trained model
results_folder = "/path/to/file/Spatial_analysis_cell2location_salivary_gland"
ref_run_name = f'{results_folder}/reference_signatures' # reference regression model
run_name = f'{results_folder}/cell2location_map'



# loading the model again
# loading the adata_pre_processed
adata_pre_processed = sc.read_h5ad(f"{run_name}/spatial_adata.h5ad")
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_pre_processed)


adata_pre_processed = mod.export_posterior(
    adata_pre_processed, sample_kwargs={'num_samples': 1000, 'batch_size': 1000, 'use_gpu': False}
)


# mod.save(f"{run_name}", overwrite=True)
adata_file = f"{run_name}/spatial_adata.h5ad"
adata_pre_processed.write(adata_file)
