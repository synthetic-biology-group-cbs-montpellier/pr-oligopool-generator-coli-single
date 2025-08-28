<img width="1281" height="396" alt="image" src="https://github.com/user-attachments/assets/777499ed-c558-4274-9f18-994f1eb9f47c" />


# pr-oligopool-generator-coli-single
This script analyzes combinatorial DPromoter-RBS (PR) library data from Kosuri et al. (https://doi.org/10.1073/pnas.1301301110) to filter, select, and prepare diverse sequence variants for oligopool ordering. It performs quality control, diversity analysis, and generates ready-to-order sequences with cloning sites.


## Filtering Steps
1. **Length Filtering**: Removes sequences outside specified length range
2. **Component Exclusion**: Filters out unwanted promoters and RBS elements
3. **BsaI Site Filtering**: Removes sequences with >2 BsaI restriction sites
4. **Sequence Reconstruction**: Filters sequences missing critical components
5. **Diversity Selection**: Selects variants using intensity binning and k-mer clustering

## Input files 
sd01.xlsx. sd02.xlsx, sd03.xlsx from Kosuri et al. (provided, put them in the target folder)

## User-Configurable Parameters

### Sequence Filtering
```python
MIN_SEQUENCE_LENGTH = 50        # Minimum sequence length (bp)
MAX_SEQUENCE_LENGTH = 70        # Maximum sequence length (bp)
EXCLUDE_RBS_NAMES = [...]       # RBS names to exclude (e.g., "DeadRBS")
EXCLUDE_PROMOTER_NAMES = [...]  # Promoter names to exclude (e.g., "pTEtO")
```

### Variant Selection
```python
TOTAL_VARIANTS_DESIRED = 400    # Total number of variants to select
MIN_VARIANTS_PER_BIN = 2        # Minimum variants per intensity bin
SELECTION_MODE = "proportional" # Selection strategy: "proportional", "equal", "hybrid"
```

### Analysis Parameters
```python
N_CLUSTERS = 8                  # Number of sequence clusters for diversity
n_bins = 12                     # Number of intensity bins
USE_LOG_SCALE_BINNING = True    # Log vs linear intensity binning
```

### Molecular Cloning Sites
```python
BSAI_SITE_5_PRIME = "GGTCTC"    # 5' BsaI site
BSAI_SITE_3_PRIME = "GAGACC"    # 3' BsaI site
PRIMER_SITE_5_PRIME = "TTGACA"  # 5' primer annealing site
PRIMER_SITE_3_PRIME = "TGTCAA"  # 3' primer annealing site
```

## Main Outputs

### CSV Files
- `01_length_filtered_constructs.csv` - Sequences after length filtering
- `02_constructs_after_exclusions.csv` - Sequences after component exclusions
- `03_constructs_after_bsai_filter.csv` - Sequences after BsaI filtering
- `04_constructs_after_sequence_reconstruction.csv` - Final filtered library
- `05_selected_variants.csv` - Selected diverse variants
- `selected_variants_with_sites.csv` - Variants with BsaI sites added
- `Ready_to_order_variants.csv` - Final sequences ready for synthesis

### Analysis Report
- `complete_analysis_report.pdf` - Comprehensive 5-page report including:
  - **Page 1**: Full vs Filtered Library PCA Analysis
  - **Page 2**: Cluster Family Composition Analysis
  - **Page 3**: Selected Variants Analysis (6 plots)
  - **Page 4**: Cumulative Distribution Analysis
  - **Page 5**: Ready-to-Order Sequence Length Analysis

### Additional Outputs
- Individual PNG plots for each analysis
- Clustering analysis CSVs
- Processing statistics and filtering summaries

## Usage
```bash
python PR_library_generator_single.py
```

The script automatically processes the library data, applies all filtering steps, performs diversity selection, adds cloning sites, and generates the complete analysis report. All parameters can be modified at the top of the script before execution.
