# %% [markdown]
# # Combinatorial Library Sampling and Analysis
# This notebook walks through the process of selecting and analyzing combinatorial libraries of genetic constructs using flow-seq data and sequence features (4-mer counts, hierarchical clustering). Data are from Kosuri et al. (2013).
# 
# ## Objectives
# - Load and preprocess flow‑seq and sequence data
# - Compute the **center‑of‑mass (CM)** expression metric
# - Cluster sequences for diversity using 4‑mer counts and hierarchical clustering
# - Sample composite subsets of 200 and 1000 constructs
# - Explore and visualize distribution and diversity
# - Validate sampling balance via summary metrics
# - Conclude with interpretation and next steps

# %% [markdown]
# ## Setup and Output Directory Creation
# Create a timestamped folder for this analysis run and copy the script.

# %%
import pandas as pd
import numpy as np
import re
import os
import shutil
from datetime import datetime
import matplotlib
matplotlib.use('Agg')  # Set non-interactive backend at the very beginning
import matplotlib.pyplot as plt
plt.ioff()  # Turn off interactive mode globally

print("=== STARTING ANALYSIS ===")
print("Imports completed successfully")

# Create timestamped output directory
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
output_dir = f"analysis_results_{timestamp}"
os.makedirs(output_dir, exist_ok=True)

print(f"Analysis results will be saved to: {output_dir}")

# %% [markdown]
# ## Data Loading and Preprocessing
# Load promoter, RBS, and construct measurements (flowsorted bins) from Excel files.
# 
# ```python
# import pandas as pd
# import numpy as np
# import re
# ```
# 
# Bin counts are in columns `bin.1` to `bin.12` in the constructs dataframe.

# %%
# Load data
print("\n=== DATA LOADING ===")
try:
    print("Loading constructs from sd03.xlsx...")
    constructs = pd.read_excel('sd03.xlsx', sheet_name='Constructs', engine='openpyxl')
    print(f"Constructs loaded: {len(constructs)} rows, {len(constructs.columns)} columns")
    
    print("Loading promoters from sd01.xlsx...")
    prom_df = pd.read_excel('sd01.xlsx', sheet_name='Promoters', engine='openpyxl')
    print(f"Promoters loaded: {len(prom_df)} rows, {len(prom_df.columns)} columns")
    
    print("Loading RBS from sd02.xlsx...")
    rbs_df = pd.read_excel('sd02.xlsx', sheet_name='RBSs', engine='openpyxl')
    print(f"RBS loaded: {len(rbs_df)} rows, {len(rbs_df.columns)} columns")
    
except Exception as e:
    print(f"ERROR loading Excel files: {e}")
    print("Please check that the following files exist:")
    print("- sd03.xlsx (with 'Constructs' sheet)")
    print("- sd01.xlsx (with 'Promoters' sheet)")  
    print("- sd02.xlsx (with 'RBSs' sheet)")
    raise

# Clean names
def clean_quotes(df, col):
    """
    Remove quotes and leading/trailing whitespace from a specific column in a DataFrame.
    
    This function cleans text data by removing quotation marks and trimming whitespace,
    which is essential for consistent data processing when reading from Excel files.
    
    Args:
        df (pd.DataFrame): The DataFrame containing the column to clean
        col (str): The name of the column to clean
    
    Returns:
        None: Modifies the DataFrame in place
    """
    if col in df.columns:
        df[col] = df[col].str.replace('"', '').str.strip()
    else:
        print(f"WARNING: Column '{col}' not found in dataframe")

print("Cleaning column names...")
clean_quotes(prom_df, 'Promoter')
clean_quotes(rbs_df, 'RBS')
clean_quotes(constructs, 'Promoter')
clean_quotes(constructs, 'RBS')

# Save processed input data
constructs.to_csv(os.path.join(output_dir, "01_raw_constructs.csv"), index=False)
prom_df.to_csv(os.path.join(output_dir, "01_promoters.csv"), index=False)
rbs_df.to_csv(os.path.join(output_dir, "01_rbs.csv"), index=False)

print("Raw data saved to output directory")
print("Constructs columns:", list(constructs.columns))
print("Promoters columns:", list(prom_df.columns))
print("RBS columns:", list(rbs_df.columns))
constructs.head()

# %% [markdown]
# ## GFP Fluorescence Data Processing
# Use the 'prot' column values as GFP fluorescence measurements and filter out unmeasured data points.

# %%
# Check if 'prot' column exists
print("\n=== GFP FLUORESCENCE PROCESSING ===")
print("Checking for 'prot' column...")

if 'prot' not in constructs.columns:
    print("'prot' column not found. Available columns:", list(constructs.columns))
    # Look for similar column names
    prot_cols = [col for col in constructs.columns if 'prot' in col.lower()]
    if prot_cols:
        print(f"Found potential protein columns: {prot_cols}")
        constructs['prot'] = constructs[prot_cols[0]]
        print(f"Using {prot_cols[0]} as 'prot' column")
    else:
        raise ValueError("No 'prot' column found. Please check your data structure.")
else:
    print("'prot' column found successfully")

# Convert prot column to numeric and identify valid measurements
print("Converting 'prot' column to numeric...")
constructs['prot'] = pd.to_numeric(constructs['prot'], errors='coerce')

# Remove data points that are not measured experimentally
# Filter out NaN values and potentially invalid measurements (e.g., <= 0)
print("Filtering out unmeasured/invalid data...")
initial_count = len(constructs)
constructs = constructs[constructs['prot'].notna()].copy()
constructs = constructs[constructs['prot'] > 0].copy()  # Remove zero or negative values
final_count = len(constructs)

print(f"Data filtering:")
print(f"  Initial constructs: {initial_count}")
print(f"  After removing unmeasured/invalid data: {final_count}")
print(f"  Removed: {initial_count - final_count} constructs")

if final_count == 0:
    raise ValueError("No valid protein measurements found after filtering!")

# Use prot column directly as GFP fluorescence
constructs['GFP_fluorescence'] = constructs['prot']

# Calculate additional metrics for analysis
constructs['log_GFP'] = np.log10(constructs['GFP_fluorescence'])

print(f"GFP fluorescence statistics:")
print(f"  Mean: {constructs['GFP_fluorescence'].mean():.3f}")
print(f"  Std: {constructs['GFP_fluorescence'].std():.3f}")
print(f"  Range: {constructs['GFP_fluorescence'].min():.3f} - {constructs['GFP_fluorescence'].max():.3f}")
print(f"  Log10 range: {constructs['log_GFP'].min():.3f} - {constructs['log_GFP'].max():.3f}")

# Save constructs with GFP fluorescence data
constructs.to_csv(os.path.join(output_dir, "02_constructs_with_GFP.csv"), index=False)
print("Constructs with GFP fluorescence data saved")

# %% [markdown]
# ## Sequence Reconstruction and Processing
# Reconstruct full sequences by looking up promoter and RBS sequences, then trim to functional regions.

# %%
# Step 3: Reconstitute concatenated sequences using promoter and RBS identifiers
print("\n=== SEQUENCE RECONSTRUCTION ===")

# Clean promoter and RBS names to ensure matching
def clean_name(name):
    """
    Clean individual name values by removing quotes and whitespace.
    
    This function standardizes individual name strings for consistent matching
    between different data sources (constructs, promoters, RBS datasets).
    Handles NaN values gracefully by returning them unchanged.
    
    Args:
        name (str or NaN): The name string to clean
    
    Returns:
        str or NaN: Cleaned name string, or original NaN if input was NaN
    """
    if pd.isna(name):
        return name
    return str(name).replace('"', '').strip()

print("Cleaning promoter and RBS names...")
# Clean all names
constructs['Promoter'] = constructs['Promoter'].apply(clean_name)
constructs['RBS'] = constructs['RBS'].apply(clean_name)
prom_df['Promoter'] = prom_df['Promoter'].apply(clean_name)
rbs_df['RBS'] = rbs_df['RBS'].apply(clean_name)

def clean_sequence(seq):
    """
    Clean DNA sequence data by removing quotes and whitespace.
    
    This function standardizes DNA sequence strings by removing both double and single quotes,
    as well as leading/trailing whitespace. Essential for consistent sequence processing
    when working with data imported from Excel files.
    
    Args:
        seq (str or NaN): The DNA sequence string to clean
    
    Returns:
        str or NaN: Cleaned sequence string, or original NaN if input was NaN
    """
    if pd.isna(seq):
        return seq
    return str(seq).replace('"', '').replace("'", '').strip()

# Clean sequences before creating lookup dictionaries
print("Cleaning sequence data...")
prom_df['Sequence'] = prom_df['Sequence'].apply(clean_sequence)
rbs_df['Sequence'] = rbs_df['Sequence'].apply(clean_sequence)

# Create lookup dictionaries
print("Creating lookup dictionaries...")
promoter_sequences = dict(zip(prom_df['Promoter'], prom_df['Sequence']))
rbs_sequences = dict(zip(rbs_df['RBS'], rbs_df['Sequence']))

print(f"Available promoters: {len(promoter_sequences)}")
print(f"Available RBS: {len(rbs_sequences)}")

# Check for potential sequence column issues
if 'Sequence' not in prom_df.columns:
    print("WARNING: 'Sequence' column not found in promoters dataframe")
    print("Promoter columns:", list(prom_df.columns))
if 'Sequence' not in rbs_df.columns:
    print("WARNING: 'Sequence' column not found in RBS dataframe")
    print("RBS columns:", list(rbs_df.columns))

# Reconstruct sequences
print("Reconstructing sequences...")
def reconstruct_sequence(row):
    """
    Reconstruct full DNA sequences by concatenating promoter and RBS sequences.
    
    This function takes a construct row with promoter and RBS identifiers,
    looks up their corresponding DNA sequences from dictionaries, and
    concatenates them to create the full sequence for the construct.
    
    Args:
        row (pd.Series): A DataFrame row containing 'Promoter' and 'RBS' columns
    
    Returns:
        str or None: Concatenated DNA sequence (promoter + RBS), or None if either sequence is missing
    """
    prom_seq = promoter_sequences.get(row['Promoter'])
    rbs_seq = rbs_sequences.get(row['RBS'])
    
    if pd.isna(prom_seq) or pd.isna(rbs_seq):
        return None
    
    return str(prom_seq) + str(rbs_seq)

# Apply reconstruction with progress tracking
print(f"Processing {len(constructs)} constructs...")
constructs['full_sequence'] = constructs.apply(reconstruct_sequence, axis=1)

# Remove constructs without reconstructed sequences
initial_seq_count = len(constructs)
constructs = constructs[constructs['full_sequence'].notna()].copy()
final_seq_count = len(constructs)

# IMPORTANT: Create a simple sequential index for PCA mapping
# Since PCA is performed on this exact dataset (12048 constructs), 
# we just need to track positions in this dataset
constructs = constructs.reset_index(drop=True)
constructs['pca_index'] = constructs.index  # This will be 0, 1, 2, ..., 12047

print(f"Sequence reconstruction:")
print(f"  Successfully reconstructed: {final_seq_count}")
print(f"  Failed to reconstruct: {initial_seq_count - final_seq_count}")

if final_seq_count == 0:
    raise ValueError("No sequences could be reconstructed! Check promoter/RBS matching.")

# Step 4: Trim sequences between GGCGCGCC (5') and CATATG (3')
print("Trimming sequences...")
def trim_sequence(seq):
    """
    Trim DNA sequences to extract functional regions between specific motifs.
    
    This function extracts the functional region of a DNA sequence by trimming
    everything before the GGCGCGCC motif (5' end) and everything after the 
    CATATG motif (3' end). If either motif is not found, it uses the sequence
    boundaries as fallback positions.
    
    Args:
        seq (str or NaN): The full DNA sequence to trim
    
    Returns:
        str or None: Trimmed sequence containing only the functional region,
                     or None if input was NaN or result is empty
    """
    if pd.isna(seq):
        return None
    
    seq = str(seq).upper()
    start_motif = "GGCGCGCC"
    end_motif = "CATATG"
    
    # Find start position (after the motif)
    start_pos = seq.find(start_motif)
    if start_pos == -1:
        # If exact motif not found, skip trimming from 5'
        start_pos = 0
    else:
        start_pos += len(start_motif)
    
    # Find end position (before the motif)
    end_pos = seq.find(end_motif, start_pos)
    if end_pos == -1:
        # If exact motif not found, keep until end
        end_pos = len(seq)
    
    trimmed = seq[start_pos:end_pos]
    return trimmed if trimmed else None

constructs['trimmed_sequence'] = constructs['full_sequence'].apply(trim_sequence)

# Remove constructs with empty trimmed sequences
constructs = constructs[constructs['trimmed_sequence'].notna()].copy()
constructs = constructs[constructs['trimmed_sequence'].str.len() > 0].copy()
trimmed_count = len(constructs)

print(f"Sequence trimming:")
print(f"  Sequences after trimming: {trimmed_count}")
if trimmed_count > 0:
    print(f"  Average trimmed length: {constructs['trimmed_sequence'].str.len().mean():.1f} bp")
    print(f"  Length range: {constructs['trimmed_sequence'].str.len().min()} - {constructs['trimmed_sequence'].str.len().max()} bp")

# Add sequence length column
constructs['sequence_length'] = constructs['trimmed_sequence'].str.len()

# Step 6: K-mer and PCA analysis for diversity assessment (BEFORE length filtering)
print("\n=== K-MER DIVERSITY ANALYSIS (FULL LIBRARY) ===")

# User-configurable k-mer analysis parameters
KMER_SIZE = 5  # User can change this to 5 for 5-mer analysis (4=256 features, 5=1024 features)
print(f"Using k-mer size: {KMER_SIZE}")
print("Performing k-mer analysis on FULL LIBRARY before length filtering...")

# Function to count k-mers
def kmer_counts(seq, k=KMER_SIZE):
    counts = Counter([seq[i:i+k] for i in range(len(seq)-k+1)])
    return counts

# Create k-mer feature matrix
from itertools import product
from collections import Counter
kmers = [''.join(p) for p in product('ACGT', repeat=KMER_SIZE)]
print(f"Analyzing {len(kmers)} {KMER_SIZE}-mers for {len(constructs)} constructs...")

# Performance warning for 5-mers
if KMER_SIZE == 5:
    print(f"WARNING: Using 5-mers will create {len(kmers)} features and may be slow for large datasets")
    print("Consider using k=4 (256 features) for faster processing, or limit dataset size")

kmat_full = np.zeros((len(constructs), len(kmers)))
print("Building k-mer matrix for full library...")

# Process in batches to show progress
batch_size = 100
try:
    for i in range(0, len(constructs), batch_size):
        end_i = min(i + batch_size, len(constructs))
        if i % 1000 == 0 or i == 0:
            print(f"  Processing constructs {i+1}-{end_i} ({i/len(constructs)*100:.1f}%)")
        
        for j in range(i, end_i):
            seq = constructs.iloc[j]['trimmed_sequence']
            if pd.notna(seq) and len(str(seq)) > 0:
                cnts = kmer_counts(str(seq))
                for k_idx, k in enumerate(kmers):
                    kmat_full[j, k_idx] = cnts.get(k, 0)
    
    print("K-mer matrix completed successfully")
    
except Exception as e:
    print(f"Error during k-mer analysis: {e}")
    print("Attempting to continue with available data...")

# Normalize k-mer matrix
print("Normalizing k-mer matrix...")
try:
    # Check for empty rows (sequences with no valid k-mers)
    row_sums = kmat_full.sum(axis=1)
    valid_rows = row_sums > 0
    
    if valid_rows.sum() < len(constructs):
        print(f"Warning: {len(constructs) - valid_rows.sum()} sequences have no valid k-mers")
    
    kmat_full_norm = np.zeros_like(kmat_full)
    kmat_full_norm[valid_rows] = kmat_full[valid_rows] / (row_sums[valid_rows, np.newaxis] + 1e-10)
    
    print("K-mer normalization completed")
except Exception as e:
    print(f"Error during normalization: {e}")
    # Fallback: use simple normalization
    kmat_full_norm = kmat_full / (kmat_full.sum(axis=1, keepdims=True) + 1e-10)

# PCA Analysis for Sequence Diversity Assessment (FULL LIBRARY)
print("\n=== PCA ANALYSIS FOR SEQUENCE DIVERSITY (FULL LIBRARY) ===")
try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.cluster import KMeans
    
    # Perform PCA on k-mer features
    print("Performing PCA on k-mer features (FULL LIBRARY)...")
    scaler_full = StandardScaler()
    kmat_full_scaled = scaler_full.fit_transform(kmat_full_norm)
    
    # PCA with components explaining 95% variance
    pca_full = PCA(n_components=0.95)
    pca_features_full = pca_full.fit_transform(kmat_full_scaled)
    
    print(f"PCA completed: {pca_full.n_components_} components explain {pca_full.explained_variance_ratio_.sum():.3f} of variance")
    
    # Add PCA coordinates to constructs (FULL LIBRARY)
    for i in range(min(5, pca_full.n_components_)):  # Store first 5 PC coordinates
        constructs[f'PC{i+1}'] = pca_features_full[:, i]
    
    # K-means clustering for sequence grouping (FULL LIBRARY)
    print("\n=== SEQUENCE CLUSTERING ANALYSIS (FULL LIBRARY) ===")
    N_CLUSTERS = 8  # User can modify this - number of sequence clusters to identify
    print(f"Performing k-means clustering with {N_CLUSTERS} clusters on FULL LIBRARY...")
    
    # Use first few PC components for clustering (or all if fewer available)
    n_components_for_clustering = min(10, pca_full.n_components_)
    clustering_features_full = pca_features_full[:, :n_components_for_clustering]
    
    # Perform k-means clustering
    kmeans_full = KMeans(n_clusters=N_CLUSTERS, random_state=42, n_init=10)
    cluster_labels_full = kmeans_full.fit_predict(clustering_features_full)
    
    # Add cluster labels to constructs (FULL LIBRARY)
    constructs['sequence_cluster'] = cluster_labels_full
    
    print(f"Clustering completed. Cluster distribution (FULL LIBRARY):")
    cluster_counts = pd.Series(cluster_labels_full).value_counts().sort_index()
    for cluster_id, count in cluster_counts.items():
        print(f"  Cluster {cluster_id}: {count} sequences ({count/len(constructs)*100:.1f}%)")
    
    # Analyze cluster characteristics (FULL LIBRARY)
    print("\nCluster characteristics (FULL LIBRARY):")
    cluster_stats_full = []
    cluster_families_analysis = []
    
    for cluster_id in range(N_CLUSTERS):
        cluster_mask = constructs['sequence_cluster'] == cluster_id
        cluster_data = constructs[cluster_mask]
        
        if len(cluster_data) > 0:
            stats = {
                'Cluster': cluster_id,
                'Count': len(cluster_data),
                'Mean_Intensity': cluster_data['GFP_fluorescence'].mean(),
                'Std_Intensity': cluster_data['GFP_fluorescence'].std(),
                'Mean_Length': cluster_data['sequence_length'].mean(),
                'Std_Length': cluster_data['sequence_length'].std(),
                'Intensity_Range': f"{cluster_data['GFP_fluorescence'].min():.0f}-{cluster_data['GFP_fluorescence'].max():.0f}"
            }
            cluster_stats_full.append(stats)
            
            print(f"  Cluster {cluster_id}: {len(cluster_data):4d} seqs, "
                  f"intensity {cluster_data['GFP_fluorescence'].mean():8.0f}±{cluster_data['GFP_fluorescence'].std():6.0f}, "
                  f"length {cluster_data['sequence_length'].mean():.1f}±{cluster_data['sequence_length'].std():.1f}")
            
            # SEQUENCE FAMILY ANALYSIS
            print(f"    Analyzing sequence families in Cluster {cluster_id}...")
            
            # Function to extract family patterns from names
            def extract_family_pattern(name, component_type='unknown'):
                """
                Extract standardized family patterns from promoter and RBS component names.
                
                This function categorizes genetic components into families based on naming conventions
                and known biological part collections. It's used for analyzing the diversity and
                distribution of different component families within sequence clusters.
                
                Args:
                    name (str or NaN): The component name to analyze
                    component_type (str): Type of component ('rbs', 'promoter', or 'unknown')
                
                Returns:
                    str: Standardized family name for the component
                """
                if pd.isna(name):
                    return "Unknown"
                name_str = str(name).upper()
                
                # RBS-specific patterns  
                if component_type == 'rbs':
                    # BBa_J611XX and Anderson are Anderson RBS
                    if 'J611' in name_str or 'ANDERSON' in name_str:
                        return 'Anderson RBS'
                    # B003X are BioBrick RBS
                    elif name_str.startswith('B003') or 'B003' in name_str:
                        return 'BioBrick RBS'
                    # apFAB are BioFAB RBS
                    elif 'APFAB' in name_str:
                        return 'BioFAB RBS'
                    # Salis RBS (contains "salis" in name)
                    elif 'SALIS' in name_str:
                        return 'Salis RBS'
                    # All others classified as Others
                    else:
                        return 'Other RBS'
                
                # Promoter-specific patterns
                elif component_type == 'promoter':
                    # BBa_J231XX are Anderson promoters
                    if 'J231' in name_str or ('J23' in name_str and '1' in name_str):
                        return 'Anderson Promoters'
                    # apFAB are BioFAB promoters
                    elif 'APFAB' in name_str:
                        return 'BioFAB Promoters'
                    # For other classes of promoters: use the name in the column
                    else:
                        # Clean up the name and return as-is for other promoter classes
                        cleaned_name = str(name).strip().replace('"', '')
                        return cleaned_name
                
                # Generic pattern (fallback)
                else:
                    if len(name_str) >= 3:
                        return name_str[:3]
                    else:
                        return name_str
            
            # Analyze promoter families
            if 'Promoter' in cluster_data.columns:
                promoter_families = cluster_data['Promoter'].apply(lambda x: extract_family_pattern(x, 'promoter'))
                promoter_counts = promoter_families.value_counts()
                total_promoters = len(cluster_data)
                
                print(f"    Promoter families (≥2%):")
                # Show all families with ≥2% representation
                significant_promoters = promoter_counts[promoter_counts / total_promoters >= 0.02]
                for family, count in significant_promoters.items():
                    pct = (count / total_promoters) * 100
                    print(f"      {family}: {count:3d} ({pct:5.1f}%)")
                
                # If no families reach 2%, show top 3
                if len(significant_promoters) == 0:
                    print(f"    Top promoter families (none ≥2%):")
                    for family, count in promoter_counts.head(3).items():
                        pct = (count / total_promoters) * 100
                        print(f"      {family}: {count:3d} ({pct:5.1f}%)")
                
                # Store top promoter family
                top_promoter_family = promoter_counts.index[0] if len(promoter_counts) > 0 else "None"
                top_promoter_pct = (promoter_counts.iloc[0] / total_promoters * 100) if len(promoter_counts) > 0 else 0
                
                # Store all significant promoter families for detailed analysis
                promoter_families_detail = []
                for family, count in significant_promoters.items():
                    pct = (count / total_promoters) * 100
                    promoter_families_detail.append(f"{family}: {pct:.1f}%")
                if not promoter_families_detail:  # If none ≥2%, use top family
                    promoter_families_detail.append(f"{top_promoter_family}: {top_promoter_pct:.1f}%")
                
            else:
                top_promoter_family = "N/A"
                top_promoter_pct = 0
                promoter_families_detail = ["N/A"]
            
            # Analyze RBS families
            if 'RBS' in cluster_data.columns:
                rbs_families = cluster_data['RBS'].apply(lambda x: extract_family_pattern(x, 'rbs'))
                rbs_counts = rbs_families.value_counts()
                total_rbs = len(cluster_data)
                
                print(f"    RBS families (≥2%):")
                # Show all families with ≥2% representation
                significant_rbs = rbs_counts[rbs_counts / total_rbs >= 0.02]
                for family, count in significant_rbs.items():
                    pct = (count / total_rbs) * 100
                    print(f"      {family}: {count:3d} ({pct:5.1f}%)")
                
                # If no families reach 2%, show top 3
                if len(significant_rbs) == 0:
                    print(f"    Top RBS families (none ≥2%):")
                    for family, count in rbs_counts.head(3).items():
                        pct = (count / total_rbs) * 100
                        print(f"      {family}: {count:3d} ({pct:5.1f}%)")
                
                # Store top RBS family
                top_rbs_family = rbs_counts.index[0] if len(rbs_counts) > 0 else "None"
                top_rbs_pct = (rbs_counts.iloc[0] / total_rbs * 100) if len(rbs_counts) > 0 else 0
                
                # Store all significant RBS families for detailed analysis
                rbs_families_detail = []
                for family, count in significant_rbs.items():
                    pct = (count / total_rbs) * 100
                    rbs_families_detail.append(f"{family}: {pct:.1f}%")
                if not rbs_families_detail:  # If none ≥2%, use top family
                    rbs_families_detail.append(f"{top_rbs_family}: {top_rbs_pct:.1f}%")
                
            else:
                top_rbs_family = "N/A"
                top_rbs_pct = 0
                rbs_families_detail = ["N/A"]
            
            # Store family analysis
            family_analysis = {
                'Cluster': cluster_id,
                'Count': len(cluster_data),
                'Top_Promoter_Family': top_promoter_family,
                'Top_Promoter_Percentage': top_promoter_pct,
                'Top_RBS_Family': top_rbs_family,
                'Top_RBS_Percentage': top_rbs_pct,
                'Mean_Intensity': cluster_data['GFP_fluorescence'].mean(),
                'Promoter_Diversity': len(promoter_counts) if 'Promoter' in cluster_data.columns else 0,
                'RBS_Diversity': len(rbs_counts) if 'RBS' in cluster_data.columns else 0,
                'Promoter_Families_Detail': '; '.join(promoter_families_detail),
                'RBS_Families_Detail': '; '.join(rbs_families_detail)
            }
            cluster_families_analysis.append(family_analysis)
            
            print(f"    Summary:")
            print(f"      Promoters: {'; '.join(promoter_families_detail)}")
            print(f"      RBS: {'; '.join(rbs_families_detail)}")
            print()
    
    # Save cluster analysis (FULL LIBRARY)
    cluster_stats_full_df = pd.DataFrame(cluster_stats_full)
    cluster_stats_full_df.to_csv(os.path.join(output_dir, "04_sequence_clusters_full_library.csv"), index=False)
    
    # Save family analysis
    cluster_families_df = pd.DataFrame(cluster_families_analysis)
    cluster_families_df.to_csv(os.path.join(output_dir, "04_cluster_family_analysis.csv"), index=False)
    
    print("Sequence clustering analysis saved (FULL LIBRARY)")
    print("Cluster family analysis saved")
    
    print("Full library k-mer and PCA analysis completed")
    
except ImportError:
    print("scikit-learn not available - skipping PCA analysis")
    pca_full = None
    cluster_stats_full_df = None
except Exception as e:
    print(f"Error in PCA analysis: {e}")
    pca_full = None
    cluster_stats_full_df = None

# Save full library data BEFORE length filtering for visualization
print("Saving full library data before length filtering...")
full_library_constructs = constructs.copy()  # This contains all ~12,000 constructs with PCA coordinates
full_library_constructs.to_csv(os.path.join(output_dir, "03_full_library_with_pca.csv"), index=False)
print(f"Full library saved: {len(full_library_constructs)} constructs with PCA coordinates")

# User-configurable sequence length filtering
MIN_SEQUENCE_LENGTH = 50  # Your minimum length
MAX_SEQUENCE_LENGTH = 70  # Your maximum length

# User-configurable sequence exclusion by component names
EXCLUDE_RBS_NAMES = [
   # "DeadRBS"
]

EXCLUDE_PROMOTER_NAMES = [
   # "No promoter",
    "PLTETo1", 
    "pLlacO1",
    "pBAD"
]

# Set to True to enable name-based exclusion, False to disable
ENABLE_NAME_BASED_EXCLUSION = True

print(f"Filtering sequences to length range {MIN_SEQUENCE_LENGTH}-{MAX_SEQUENCE_LENGTH}bp...")
initial_count_before_length_filter = len(constructs)

# Apply length range filter
constructs = constructs[
    (constructs['sequence_length'] >= MIN_SEQUENCE_LENGTH) & 
    (constructs['sequence_length'] <= MAX_SEQUENCE_LENGTH)
].copy()

final_count_after_length_filter = len(constructs)

print(f"Length filtering:")
print(f"  Before length filter: {initial_count_before_length_filter}")
print(f"  After length filter ({MIN_SEQUENCE_LENGTH}-{MAX_SEQUENCE_LENGTH}bp): {final_count_after_length_filter}")
print(f"  Removed: {initial_count_before_length_filter - final_count_after_length_filter} sequences")

if final_count_after_length_filter == 0:
    raise ValueError(f"No sequences remain after filtering for length {MIN_SEQUENCE_LENGTH}-{MAX_SEQUENCE_LENGTH}bp!")

# Name-based exclusion filtering
if ENABLE_NAME_BASED_EXCLUSION:
    print(f"\n=== NAME-BASED EXCLUSION FILTERING ===")
    print(f"Excluding sequences with specific component names...")
    
    initial_count_before_name_filter = len(constructs)
    
    # Function to check if a name should be excluded (case-insensitive partial matching)
    def should_exclude_sequence(row):
        """
        Determine if a construct should be excluded based on component names.
        
        This function checks if a construct contains promoter or RBS components that
        are on the exclusion lists. It performs case-insensitive partial matching
        to catch variations in naming conventions.
        
        Args:
            row (pd.Series): A DataFrame row containing 'Promoter' and 'RBS' columns
        
        Returns:
            tuple: (should_exclude: bool, reason: str) where should_exclude indicates
                   if the sequence should be excluded and reason provides the specific
                   component and exclusion match that triggered the exclusion
        """
        promoter_name = str(row.get('Promoter', '')).strip()
        rbs_name = str(row.get('RBS', '')).strip()
        
        # Check promoter exclusions (case-insensitive, partial matching)
        for exclude_name in EXCLUDE_PROMOTER_NAMES:
            if exclude_name.lower() in promoter_name.lower():
                return True, f"Promoter: {promoter_name} (matches {exclude_name})"
        
        # Check RBS exclusions (case-insensitive, partial matching)
        for exclude_name in EXCLUDE_RBS_NAMES:
            if exclude_name.lower() in rbs_name.lower():
                return True, f"RBS: {rbs_name} (matches {exclude_name})"
        
        return False, ""
    
    # Apply exclusion filtering
    exclusion_reasons = []
    keep_mask = []
    
    for idx, row in constructs.iterrows():
        should_exclude, reason = should_exclude_sequence(row)
        keep_mask.append(not should_exclude)
        if should_exclude:
            exclusion_reasons.append(reason)
    
    # Filter constructs
    constructs = constructs[keep_mask].copy()
    final_count_after_name_filter = len(constructs)
    
    print(f"Name-based filtering:")
    print(f"  Before name filter: {initial_count_before_name_filter}")
    print(f"  After name filter: {final_count_after_name_filter}")
    print(f"  Removed: {initial_count_before_name_filter - final_count_after_name_filter} sequences")
    
    if initial_count_before_name_filter - final_count_after_name_filter > 0:
        print(f"Exclusion breakdown:")
        # Count exclusion reasons
        from collections import Counter
        reason_counts = Counter(exclusion_reasons)
        for reason, count in reason_counts.items():
            print(f"  - {reason}: {count} sequences")
    
    if final_count_after_name_filter == 0:
        raise ValueError("No sequences remain after name-based exclusion filtering!")
    
    print(f"Excluded components:")
    print(f"  Promoters: {', '.join(EXCLUDE_PROMOTER_NAMES)}")
    print(f"  RBS: {', '.join(EXCLUDE_RBS_NAMES)}")
else:
    print(f"\n=== NAME-BASED EXCLUSION DISABLED ===")
    print("Set ENABLE_NAME_BASED_EXCLUSION = True to enable component name filtering")

# Note: pca_index already stored after sequence reconstruction and maps directly to PCA space

# Show length distribution after filtering
length_counts = constructs['sequence_length'].value_counts().sort_index()
print(f"Length distribution after filtering:")
for length, count in length_counts.items():
    print(f"  {length:2d}bp: {count:4d} sequences ({count/len(constructs)*100:5.1f}%)")

# Update sequence length statistics after filtering
print(f"Updated sequence length range: {constructs['sequence_length'].min()} - {constructs['sequence_length'].max()} bp")

# =============================================================================
# COMPREHENSIVE FILTERING SUMMARY
# =============================================================================
print(f"\n" + "="*80)
print(f"COMPREHENSIVE FILTERING SUMMARY")
print(f"="*80)

# Calculate filtering steps from stored counts
data_filter_removed = initial_count - final_count  # From line 148-156
seq_reconstruction_removed = initial_seq_count - final_seq_count  # From line 279-291
trimming_removed = final_seq_count - trimmed_count  # From line 340-343
length_filter_removed = initial_count_before_length_filter - final_count_after_length_filter  # From line 693-710
if ENABLE_NAME_BASED_EXCLUSION:
    name_filter_removed = initial_count_before_name_filter - final_count_after_name_filter  # From line 712-787
else:
    name_filter_removed = 0

current_count = len(constructs)
total_removed = initial_count - current_count

print(f"1. Initial constructs loaded:                     {initial_count:,}")
print(f"2. After removing unmeasured/invalid data:        {final_count:,} (removed: {data_filter_removed:,})")
print(f"3. After sequence reconstruction:                 {final_seq_count:,} (removed: {seq_reconstruction_removed:,})")
print(f"4. After trimming to functional regions:          {trimmed_count:,} (removed: {trimming_removed:,})")
print(f"5. After length filtering ({MIN_SEQUENCE_LENGTH}-{MAX_SEQUENCE_LENGTH}bp):              {final_count_after_length_filter:,} (removed: {length_filter_removed:,})")
if ENABLE_NAME_BASED_EXCLUSION:
    print(f"6. After name-based exclusion (dead RBS, etc.):   {final_count_after_name_filter:,} (removed: {name_filter_removed:,})")
    print(f"   - Excluded dead RBS: {EXCLUDE_RBS_NAMES}")
    print(f"   - Excluded promoters: {EXCLUDE_PROMOTER_NAMES}")
else:
    print(f"6. Name-based exclusion: DISABLED")

print(f"\nFINAL COUNT READY FOR BINNING & SELECTION:        {current_count:,}")
print(f"TOTAL FILTERED OUT:                               {total_removed:,} ({(total_removed/initial_count)*100:.1f}%)")
print(f"SURVIVAL RATE:                                    {(current_count/initial_count)*100:.1f}%")

# Breakdown by exclusion category (if enabled)
if ENABLE_NAME_BASED_EXCLUSION and 'exclusion_reasons' in locals() and len(exclusion_reasons) > 0:
    print(f"\nDEAD/UNWANTED COMPONENT BREAKDOWN:")
    reason_counts = Counter(exclusion_reasons)
    for reason, count in reason_counts.items():
        print(f"  - {reason}: {count:,} sequences")

print(f"="*80)

# Save processed sequences
constructs.to_csv(os.path.join(output_dir, "03_sequences_reconstructed.csv"), index=False)
print("Reconstructed, trimmed, and length-filtered sequences saved")

# %% [markdown]
# ## Intensity-Based Binning and Diverse Selection
# Create bins based on protein fluorescence intensity and select diverse variants within each bin.

# %%
from itertools import product
from collections import Counter
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
import random

# Step 5: Establish 12 bins based on protein fluorescence intensity
print("=== INTENSITY BINNING ===")

n_bins = 12
USE_LOG_SCALE_BINNING = True  # User can change this - True for log10 scale, False for linear scale

if USE_LOG_SCALE_BINNING:
    print("Using log10 scale binning for better dynamic range coverage")
    # Create bins on log10 scale
    log_min = np.log10(constructs['GFP_fluorescence'].min())
    log_max = np.log10(constructs['GFP_fluorescence'].max())
    log_bins = np.logspace(log_min, log_max, n_bins + 1)
    constructs['intensity_bin'] = pd.cut(constructs['GFP_fluorescence'], 
                                       bins=log_bins, 
                                       labels=range(1, n_bins+1),
                                       include_lowest=True)
    print(f"Log10 bins range: {log_min:.2f} to {log_max:.2f}")
else:
    print("Using linear scale binning")
    constructs['intensity_bin'] = pd.cut(constructs['GFP_fluorescence'], 
                                       bins=n_bins, 
                                       labels=range(1, n_bins+1),
                                       include_lowest=True)

# Print bin statistics
print("Bin boundaries and counts:")
bin_info = []
total_constructs = len(constructs)
for i in range(1, n_bins+1):
    bin_data = constructs[constructs['intensity_bin'] == i]
    if len(bin_data) > 0:
        min_val = bin_data['GFP_fluorescence'].min()
        max_val = bin_data['GFP_fluorescence'].max()
        count = len(bin_data)
        percentage = (count / total_constructs) * 100
        bin_info.append({
            'Bin': i,
            'Count': count,
            'Percentage': percentage,
            'Min_Intensity': min_val,
            'Max_Intensity': max_val,
            'Mean_Intensity': bin_data['GFP_fluorescence'].mean(),
            'Binning_Scale': 'Log10' if USE_LOG_SCALE_BINNING else 'Linear'
        })
        if USE_LOG_SCALE_BINNING:
            print(f"  Bin {i:2d}: {count:4d} variants ({percentage:5.1f}%), intensity {min_val:8.1f}-{max_val:8.1f} (log scale)")
        else:
            print(f"  Bin {i:2d}: {count:4d} variants ({percentage:5.1f}%), intensity {min_val:8.1f}-{max_val:8.1f} (linear scale)")

# Analyze the overall distribution
print(f"\n=== INTENSITY DISTRIBUTION ANALYSIS ===")
print(f"Total constructs: {total_constructs}")
print(f"Intensity range: {constructs['GFP_fluorescence'].min():.1f} - {constructs['GFP_fluorescence'].max():.1f}")
print(f"Median intensity: {constructs['GFP_fluorescence'].median():.1f}")
print(f"Mean intensity: {constructs['GFP_fluorescence'].mean():.1f}")

# Analyze specific ranges
ranges_of_interest = [
    (0, 20000, "Very Low"),
    (20000, 50000, "Low"), 
    (50000, 100000, "Medium"),
    (100000, 150000, "High"),
    (150000, float('inf'), "Very High")
]

print(f"\n=== VARIANTS BY INTENSITY RANGE ===")
for min_val, max_val, label in ranges_of_interest:
    if max_val == float('inf'):
        mask = constructs['GFP_fluorescence'] >= min_val
        range_str = f"≥{min_val:,}"
    else:
        mask = (constructs['GFP_fluorescence'] >= min_val) & (constructs['GFP_fluorescence'] < max_val)
        range_str = f"{min_val:,}-{max_val:,}"
    
    count = mask.sum()
    percentage = (count / total_constructs) * 100
    print(f"  {label:10s} ({range_str:15s}): {count:4d} variants ({percentage:5.1f}%)")

bin_info_df = pd.DataFrame(bin_info)
bin_info_df.to_csv(os.path.join(output_dir, "04_bin_information.csv"), index=False)

# Create early visualization of intensity distribution
plt.figure(figsize=(15, 10))

# Log-scale histogram to better see the distribution
plt.subplot(2, 2, 1)
plt.hist(constructs['GFP_fluorescence'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
plt.xlabel('GFP Fluorescence')
plt.ylabel('Count')
plt.title('Intensity Distribution (Linear Scale)')
plt.yscale('log')

# Linear scale
plt.subplot(2, 2, 2)
plt.hist(constructs['GFP_fluorescence'], bins=50, alpha=0.7, color='lightcoral', edgecolor='black')
plt.xlabel('GFP Fluorescence')
plt.ylabel('Count')
plt.title('Intensity Distribution (Linear Scale)')

# Log intensity histogram
plt.subplot(2, 2, 3)
plt.hist(constructs['log_GFP'], bins=50, alpha=0.7, color='lightgreen', edgecolor='black')
plt.xlabel('Log10(GFP Fluorescence)')
plt.ylabel('Count')
plt.title('Log10 Intensity Distribution')

# Binning visualization
plt.subplot(2, 2, 4)
bin_counts = [len(constructs[constructs['intensity_bin'] == i]) for i in range(1, n_bins+1)]
plt.bar(range(1, n_bins+1), bin_counts, alpha=0.7, color='orange', edgecolor='black')
plt.xlabel('Intensity Bin')
plt.ylabel('Count')
plt.title('Variants per Intensity Bin')
plt.xticks(range(1, n_bins+1))

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "04_intensity_distribution_analysis.png"), dpi=300, bbox_inches='tight')
plt.close()  # Close the figure instead of showing it

print("Intensity distribution analysis saved")

# Note: K-mer and PCA analysis now performed on FULL LIBRARY before length filtering
# This ensures PCA space represents the complete sequence diversity

# Analyze sequence diversity within intensity bins (using filtered data)
print("\n=== ANALYZING FILTERED DATA IN FULL LIBRARY PCA SPACE ===")
if 'PC1' in constructs.columns:
    print("Analyzing sequence diversity within intensity bins (filtered data in full PCA space)...")
    bin_diversity_metrics = []
    
    # Get PCA indices of filtered constructs for mapping to full library PCA space
    filtered_pca_indices = constructs['pca_index'].tolist()
    print(f"Mapping {len(filtered_pca_indices)} filtered constructs to full library PCA space")
    print(f"Debug: PCA features shape: {pca_features_full.shape}")
    print(f"Debug: Min PCA index: {min(filtered_pca_indices)}")
    print(f"Debug: Max PCA index: {max(filtered_pca_indices)}")
    print(f"Debug: First few PCA indices: {filtered_pca_indices[:10]}")
    
    for bin_num in range(1, n_bins+1):
        bin_mask = constructs['intensity_bin'] == bin_num
        bin_data = constructs[bin_mask]
        
        if len(bin_data) < 2:
            continue
            
        # Get PCA features for this bin from the full library PCA using PCA indices
        bin_pca_indices = bin_data['pca_index'].tolist()
        bin_pca = pca_features_full[bin_pca_indices]
        
        # Calculate diversity metrics
        from scipy.spatial.distance import pdist
        
        # Mean pairwise distance in PCA space
        pairwise_distances = pdist(bin_pca[:, :min(10, pca_full.n_components_)], metric='euclidean')
        mean_distance = np.mean(pairwise_distances)
        std_distance = np.std(pairwise_distances)
        
        # Variance in each PC dimension
        pc_variances = np.var(bin_pca[:, :min(5, pca_full.n_components_)], axis=0)
        
        bin_diversity_metrics.append({
            'Bin': bin_num,
            'N_Variants': len(bin_data),
            'Mean_PCA_Distance': mean_distance,
            'Std_PCA_Distance': std_distance,
            'PC1_Variance': pc_variances[0] if len(pc_variances) > 0 else 0,
            'PC2_Variance': pc_variances[1] if len(pc_variances) > 1 else 0,
            'PC3_Variance': pc_variances[2] if len(pc_variances) > 2 else 0,
            'Min_Intensity': bin_data['GFP_fluorescence'].min(),
            'Max_Intensity': bin_data['GFP_fluorescence'].max()
        })
        
        print(f"  Bin {bin_num}: {len(bin_data):4d} variants, mean PCA distance: {mean_distance:.3f}")
    
    # Save diversity metrics
    diversity_df = pd.DataFrame(bin_diversity_metrics)
    diversity_df.to_csv(os.path.join(output_dir, "04_sequence_diversity_by_bin.csv"), index=False)
    
    print("Sequence diversity analysis completed (filtered data in full PCA space)")
else:
    print("PCA coordinates not available - skipping diversity analysis")
    diversity_df = None

# Step 7: User-configurable variant selection with distribution-based sampling
# Default values - can be modified by user
TOTAL_VARIANTS_DESIRED = 400 # User can change this
MIN_VARIANTS_PER_BIN = 2       # Minimum variants per bin (to ensure coverage)
SELECTION_MODE = "proportional"  # Options: "proportional", "equal", "hybrid"

# Note: Additional configurable parameters throughout the script:
# - KMER_SIZE: K-mer analysis (4=256 features, faster | 5=1024 features, more detailed but slower)
# - MIN_SEQUENCE_LENGTH & MAX_SEQUENCE_LENGTH: Sequence length range filter (e.g., 50-65bp)
# - N_CLUSTERS: Number of sequence clusters for diversity analysis (default: 8)
# - USE_LOG_SCALE_BINNING: Log vs linear intensity binning (default: True for log scale)

print(f"\n=== VARIANT SELECTION ===")
print(f"Target total variants: {TOTAL_VARIANTS_DESIRED}")
print(f"Selection mode: {SELECTION_MODE}")

# Calculate variants per bin based on selection mode
if SELECTION_MODE == "equal":
    # Original approach: equal variants per bin
    variants_per_bin = TOTAL_VARIANTS_DESIRED // n_bins
    remainder = TOTAL_VARIANTS_DESIRED % n_bins
    variants_per_bin_dict = {}
    for i in range(1, n_bins+1):
        variants_per_bin_dict[i] = variants_per_bin + (1 if remainder > 0 else 0)
        if remainder > 0:
            remainder -= 1
    
elif SELECTION_MODE == "proportional":
    # New approach: proportional to natural distribution
    variants_per_bin_dict = {}
    total_constructs_in_bins = sum([len(constructs[constructs['intensity_bin'] == i]) for i in range(1, n_bins+1)])
    
    for i in range(1, n_bins+1):
        bin_size = len(constructs[constructs['intensity_bin'] == i])
        if bin_size > 0:
            proportion = bin_size / total_constructs_in_bins
            calculated_variants = int(TOTAL_VARIANTS_DESIRED * proportion)
            # Ensure minimum representation
            variants_per_bin_dict[i] = max(MIN_VARIANTS_PER_BIN, calculated_variants)
        else:
            variants_per_bin_dict[i] = 0
    
elif SELECTION_MODE == "hybrid":
    # Hybrid approach: balance between equal and proportional
    variants_per_bin_dict = {}
    total_constructs_in_bins = sum([len(constructs[constructs['intensity_bin'] == i]) for i in range(1, n_bins+1)])
    equal_share = TOTAL_VARIANTS_DESIRED // n_bins
    
    for i in range(1, n_bins+1):
        bin_size = len(constructs[constructs['intensity_bin'] == i])
        if bin_size > 0:
            proportion = bin_size / total_constructs_in_bins
            proportional_share = int(TOTAL_VARIANTS_DESIRED * proportion)
            # Average of equal and proportional
            hybrid_share = int((equal_share + proportional_share) / 2)
            variants_per_bin_dict[i] = max(MIN_VARIANTS_PER_BIN, hybrid_share)
        else:
            variants_per_bin_dict[i] = 0

# Adjust to hit target total (may go slightly over due to minimums)
actual_total = sum(variants_per_bin_dict.values())
print(f"Calculated selection: {actual_total} variants")

# Show selection breakdown
print("Selection by bin:")
for i in range(1, n_bins+1):
    bin_size = len(constructs[constructs['intensity_bin'] == i])
    if bin_size > 0:
        selection_pct = (variants_per_bin_dict[i] / bin_size) * 100
        print(f"  Bin {i:2d}: {variants_per_bin_dict[i]:3d} from {bin_size:4d} available ({selection_pct:5.1f}%)")
    else:
        print(f"  Bin {i:2d}: {variants_per_bin_dict[i]:3d} from {bin_size:4d} available (N/A%)")

def select_diverse_variants_from_bin(bin_data, bin_kmer_matrix, n_select):
    """
    Select diverse variants from an intensity bin using k-mer distance analysis.
    
    This function implements a sophisticated selection strategy that:
    1. Filters out sequences containing BsaI restriction sites (GGTCTC/GAGACC)
    2. Ensures extreme intensity values (min/max) are always included
    3. Uses k-mer distance analysis to select maximally diverse sequences
    4. Falls back to random sampling for very large bins (>500 variants)
    
    The selection prioritizes functional diversity by choosing sequences that are
    maximally different from each other in their k-mer composition while maintaining
    representation across the full intensity range of the bin.
    
    Args:
        bin_data (pd.DataFrame): DataFrame containing sequences from one intensity bin
        bin_kmer_matrix (np.array): K-mer feature matrix for the sequences in bin_data
        n_select (int): Number of variants to select from this bin
    
    Returns:
        pd.DataFrame: Selected variants from the bin, or empty DataFrame if no valid sequences
    """
    # Filter out sequences containing >2 BsaI sites (GGTCTC and reverse complement GAGACC)
    print(f"    Checking for BsaI sites (GGTCTC/GAGACC) in {len(bin_data)} sequences...")
    print(f"    Filtering criteria: Remove sequences with >2 total BsaI sites")
    bsai_site_forward = "GGTCTC"
    bsai_site_reverse = "GAGACC"  # Reverse complement of GGTCTC
    
    # Count BsaI sites in each sequence (both orientations)
    def count_bsai_sites(sequence):
        if pd.isna(sequence):
            return 0
        seq_upper = sequence.upper()
        forward_count = seq_upper.count(bsai_site_forward)
        reverse_count = seq_upper.count(bsai_site_reverse)
        return forward_count + reverse_count
    
    bin_data['bsai_site_count'] = bin_data['trimmed_sequence'].apply(count_bsai_sites)
    
    # Filter out sequences with >2 BsaI sites
    bin_data_filtered = bin_data[bin_data['bsai_site_count'] <= 2].copy()
    
    if len(bin_data_filtered) < len(bin_data):
        removed_count = len(bin_data) - len(bin_data_filtered)
        print(f"    Removed {removed_count} sequences with >2 BsaI sites")
        
        # Show distribution of BsaI site counts
        site_counts = bin_data['bsai_site_count'].value_counts().sort_index()
        print(f"    BsaI site distribution:")
        for count, freq in site_counts.items():
            status = "✓ kept" if count <= 2 else "✗ removed"
            print(f"      {count} sites: {freq} sequences ({status})")
        
        # Update k-mer matrix to match filtered data
        if len(bin_data_filtered) > 0:
            # Get indices of filtered sequences in the original bin_data
            filtered_indices = [i for i, idx in enumerate(bin_data.index) if idx in bin_data_filtered.index]
            bin_kmer_matrix = bin_kmer_matrix[filtered_indices]
        else:
            print(f"    WARNING: No sequences remain after BsaI filtering!")
            return pd.DataFrame()  # Return empty dataframe
    else:
        print(f"    All sequences have ≤2 BsaI sites - all sequences are compatible")
    
    # Use filtered data for selection
    bin_data = bin_data_filtered
    
    if len(bin_data) <= n_select:
        return bin_data.copy()
    
    # Always ensure we include the minimum and maximum intensity values from each bin
    bin_sorted = bin_data.sort_values('GFP_fluorescence')
    
    if n_select <= 2:
        # For very small selections, prioritize extreme values
        if n_select == 1:
            return bin_sorted.iloc[[0]].copy()  # Take minimum
        else:  # n_select == 2
            return bin_sorted.iloc[[0, -1]].copy()  # Take min and max
    
    # For larger bins, use faster sampling methods but include extremes
    if len(bin_data) > 500:
        print(f"    Large bin ({len(bin_data)} variants) - using smart sampling with extremes")
        # Always include min and max
        extreme_indices = [0, len(bin_sorted)-1]
        selected_indices = set(extreme_indices)
        
        # Add random middle values
        remaining_needed = n_select - len(selected_indices)
        if remaining_needed > 0:
            middle_indices = list(range(1, len(bin_sorted)-1))
            if len(middle_indices) > 0:
                np.random.seed(42)
                additional_indices = np.random.choice(middle_indices, 
                                                    size=min(remaining_needed, len(middle_indices)), 
                                                    replace=False)
                selected_indices.update(additional_indices)
        
        return bin_sorted.iloc[list(selected_indices)].copy()
    
    try:
        # Calculate pairwise distances based on k-mer profiles
        print(f"    Computing distances for {len(bin_data)} variants...")
        distances = pairwise_distances(bin_kmer_matrix, metric='euclidean')
        
        # Greedy selection for maximum diversity, but ensure extremes are included
        selected_indices = []
        
        # Always start with extreme intensity values
        selected_indices.extend([0, len(bin_data)-1])  # Min and max from sorted data
        
        # Select remaining variants to maximize minimum distance to already selected
        for selection_round in range(n_select - 2):  # -2 because we already have 2 extremes
            if selection_round % 5 == 0 and selection_round > 0:
                print(f"      Selection round {selection_round+3}/{n_select}")  # +3 because we start with 2
            
            max_min_dist = -1
            best_idx = -1
            
            # Map back to original bin_data indices (since we're working with sorted data)
            original_indices = [bin_data.index[bin_sorted.index[i]] for i in selected_indices]
            bin_data_reset = bin_data.reset_index(drop=True)
            
            for i in range(len(bin_data)):
                if i in selected_indices:
                    continue
                
                # Calculate minimum distance to already selected variants
                min_dist = min([distances[i][j] for j in selected_indices])
                
                if min_dist > max_min_dist:
                    max_min_dist = min_dist
                    best_idx = i
            
            if best_idx != -1:
                selected_indices.append(best_idx)
        
        print(f"    Diversity selection with extremes completed")
        return bin_sorted.iloc[selected_indices].copy()
        
    except Exception as e:
        print(f"    Error in diversity selection: {e}")
        print(f"    Falling back to smart sampling with extremes")
        # Fallback: include extremes + random middle values
        extreme_indices = [0, len(bin_sorted)-1]
        selected_indices = set(extreme_indices)
        
        remaining_needed = n_select - len(selected_indices)
        if remaining_needed > 0:
            middle_indices = list(range(1, len(bin_sorted)-1))
            if len(middle_indices) > 0:
                np.random.seed(42)
                additional_indices = np.random.choice(middle_indices, 
                                                    size=min(remaining_needed, len(middle_indices)), 
                                                    replace=False)
                selected_indices.update(additional_indices)
        
        return bin_sorted.iloc[list(selected_indices)].copy()

# Select variants from each bin
selected_variants = []
selection_summary = []

print(f"\n=== SELECTING VARIANTS FROM BINS ===")
for bin_num in range(1, n_bins+1):
    print(f"\nProcessing Bin {bin_num}/{n_bins}...")
    bin_mask = constructs['intensity_bin'] == bin_num
    bin_data = constructs[bin_mask].copy()
    
    # Get target number for this bin
    n_select = variants_per_bin_dict[bin_num]
    
    if len(bin_data) == 0:
        print(f"Bin {bin_num}: No variants available (target: {n_select})")
        continue
    
    if n_select == 0:
        print(f"Bin {bin_num}: Skipping (target: 0 variants)")
        continue
    
    # Get k-mer matrix for this bin from full library matrix using PCA indices
    bin_pca_indices = bin_data['pca_index'].tolist()
    bin_kmer_matrix = kmat_full_norm[bin_pca_indices]
    
    # Don't select more than available
    n_select = min(n_select, len(bin_data))
    
    # Select diverse variants
    selected_bin_variants = select_diverse_variants_from_bin(bin_data, bin_kmer_matrix, n_select)
    selected_variants.append(selected_bin_variants)
    
    selection_pct = (n_select / len(bin_data)) * 100
    print(f"Bin {bin_num}: Selected {len(selected_bin_variants)} from {len(bin_data)} available variants ({selection_pct:.1f}%)")
    
    selection_summary.append({
        'Bin': bin_num,
        'Available': len(bin_data),
        'Selected': len(selected_bin_variants),
        'Selection_Percentage': selection_pct,
        'Min_Intensity': bin_data['GFP_fluorescence'].min(),
        'Max_Intensity': bin_data['GFP_fluorescence'].max()
    })

# Combine all selected variants
if selected_variants:
    final_selection = pd.concat(selected_variants, ignore_index=True)
    print(f"\nTotal variants selected: {len(final_selection)}")
    
    # Verify no BsaI sites in final selection (check both orientations)
    bsai_forward_check = final_selection['trimmed_sequence'].str.upper().str.contains("GGTCTC", na=False).sum()
    bsai_reverse_check = final_selection['trimmed_sequence'].str.upper().str.contains("GAGACC", na=False).sum()
    total_bsai_check = bsai_forward_check + bsai_reverse_check
    
    if total_bsai_check > 0:
        print(f"WARNING: {total_bsai_check} selected sequences still contain BsaI sites!")
        if bsai_forward_check > 0:
            print(f"  - {bsai_forward_check} contain GGTCTC (forward)")
        if bsai_reverse_check > 0:
            print(f"  - {bsai_reverse_check} contain GAGACC (reverse complement)")
    else:
        print("✓ Confirmed: No BsaI sites (GGTCTC/GAGACC) in selected variants")
else:
    final_selection = pd.DataFrame()
    print("No variants selected!")

# Save selection summary
selection_summary_df = pd.DataFrame(selection_summary)
selection_summary_df.to_csv(os.path.join(output_dir, "04_selection_summary.csv"), index=False)

# Step 8: Create output with required columns
if len(final_selection) > 0:
    # Sort final selection by increasing protein fluorescence
    final_selection = final_selection.sort_values('GFP_fluorescence', ascending=True).reset_index(drop=True)
    
    # Create output dataframe with required columns
    output_columns = ['target', 'trimmed_sequence', 'GFP_fluorescence', 'sequence_length']
    
    # Check if target column exists, if not use index or create
    if 'target' not in final_selection.columns:
        final_selection['target'] = 'construct_' + final_selection.index.astype(str)
    
    output_data = final_selection[output_columns].copy()
    output_data.columns = ['Target_Name', 'Sequence', 'Protein_Fluorescence', 'Sequence_Length']
    
    # Save final selection (already sorted by increasing protein fluorescence)
    output_data.to_csv(os.path.join(output_dir, "05_selected_variants.csv"), index=False)
    
    print(f"\nFinal selection saved with {len(output_data)} variants")
    print(f"Variants sorted by increasing protein fluorescence")
    print(f"Intensity range: {output_data['Protein_Fluorescence'].min():.1f} - {output_data['Protein_Fluorescence'].max():.1f}")
    print(f"Sequence length range: {output_data['Sequence_Length'].min()} - {output_data['Sequence_Length'].max()} bp")
else:
    print("No variants to save!")

# Step 9: Add BsaI sites and primer annealing sites to selected variants
# User-configurable BsaI sites for 5' and 3' ends
# These sites will be added to the selected variants for cloning purposes
BSAI_SITE_5_PRIME = "GGTCTCTgtac"  # User can change this
BSAI_SITE_3_PRIME = "ggtgTGAGACC"  # User can change this

# User-configurable primer annealing sites for PCR amplification
# These sites will be added outside the BsaI sites for primer binding
# Structure: [PRIMER_5] + [BSAI_5] + [SEQUENCE] + [BSAI_3] + [PRIMER_3]
PRIMER_SITE_5_PRIME = "GTAAAACGACGGCCAGT"     # M13 forward primer site
PRIMER_SITE_3_PRIME = "CAGGAAACAGCTATGAC"     # M13 reverse primer site


def add_bsai_sites_to_sequences(df, site_5_prime, site_3_prime):
    """
    Add user-defined BsaI sites to the 5' and 3' ends of sequences with junction analysis.
    
    This function adds BsaI sites to sequences and performs comprehensive junction analysis
    to identify and filter out sequences that would create unwanted restriction sites at
    the 5' and 3' junctions after concatenation.
    
    Args:
        df (pd.DataFrame): DataFrame containing sequences in 'Sequence' column
        site_5_prime (str): BsaI site to add to 5' end
        site_3_prime (str): BsaI site to add to 3' end
    
    Returns:
        pd.DataFrame: DataFrame with sequences that pass junction analysis, 
                     includes columns for sequences with BsaI sites and junction check results
    """
    # Define common restriction sites to check for at junctions
    # These are sites that could interfere with cloning
    problematic_sites = {
        'BsaI_forward': 'GGTCTC',
        'BsaI_reverse': 'GAGACC', 
        # 'EcoRI': 'GAATTC',
        # 'BamHI': 'GGATCC',
        # 'HindIII': 'AAGCTT',
        # 'XhoI': 'CTCGAG',
        # 'SalI': 'GTCGAC',
        # 'NotI': 'GCGGCCGC',
        # 'SpeI': 'ACTAGT',
        # 'XbaI': 'TCTAGA'
    }
    
    print(f"    Analyzing junctions for {len(df)} sequences...")
    print(f"    5' site: {site_5_prime}")
    print(f"    3' site: {site_3_prime}")
    print(f"    Checking for {len(problematic_sites)} restriction sites at junctions")
    
    # Create a copy of the input DataFrame
    result_df = df.copy()
    
    # Add initial concatenated sequences (strip whitespace to ensure clean concatenation)
    result_df['Sequence_with_BsaI'] = site_5_prime.strip() + result_df['Sequence'].str.strip() + site_3_prime.strip()
    result_df['Sequence_Length_with_BsaI'] = result_df['Sequence_with_BsaI'].str.len()
    result_df['BsaI_5_Prime'] = site_5_prime
    result_df['BsaI_3_Prime'] = site_3_prime
    
    # Initialize junction analysis columns
    result_df['Junction_5_Prime_Clean'] = True
    result_df['Junction_3_Prime_Clean'] = True
    result_df['Problematic_Sites_Found'] = ''
    
    # Analyze junctions for each sequence
    junction_length = 12  # Check 6 bp on each side of junction
    problematic_sequences = []
    
    for idx, row in result_df.iterrows():
        sequence = row['Sequence']
        problems_found = []
        
        # Check 5' junction: only flag truly problematic sites
        if len(sequence) >= 6:
            junction_5 = site_5_prime + sequence[:6]
            for site_name, site_seq in problematic_sites.items():
                if site_seq in junction_5:
                    # Skip BsaI sites that are part of our intentionally added sites
                    if site_name in ['BsaI_forward', 'BsaI_reverse'] and site_seq in site_5_prime:
                        continue  # This is intentional, don't flag it
                    # Only flag if it's a truly accidental site
                    result_df.at[idx, 'Junction_5_Prime_Clean'] = False
                    problems_found.append(f"5'_{site_name}")
        
        # Check 3' junction: only flag truly problematic sites  
        if len(sequence) >= 6:
            junction_3 = sequence[-6:] + site_3_prime
            for site_name, site_seq in problematic_sites.items():
                if site_seq in junction_3:
                    # Skip BsaI sites that are part of our intentionally added sites
                    if site_name in ['BsaI_forward', 'BsaI_reverse'] and site_seq in site_3_prime:
                        continue  # This is intentional, don't flag it
                    # Only flag if it's a truly accidental site
                    result_df.at[idx, 'Junction_3_Prime_Clean'] = False
                    problems_found.append(f"3'_{site_name}")
        
        # Record problems found
        if problems_found:
            result_df.at[idx, 'Problematic_Sites_Found'] = ';'.join(problems_found)
            problematic_sequences.append(idx)
    
    # Filter to keep only clean sequences
    clean_sequences = result_df[
        (result_df['Junction_5_Prime_Clean'] == True) & 
        (result_df['Junction_3_Prime_Clean'] == True)
    ].copy()
    
    # Report filtering statistics
    total_input = len(df)
    total_problematic = len(problematic_sequences)
    total_clean = len(clean_sequences)
    
    print(f"    Junction Analysis Results:")
    print(f"      Input sequences: {total_input}")
    print(f"      Problematic junctions: {total_problematic}")
    print(f"      Clean sequences retained: {total_clean}")
    
    if total_problematic > 0:
        print(f"      Filtering efficiency: {(total_clean/total_input)*100:.1f}%")
        
        # Show breakdown of problems found
        problem_counts = {}
        for idx in problematic_sequences:
            problems = result_df.at[idx, 'Problematic_Sites_Found'].split(';')
            for problem in problems:
                problem_counts[problem] = problem_counts.get(problem, 0) + 1
        
        print(f"      Problem breakdown:")
        for problem, count in sorted(problem_counts.items()):
            print(f"        {problem}: {count} sequences")
    else:
        print(f"      ✓ All sequences passed junction analysis!")
    
    return clean_sequences

def add_primer_sites_to_sequences(df, primer_5_prime, primer_3_prime):
    """
    Add primer annealing sites to sequences that already have BsaI sites.
    
    This function adds primer binding sites outside the BsaI sites to enable
    PCR amplification of the constructs. The final structure is:
    [PRIMER_5] + [BSAI_5] + [SEQUENCE] + [BSAI_3] + [PRIMER_3]
    
    Args:
        df (pd.DataFrame): DataFrame containing sequences with BsaI sites
        primer_5_prime (str): 5' primer annealing site
        primer_3_prime (str): 3' primer annealing site
    
    Returns:
        pd.DataFrame: DataFrame with primer sites added, ready for ordering
    """
    print(f"    Adding primer annealing sites to {len(df)} sequences...")
    print(f"    5' primer site: {primer_5_prime}")
    print(f"    3' primer site: {primer_3_prime}")
    
    # Check if dataframe is empty
    if len(df) == 0:
        print("    WARNING: No sequences provided to add primer sites!")
        print("    This likely means all sequences were filtered out during BsaI junction analysis.")
        return df.copy()  # Return empty dataframe with same structure
    
    # Create a copy of the input DataFrame
    result_df = df.copy()
    
    # Add primer sites to sequences with BsaI sites (strip whitespace to ensure clean concatenation)
    result_df['Final_Sequence'] = primer_5_prime.strip() + result_df['Sequence_with_BsaI'].str.strip() + primer_3_prime.strip()
    result_df['Final_Sequence_Length'] = result_df['Final_Sequence'].str.len()
    result_df['Primer_5_Prime'] = primer_5_prime
    result_df['Primer_3_Prime'] = primer_3_prime
    
    # Calculate total additions (now safe since we checked df is not empty)
    bsai_addition = len(result_df['BsaI_5_Prime'].iloc[0]) + len(result_df['BsaI_3_Prime'].iloc[0])
    primer_addition = len(primer_5_prime) + len(primer_3_prime)
    total_addition = bsai_addition + primer_addition
    
    print(f"    Primer sites added successfully!")
    print(f"    BsaI sites addition: {bsai_addition} bp")
    print(f"    Primer sites addition: {primer_addition} bp") 
    print(f"    Total addition per sequence: {total_addition} bp")
    print(f"    Final sequence length range: {result_df['Final_Sequence_Length'].min()} - {result_df['Final_Sequence_Length'].max()} bp")
    
    return result_df

if len(final_selection) > 0:
    print(f"\n=== ADDING BsaI SITES TO SELECTED VARIANTS ===")
    print(f"5' BsaI site: {BSAI_SITE_5_PRIME}")
    print(f"3' BsaI site: {BSAI_SITE_3_PRIME}")
    
    # Add BsaI sites to the output data
    output_data_with_sites = add_bsai_sites_to_sequences(output_data, BSAI_SITE_5_PRIME, BSAI_SITE_3_PRIME)
    
    # Save variants with BsaI sites
    output_data_with_sites.to_csv(os.path.join(output_dir, "selected_variants_with_sites.csv"), index=False)
    
    print(f"Selected variants with BsaI sites saved: {len(output_data_with_sites)} variants")
    print(f"Original sequence length range: {output_data_with_sites['Sequence_Length'].min()} - {output_data_with_sites['Sequence_Length'].max()} bp")
    print(f"With BsaI sites length range: {output_data_with_sites['Sequence_Length_with_BsaI'].min()} - {output_data_with_sites['Sequence_Length_with_BsaI'].max()} bp")
    print(f"Length increase per sequence: {len(BSAI_SITE_5_PRIME) + len(BSAI_SITE_3_PRIME)} bp")
    
    # Step 10: Add primer annealing sites for PCR amplification
    if len(output_data_with_sites) > 0:
        print(f"\n=== ADDING PRIMER ANNEALING SITES ===")
        print(f"5' primer site: {PRIMER_SITE_5_PRIME}")
        print(f"3' primer site: {PRIMER_SITE_3_PRIME}")
        
        # Add primer sites to sequences with BsaI sites
        ready_to_order_data = add_primer_sites_to_sequences(output_data_with_sites, PRIMER_SITE_5_PRIME, PRIMER_SITE_3_PRIME)
        
        if len(ready_to_order_data) > 0:
            # Create final output dataframe with key columns for ordering
            final_output_columns = [
                'Target_Name', 'Final_Sequence', 'Final_Sequence_Length',
                'Protein_Fluorescence', 'Sequence_Length',
                'BsaI_5_Prime', 'BsaI_3_Prime', 
                'Primer_5_Prime', 'Primer_3_Prime'
            ]
            ready_to_order_final = ready_to_order_data[final_output_columns].copy()
            
            # Save ready-to-order variants
            ready_to_order_final.to_csv(os.path.join(output_dir, "Ready_to_order_variants.csv"), index=False)
            
            print(f"\n=== READY-TO-ORDER VARIANTS SAVED ===")
            print(f"File: Ready_to_order_variants.csv")
            print(f"Total variants ready for ordering: {len(ready_to_order_final)}")
            print(f"Final sequence structure: [PRIMER_5] + [BSAI_5] + [SEQUENCE] + [BSAI_3] + [PRIMER_3]")
            print(f"Ready-to-order sequence length range: {ready_to_order_final['Final_Sequence_Length'].min()} - {ready_to_order_final['Final_Sequence_Length'].max()} bp")
            
            # =============================================================================
            # READY-TO-ORDER SEQUENCE LENGTH ANALYSIS
            # =============================================================================
            print(f"\n" + "="*60)
            print(f"READY-TO-ORDER SEQUENCE LENGTH ANALYSIS")
            print(f"="*60)
            
            # Calculate detailed length statistics
            final_lengths = ready_to_order_final['Final_Sequence_Length']
            original_lengths = ready_to_order_final['Sequence_Length']
            
            # Component length contributions
            primer_5_len = len(PRIMER_SITE_5_PRIME.strip())
            bsai_5_len = len(BSAI_SITE_5_PRIME.strip())
            bsai_3_len = len(BSAI_SITE_3_PRIME.strip())
            primer_3_len = len(PRIMER_SITE_3_PRIME.strip())
            total_added = primer_5_len + bsai_5_len + bsai_3_len + primer_3_len
            
            print(f"LENGTH BREAKDOWN:")
            print(f"  5' Primer site:        {primer_5_len:3d} bp")
            print(f"  5' BsaI site:          {bsai_5_len:3d} bp")
            print(f"  Original sequence:     {original_lengths.min():3d}-{original_lengths.max():3d} bp (range)")
            print(f"  3' BsaI site:          {bsai_3_len:3d} bp")
            print(f"  3' Primer site:        {primer_3_len:3d} bp")
            print(f"  Total additions:       {total_added:3d} bp")
            
            print(f"\nFINAL LENGTH STATISTICS:")
            print(f"  Mean length:           {final_lengths.mean():.1f} bp")
            print(f"  Median length:         {final_lengths.median():.0f} bp")
            print(f"  Standard deviation:    {final_lengths.std():.1f} bp")
            print(f"  Range:                 {final_lengths.min()}-{final_lengths.max()} bp")
            
            # Length distribution
            length_counts = final_lengths.value_counts().sort_index()
            print(f"\nLENGTH DISTRIBUTION:")
            for length, count in length_counts.head(10).items():  # Show top 10 most common lengths
                pct = (count / len(final_lengths)) * 100
                print(f"  {length:3d} bp: {count:3d} sequences ({pct:4.1f}%)")
            if len(length_counts) > 10:
                print(f"  ... and {len(length_counts) - 10} other lengths")
            
            # Create length distribution plot
            print(f"\nCreating ready-to-order sequence length distribution plot...")
            
            plt.figure(figsize=(14, 6))
            
            # Plot 1: Bar plot of exact length counts with labels (all consecutive lengths)
            plt.subplot(1, 2, 1)
            length_counts_plot = final_lengths.value_counts().sort_index()
            
            # Create complete range of lengths from min to max
            max_length = final_lengths.max()
            min_length = final_lengths.min()
            all_lengths = list(range(min_length, max_length + 1))
            
            # Get counts for all lengths (0 for missing lengths)
            all_counts = [length_counts_plot.get(length, 0) for length in all_lengths]
            
            # Create colors for bars (special colors for min/max, default for others)
            bar_colors = []
            for length in all_lengths:
                if length == max_length:
                    bar_colors.append('darkred')
                elif length == min_length:
                    bar_colors.append('darkblue')
                else:
                    bar_colors.append('darkgreen')
            
            # Plot all bars
            bars = plt.bar(all_lengths, all_counts, alpha=0.8, color=bar_colors, edgecolor='black')
            
            # Add value labels on each bar (only for bars with count > 0)
            max_count = max(all_counts) if all_counts else 1
            for length, count in zip(all_lengths, all_counts):
                if count > 0:  # Only label bars that have sequences
                    plt.text(length, count + max_count * 0.01, str(count), 
                            ha='center', va='bottom', fontweight='bold', fontsize=9)
            
            # Create legend for min/max highlighting
            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor='darkgreen', label='Standard lengths'),
                             Patch(facecolor='darkblue', label=f'Min length: {min_length} bp'),
                             Patch(facecolor='darkred', label=f'Max length: {max_length} bp')]
            
            plt.title('Exact Length Value Counts')
            plt.xlabel('Final Sequence Length (bp)')
            plt.ylabel('Number of Sequences')
            plt.xticks(all_lengths, rotation=45)  # Show all consecutive length values
            plt.legend(handles=legend_elements)
            plt.grid(True, alpha=0.3, axis='y')
            
            # Plot 2: Component contribution breakdown
            plt.subplot(1, 2, 2)
            components = ['5\' Primer', '5\' BsaI', 'Original\nSequence', '3\' BsaI', '3\' Primer']
            component_lengths = [primer_5_len, bsai_5_len, original_lengths.mean(), bsai_3_len, primer_3_len]
            colors = ['lightblue', 'lightgreen', 'gold', 'lightcoral', 'plum']
            
            # Calculate error bars for original sequence (min to max range)
            original_min = original_lengths.min()
            original_max = original_lengths.max()
            original_mean = original_lengths.mean()
            
            # Create error bar data (only for original sequence)
            error_lower = [0, 0, original_mean - original_min, 0, 0]  # Distance from mean to min
            error_upper = [0, 0, original_max - original_mean, 0, 0]  # Distance from mean to max
            
            bars = plt.bar(components, component_lengths, color=colors, edgecolor='black', alpha=0.8,
                          yerr=[error_lower, error_upper], capsize=5, error_kw={'linewidth': 2, 'capthick': 2})
            plt.title('Sequence Component Length Breakdown')
            plt.ylabel('Length (bp)')
            plt.xticks(rotation=45)
            
            # Add value labels on bars
            for i, (bar, length) in enumerate(zip(bars, component_lengths)):
                if i == 2:  # Original sequence bar - show range
                    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + (original_max - original_mean) + 2, 
                            f'{original_min:.0f}-{original_max:.0f}', ha='center', va='bottom', fontweight='bold')
                    # Also show mean value
                    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() - 5, 
                            f'μ={original_mean:.0f}', ha='center', va='top', fontweight='bold', fontsize=8, color='black')
                else:
                    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, 
                            f'{length:.0f}', ha='center', va='bottom', fontweight='bold')
            
            plt.tight_layout()
            
            # Save the plot
            length_plot_filename = os.path.join(output_dir, "Ready_to_order_length_analysis.png")
            plt.savefig(length_plot_filename, dpi=300, bbox_inches='tight')
            plt.show()
            
            print(f"Ready-to-order length analysis plot saved: Ready_to_order_length_analysis.png")
            print(f"="*60)
        else:
            print(f"\n⚠️  NO SEQUENCES READY FOR ORDERING!")
            print(f"All sequences were filtered out during processing.")
    else:
        print(f"\n⚠️  NO SEQUENCES SURVIVED BsaI JUNCTION ANALYSIS!")
        print(f"Reason: All selected variants create problematic restriction sites when BsaI sites are added.")
        print(f"Solutions:")
        print(f"  1. Try different BsaI site sequences (lines 1314-1315)")
        print(f"  2. Relax restriction site checking in junction analysis")
        print(f"  3. Select more variants to increase chances of finding compatible ones")
    
else:
    print("No variants available to add BsaI sites!")

# Save k-mer analysis results (full library)
if 'kmat_full' in locals():
    np.savetxt(os.path.join(output_dir, "04_kmer_matrix.csv"), kmat_full, delimiter=',')
    np.savetxt(os.path.join(output_dir, "04_kmer_matrix_normalized.csv"), kmat_full_norm, delimiter=',')
    print("\nK-mer matrices saved for future analysis (full library)")
else:
    print("\nK-mer matrices not available to save")

# Save checkpoint before selection
print("Saving checkpoint before variant selection...")
constructs.to_csv(os.path.join(output_dir, "checkpoint_before_selection.csv"), index=False)

# Monitor memory usage
import psutil
import os
process = psutil.Process(os.getpid())
memory_usage = process.memory_info().rss / 1024 / 1024  # MB
print(f"Current memory usage: {memory_usage:.1f} MB")

# %% [markdown]
# ## Visualization and Summary

# %%
from matplotlib.backends.backend_pdf import PdfPages

# matplotlib backend already set at top of script

# Create comprehensive visualizations and PDF report
print("\n=== CREATING VISUALIZATIONS AND PDF REPORT ===")
print("This may take a few minutes for large datasets...")

if len(final_selection) > 0:
    # Create PDF report with all plots
    pdf_filename = os.path.join(output_dir, "complete_analysis_report.pdf")
    print(f"Creating PDF report: {pdf_filename}")
    
    with PdfPages(pdf_filename) as pdf:
        # Page 1: Full Library vs Filtered Library PCA Analysis
        print("Creating Page 1: Full Library vs Filtered Library PCA Analysis...")
        fig = plt.figure(figsize=(18, 12))
        
        # Title page info
        fig.suptitle(f'Combinatorial Library Analysis Report\nPage 1: Full Library vs Filtered Library PCA Analysis\nGenerated: {timestamp}', 
                     fontsize=16, y=0.95)
        
        # We need both full library and filtered library data
        # Get full library data (before filtering)
        if 'full_library_constructs' in locals() and 'pca_features_full' in locals():
            # Add PCA coordinates to full library if available
            if len(full_library_constructs) == pca_features_full.shape[0]:
                for i in range(min(10, pca_features_full.shape[1])):
                    full_library_constructs[f'PC{i+1}'] = pca_features_full[:, i]
                full_sample = full_library_constructs
            else:
                full_sample = constructs  # Fallback to filtered data
        else:
            full_sample = constructs  # Fallback to filtered data
        
        filtered_sample = constructs  # Filtered library data
        
        has_pc2_full = 'PC2' in full_sample.columns
        has_pc2_filtered = 'PC2' in filtered_sample.columns
        has_clusters = 'sequence_cluster' in filtered_sample.columns
        
        # Plot 1 (2,3,1): Full Library PC1 vs PC2 colored by intensity
        plt.subplot(2, 3, 1)
        if has_pc2_full:
            scatter = plt.scatter(full_sample['PC1'], full_sample['PC2'], 
                                c=full_sample['GFP_fluorescence'], cmap='viridis', 
                                alpha=0.6, s=8, marker='o')
            plt.colorbar(scatter, label='GFP Fluorescence')
            plt.xlabel('PC1')
            plt.ylabel('PC2')
            plt.title(f'Full Library: PC1 vs PC2 by Intensity\n(n={len(full_sample):,})')
        else:
            plt.text(0.5, 0.5, 'PC2 not available', ha='center', va='center', transform=plt.gca().transAxes)
            plt.title('Full Library: PC1 vs PC2 (PC2 not available)')
        
        
        # Plot 2 (2,3,2): Full Library PC1 vs PC2 colored by cluster
        plt.subplot(2, 3, 2)
        if has_pc2_full and has_clusters and 'sequence_cluster' in full_sample.columns:
            cluster_colors = plt.cm.Set3(np.linspace(0, 1, N_CLUSTERS))
            for cluster_id in range(N_CLUSTERS):
                cluster_data = full_sample[full_sample['sequence_cluster'] == cluster_id]
                if len(cluster_data) > 0:
                    plt.scatter(cluster_data['PC1'], cluster_data['PC2'], 
                              color=cluster_colors[cluster_id], alpha=0.7, s=8, 
                              label=f'Cluster {cluster_id}')
            plt.xlabel('PC1')
            plt.ylabel('PC2')
            plt.title(f'Full Library: PC1 vs PC2 by Cluster\n(n={len(full_sample):,})')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            plt.text(0.5, 0.5, 'PC2 or clusters not available', ha='center', va='center', transform=plt.gca().transAxes)
            plt.title('Full Library: PC1 vs PC2 by Cluster (not available)')
        
        # Plot 3 (2,3,3): Full Library PC1 vs Intensity (colored by cluster)
        plt.subplot(2, 3, 3)
        if 'PC1' in full_sample.columns:
            if has_clusters and 'sequence_cluster' in full_sample.columns:
                cluster_colors = plt.cm.Set3(np.linspace(0, 1, N_CLUSTERS))
                for cluster_id in range(N_CLUSTERS):
                    cluster_data = full_sample[full_sample['sequence_cluster'] == cluster_id]
                    if len(cluster_data) > 0:
                        plt.scatter(cluster_data['PC1'], cluster_data['GFP_fluorescence'], 
                                  color=cluster_colors[cluster_id], alpha=0.7, s=8, 
                                  label=f'Cluster {cluster_id}')
                plt.xlabel('PC1')
                plt.ylabel('GFP Fluorescence')
                plt.title(f'Full Library: PC1 vs Intensity by Cluster\n(n={len(full_sample):,})')
                plt.yscale('log')
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            else:
                plt.scatter(full_sample['PC1'], full_sample['GFP_fluorescence'], 
                           alpha=0.6, s=8, c='blue')
                plt.xlabel('PC1')
                plt.ylabel('GFP Fluorescence')
                plt.title(f'Full Library: PC1 vs Intensity\n(n={len(full_sample):,})')
                plt.yscale('log')
        else:
            plt.text(0.5, 0.5, 'PC1 not available', ha='center', va='center', transform=plt.gca().transAxes)
            plt.title('Full Library: PC1 vs Intensity (not available)')
        
        # Plot 4 (2,3,4): Filtered Library PC1 vs PC2 colored by intensity
        plt.subplot(2, 3, 4)
        if has_pc2_filtered:
            scatter = plt.scatter(filtered_sample['PC1'], filtered_sample['PC2'], 
                                c=filtered_sample['GFP_fluorescence'], cmap='viridis', 
                                alpha=0.6, s=8, marker='o')
            plt.colorbar(scatter, label='GFP Fluorescence')
            plt.xlabel('PC1')
            plt.ylabel('PC2')
            plt.title(f'Filtered Library: PC1 vs PC2 by Intensity\n(n={len(filtered_sample):,})')
        else:
            plt.text(0.5, 0.5, 'PC2 not available', ha='center', va='center', transform=plt.gca().transAxes)
            plt.title('Filtered Library: PC1 vs PC2 (PC2 not available)')
        
        # Plot 5 (2,3,5): Filtered Library PC1 vs PC2 colored by cluster
        plt.subplot(2, 3, 5)
        if has_pc2_filtered and has_clusters:
            cluster_colors = plt.cm.Set3(np.linspace(0, 1, N_CLUSTERS))
            for cluster_id in range(N_CLUSTERS):
                cluster_data = filtered_sample[filtered_sample['sequence_cluster'] == cluster_id]
                if len(cluster_data) > 0:
                    plt.scatter(cluster_data['PC1'], cluster_data['PC2'], 
                              color=cluster_colors[cluster_id], alpha=0.7, s=8, 
                              label=f'Cluster {cluster_id}')
            plt.xlabel('PC1')
            plt.ylabel('PC2')
            plt.title(f'Filtered Library: PC1 vs PC2 by Cluster\n(n={len(filtered_sample):,})')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            plt.text(0.5, 0.5, 'PC2 or clusters not available', ha='center', va='center', transform=plt.gca().transAxes)
            plt.title('Filtered Library: PC1 vs PC2 by Cluster (not available)')
        
        # Plot 6 (2,3,6): NEW - Filtered Library PC1 vs Intensity (colored by cluster)
        plt.subplot(2, 3, 6)
        if 'PC1' in filtered_sample.columns:
            if has_clusters:
                cluster_colors = plt.cm.Set3(np.linspace(0, 1, N_CLUSTERS))
                for cluster_id in range(N_CLUSTERS):
                    cluster_data = filtered_sample[filtered_sample['sequence_cluster'] == cluster_id]
                    if len(cluster_data) > 0:
                        plt.scatter(cluster_data['PC1'], cluster_data['GFP_fluorescence'], 
                                  color=cluster_colors[cluster_id], alpha=0.7, s=8, 
                                  label=f'Cluster {cluster_id}')
                plt.xlabel('PC1')
                plt.ylabel('GFP Fluorescence')
                plt.title(f'Filtered Library: PC1 vs Intensity by Cluster\n(n={len(filtered_sample):,})')
                plt.yscale('log')
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            else:
                plt.scatter(filtered_sample['PC1'], filtered_sample['GFP_fluorescence'], 
                           alpha=0.6, s=8, c='green')
                plt.xlabel('PC1')
                plt.ylabel('GFP Fluorescence')
                plt.title(f'Filtered Library: PC1 vs Intensity\n(n={len(filtered_sample):,})')
                plt.yscale('log')
        else:
            plt.text(0.5, 0.5, 'PC1 not available', ha='center', va='center', transform=plt.gca().transAxes)
            plt.title('Filtered Library: PC1 vs Intensity (not available)')
        

        
        # Adjust spacing for better layout
        plt.tight_layout(rect=[0, 0.03, 1, 0.93])
        print("Saving Page 1 to PDF...")
        pdf.savefig(fig, bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, "06_selection_analysis.png"), dpi=300, bbox_inches='tight')
        plt.close(fig)  # Close figure to free memory
        print("Page 1 completed.")
        
        # Page 2: Cluster Family Analysis
        print("Creating Page 2: Cluster Family Analysis...")
        fig = plt.figure(figsize=(16, 12))
        fig.suptitle('Cluster Family Analysis (Full Library)', fontsize=16, y=0.95)
        
        # Create a text-based table showing cluster composition
        plt.subplot(1, 1, 1)
        plt.axis('off')  # Turn off axes for text display
        
        # Build the analysis text
        analysis_text = []
        if 'cluster_families_analysis' in locals() and len(cluster_families_analysis) > 0:
            analysis_text.append("CLUSTER FAMILY ANALYSIS (Full Library)")
            analysis_text.append("")
            
            for family_data in cluster_families_analysis:
                cluster_id = family_data['Cluster']
                cluster_size = family_data['Count']
                promoter_families = family_data['Promoter_Families_Detail']
                rbs_families = family_data['RBS_Families_Detail']
                
                analysis_text.append(f"Cluster {cluster_id} (n={cluster_size}):")
                analysis_text.append(f"Promoters: {promoter_families}")
                analysis_text.append(f"RBS: {rbs_families}")
                analysis_text.append("")
        else:
            analysis_text.append("CLUSTER FAMILY ANALYSIS")
            analysis_text.append("")
            analysis_text.append("Cluster analysis data not available.")
            analysis_text.append("Please ensure the script completed the clustering analysis.")
        
        # Display the text with proper formatting
        full_text = '\n'.join(analysis_text)
        plt.text(0.05, 0.95, full_text, transform=plt.gca().transAxes, 
                fontsize=10, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))
        
        # Save Page 2
        print("Saving Page 2 to PDF...")
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)  # Close figure to free memory
        print("Page 2 completed.")
        
        # Page 3: Selected Variants Analysis
        print("Creating Page 3: Selected Variants Analysis...")
        fig = plt.figure(figsize=(18, 12))
        fig.suptitle('Selected Variants Analysis', fontsize=16, y=0.95)
        
        # Page 3: Selected Variants Analysis (2x3 layout)
        # Plot 1 (2,3,1): Selected Variants Colored by Intensity
        plt.subplot(2, 3, 1)
        if has_pc2_filtered:
            scatter = plt.scatter(final_selection['PC1'], final_selection['PC2'], 
                                c=final_selection['GFP_fluorescence'], cmap='viridis', 
                                alpha=0.8, s=20, marker='o')
            plt.colorbar(scatter, label='GFP Fluorescence')
            plt.xlabel('PC1')
            plt.ylabel('PC2')
            plt.title(f'Selected Variants: PC1 vs PC2 by Intensity\n(n={len(final_selection):,})')
        else:
            plt.scatter(final_selection['PC1'], np.zeros(len(final_selection)), 
                       c=final_selection['GFP_fluorescence'], cmap='viridis', 
                       alpha=0.8, s=20, marker='o')
            plt.colorbar(label='GFP Fluorescence')
            plt.xlabel('PC1')
            plt.ylabel('PC2 (not available)')
            plt.title(f'Selected Variants: PC1 by Intensity\n(n={len(final_selection):,})')
        
        # Plot 2 (2,3,2): Selected Variants Colored by cluster
        plt.subplot(2, 3, 2)
        if has_clusters:
            cluster_colors = plt.cm.Set3(np.linspace(0, 1, N_CLUSTERS))
            for cluster_id in range(N_CLUSTERS):
                cluster_data = final_selection[final_selection['sequence_cluster'] == cluster_id]
                if len(cluster_data) > 0:
                    if has_pc2_filtered:
                        plt.scatter(cluster_data['PC1'], cluster_data['PC2'], 
                                  color=cluster_colors[cluster_id], alpha=0.8, s=20, 
                                  label=f'Cluster {cluster_id}')
                    else:
                        plt.scatter(cluster_data['PC1'], np.zeros(len(cluster_data)), 
                                  color=cluster_colors[cluster_id], alpha=0.8, s=20, 
                                  label=f'Cluster {cluster_id}')
            plt.xlabel('PC1')
            plt.ylabel('PC2' if has_pc2_filtered else 'PC2 (not available)')
            plt.title(f'Selected Variants: PC1 vs PC2 by Cluster\n(n={len(final_selection):,})')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            plt.text(0.5, 0.5, 'Clusters not available', ha='center', va='center', transform=plt.gca().transAxes)
            plt.title('Selected Variants by Cluster (not available)')
        
        # Plot 3 (2,3,3): PC1 vs Intensity (filtered Clusters & Selected Variants)
        plt.subplot(2, 3, 3)
        # Plot filtered library variants colored by cluster
        if has_clusters:
            cluster_colors = plt.cm.Set3(np.linspace(0, 1, N_CLUSTERS))
            for cluster_id in range(N_CLUSTERS):
                cluster_data = filtered_sample[filtered_sample['sequence_cluster'] == cluster_id]
                if len(cluster_data) > 0:
                    plt.scatter(cluster_data['PC1'], cluster_data['GFP_fluorescence'], 
                              color=cluster_colors[cluster_id], alpha=0.4, s=8, 
                              label=f'Cluster {cluster_id}' if cluster_id < 5 else "")  # Limit legend entries
        else:
            # Fallback to gray if no clusters
            plt.scatter(filtered_sample['PC1'], filtered_sample['GFP_fluorescence'], 
                       color='gray', alpha=0.3, s=8, label=f'All variants (n={len(filtered_sample)})')
        
        # Plot selected variants as large red dots (same color for all)
        plt.scatter(final_selection['PC1'], final_selection['GFP_fluorescence'], 
                   color='red', alpha=0.8, s=30, edgecolor='black', linewidth=0.5, 
                   marker='o', label=f'Selected (n={len(final_selection)})')
        plt.xlabel('PC1')
        plt.ylabel('GFP Fluorescence')
        plt.title('PC1 vs Intensity: Filtered Library by Cluster & Selected Variants')
        plt.yscale('log')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Plot 4 (2,3,4): Selected Variants PC1 vs Intensity colored by Bin
        plt.subplot(2, 3, 4)
        if 'intensity_bin' in final_selection.columns:
            scatter = plt.scatter(final_selection['PC1'], final_selection['GFP_fluorescence'], 
                                 c=final_selection['intensity_bin'], cmap='Set3', 
                                 alpha=0.8, s=30, marker='o', edgecolor='black', linewidth=0.5,
                                 vmin=1, vmax=n_bins)
            cbar = plt.colorbar(scatter, label='Intensity Bin')
            # Set colorbar ticks to show all bins 1 through n_bins
            cbar.set_ticks(range(1, n_bins + 1))
            cbar.set_ticklabels([str(i) for i in range(1, n_bins + 1)])
            plt.xlabel('PC1')
            plt.ylabel('GFP Fluorescence')
            plt.title('Selected Variants: PC1 vs Intensity by Bin')
            plt.yscale('log')
        else:
            plt.text(0.5, 0.5, 'Intensity bins not available', ha='center', va='center', transform=plt.gca().transAxes)
            plt.title('Selected Variants by Bin (not available)')
        
        # Plot 5 (2,3,5): Selected Variants by Sequence Cluster
        plt.subplot(2, 3, 5)
        if has_clusters:
            cluster_counts = final_selection['sequence_cluster'].value_counts().sort_index()
            
            # Create complete list of all possible clusters (0 to N_CLUSTERS-1)
            all_clusters = list(range(N_CLUSTERS))
            cluster_values = [cluster_counts.get(cluster_id, 0) for cluster_id in all_clusters]
            
            # Create bar plot
            bars = plt.bar(all_clusters, cluster_values, color='lightblue', alpha=0.7, edgecolor='black')
            
            # Add value labels on top of bars
            for cluster_id, count in zip(all_clusters, cluster_values):
                plt.text(cluster_id, count + max(cluster_values) * 0.01, str(count), 
                        ha='center', va='bottom', fontweight='bold')
            
            plt.title('Selected Variants by Sequence Cluster')
            plt.xlabel('Sequence Cluster')
            plt.ylabel('Number of Selected Variants')
            plt.xticks(all_clusters)  # Show all cluster numbers
            
            # Set y-axis to start from 0 and have some headroom
            plt.ylim(0, max(cluster_values) * 1.1 if max(cluster_values) > 0 else 1)
        else:
            plt.text(0.5, 0.5, 'Clusters not available', ha='center', va='center', transform=plt.gca().transAxes)
            plt.title('Selected Variants by Cluster (not available)')
        
        # Plot 6 (2,3,6): Selected Variants Location Summary (table)
        plt.subplot(2, 3, 6)
        plt.axis('off')
        
        # Create summary statistics table
        summary_data = []
        summary_data.append(['Total Selected', len(final_selection)])
        summary_data.append(['Intensity Range', f'{final_selection["GFP_fluorescence"].min():.1f} - {final_selection["GFP_fluorescence"].max():.1f}'])
        summary_data.append(['Mean Intensity', f'{final_selection["GFP_fluorescence"].mean():.1f}'])
        summary_data.append(['Median Intensity', f'{final_selection["GFP_fluorescence"].median():.1f}'])
        
        if has_clusters:
            summary_data.append(['Clusters Represented', final_selection['sequence_cluster'].nunique()])
        
        if 'intensity_bin' in final_selection.columns:
            summary_data.append(['Bins Represented', final_selection['intensity_bin'].nunique()])
            
        if 'sequence_length' in final_selection.columns:
            summary_data.append(['Length Range (bp)', f'{final_selection["sequence_length"].min()} - {final_selection["sequence_length"].max()}'])
        
        # Create table
        table = plt.table(cellText=summary_data, 
                         colLabels=['Metric', 'Value'],
                         cellLoc='center', loc='center',
                         bbox=[0.1, 0.1, 0.8, 0.8])
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 2)
        plt.title('Selected Variants: Location Summary', fontsize=12, pad=20)

        plt.tight_layout(rect=[0, 0.03, 1, 0.93])
        print("Saving Page 3 to PDF...")
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        print("Page 3 completed.")
        
        # Page 4: Cumulative Distribution and Selection Analysis
        print("Creating Page 4: Cumulative Distribution and Selection Analysis...")
        fig = plt.figure(figsize=(16, 12))
        fig.suptitle('Cumulative Distribution and Selection Analysis', fontsize=16, y=0.95)
        
        # Plot 1: Cumulative distribution comparison log-log scale with log10 boundaries
        plt.subplot(1, 2, 1)
        
        # Get full library data (before filtering)
        if 'full_library_constructs' in locals() and len(full_library_constructs) > 0:
            full_sample = full_library_constructs
        else:
            full_sample = constructs  # Fallback to filtered data
        
        # Full library
        sorted_full = np.sort(full_sample['GFP_fluorescence'])
        cumulative_full = np.arange(1, len(sorted_full) + 1)
        plt.plot(sorted_full, cumulative_full, 'lightgray', linewidth=2, label=f'Full Library (n={len(full_sample):,})')
        
        # Filtered library
        sorted_filtered = np.sort(constructs['GFP_fluorescence'])
        cumulative_filtered = np.arange(1, len(sorted_filtered) + 1)
        plt.plot(sorted_filtered, cumulative_filtered, 'b-', linewidth=2, label=f'Filtered Library (n={len(constructs):,})')
        
        # Selected variants
        sorted_selected = np.sort(final_selection['GFP_fluorescence'])
        cumulative_selected = np.arange(1, len(sorted_selected) + 1)
        plt.plot(sorted_selected, cumulative_selected, 'r-', linewidth=3, label=f'Selected (n={len(final_selection):,})')
        
        # Add bin boundaries to show binning effect
        if USE_LOG_SCALE_BINNING:
            # Show log bin boundaries
            for boundary in log_bins[1:-1]:  # Skip first and last boundaries
                plt.axvline(boundary, color='red', linestyle=':', alpha=0.6, linewidth=1)
            title_suffix = 'with Log10 Bin Boundaries'
        else:
            # Show linear bin boundaries
            linear_bins = np.linspace(constructs['GFP_fluorescence'].min(), 
                                     constructs['GFP_fluorescence'].max(), n_bins + 1)
            for boundary in linear_bins[1:-1]:
                plt.axvline(boundary, color='blue', linestyle='--', alpha=0.6, linewidth=1)
            title_suffix = 'with Linear Bin Boundaries'
        
        plt.xlabel('GFP Fluorescence Intensity (Log Scale)')
        plt.ylabel('Cumulative Number of Variants (Log Scale)')
        plt.xscale('log')  # Add log scale to X-axis
        plt.yscale('log')
        plt.title(f'Cumulative Distribution Comparison\nFull vs Filtered vs Selected {title_suffix}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Plot 2: Cumulative distribution with bin boundaries and dual y-axes (Selected Variants)
        plt.subplot(1, 2, 2)
        
        # Prepare data - use selected variants instead of filtered library
        sorted_selected_intensities = np.sort(final_selection['GFP_fluorescence'])
        cumulative_selected_counts = np.arange(1, len(sorted_selected_intensities) + 1)
        cumulative_selected_pct = cumulative_selected_counts / len(sorted_selected_intensities) * 100
        
        # Create primary axis for cumulative count
        ax1 = plt.gca()
        line1 = ax1.plot(sorted_selected_intensities, cumulative_selected_counts, 'r-', linewidth=2, 
                        label='Cumulative Count')
        ax1.set_xlabel('GFP Fluorescence Intensity (Log Scale)')
        ax1.set_ylabel('Cumulative Number of Selected Variants', color='r')
        ax1.tick_params(axis='y', labelcolor='r')
        ax1.set_xscale('log')
        
        # Create secondary axis for cumulative percentage
        ax2 = ax1.twinx()
        line2 = ax2.plot(sorted_selected_intensities, cumulative_selected_pct, 'g-', linewidth=2, 
                        label='Cumulative Percentage')
        ax2.set_ylabel('Cumulative Percentage (%)', color='g')
        ax2.tick_params(axis='y', labelcolor='g')
        
        # Add bin boundaries as vertical lines (use same as Plot 1)
        if USE_LOG_SCALE_BINNING:
            # Use the same log bins that were created earlier
            log_min = np.log10(constructs['GFP_fluorescence'].min())
            log_max = np.log10(constructs['GFP_fluorescence'].max())
            bin_boundaries = np.logspace(log_min, log_max, n_bins + 1)
            boundary_color = 'red'
            boundary_label = 'Log Bin Boundaries'
        else:
            # Create linear bin boundaries
            bin_boundaries = np.linspace(constructs['GFP_fluorescence'].min(), 
                                       constructs['GFP_fluorescence'].max(), n_bins + 1)
            boundary_color = 'orange'
            boundary_label = 'Linear Bin Boundaries'
        
        # Plot vertical lines for bin boundaries (skip first and last which are min/max)
        boundary_lines = []
        for i, boundary in enumerate(bin_boundaries[1:-1], 1):
            line = ax1.axvline(boundary, color=boundary_color, linestyle='--', alpha=0.7, linewidth=1)
            if i == 1:  # Store first line for legend
                boundary_lines.append(line)
        
        # Create combined legend
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        all_lines = lines1 + lines2 + boundary_lines
        all_labels = labels1 + labels2 + [boundary_label]
        ax1.legend(all_lines, all_labels, loc='center right')
        
        ax1.set_title(f'Selected Variants: Cumulative Distribution\n(n={len(final_selection):,} variants)')
        ax1.grid(True, alpha=0.3)
        
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.93])
        print("Saving Page 4 to PDF...")
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        print("Page 4 completed.")
        
        # Page 5: Ready-to-Order Sequence Length Analysis
        if 'ready_to_order_final' in locals() and len(ready_to_order_final) > 0:
            print("Creating Page 5: Ready-to-Order Length Analysis...")
            fig = plt.figure(figsize=(14, 6))
            fig.suptitle('Ready-to-Order Sequence Length Analysis', fontsize=16, y=0.92)
            
            # Use the ready-to-order data for length analysis
            sequence_lengths = ready_to_order_final['Final_Sequence'].str.len()
            
            # Calculate statistics
            length_stats = {
                'count': len(sequence_lengths),
                'mean': sequence_lengths.mean(),
                'median': sequence_lengths.median(),
                'std': sequence_lengths.std(),
                'min': sequence_lengths.min(),
                'max': sequence_lengths.max()
            }
            
            print(f"Ready-to-order sequence length statistics:")
            for stat, value in length_stats.items():
                print(f"  {stat}: {value:.1f}" if isinstance(value, float) else f"  {stat}: {value}")
            
            # Plot 1: Exact Length Value Counts with all consecutive lengths
            plt.subplot(1, 2, 1)
            min_length = int(sequence_lengths.min())
            max_length = int(sequence_lengths.max())
            all_lengths = list(range(min_length, max_length + 1))
            
            # Count occurrences for each length
            length_counts = sequence_lengths.value_counts().sort_index()
            counts_for_all = [length_counts.get(length, 0) for length in all_lengths]
            
            # Create bars
            bars = plt.bar(all_lengths, counts_for_all, alpha=0.7, edgecolor='black', linewidth=0.5)
            
            # Color code bars: min=red, max=blue, others=lightblue
            for i, (length, count) in enumerate(zip(all_lengths, counts_for_all)):
                if length == min_length:
                    bars[i].set_color('red')
                elif length == max_length:
                    bars[i].set_color('blue')
                else:
                    bars[i].set_color('lightblue')
            
            # Add value labels only to bars with count > 0
            for length, count in zip(all_lengths, counts_for_all):
                if count > 0:
                    plt.text(length, count + 0.05, str(count), ha='center', va='bottom', fontsize=8)
            
            plt.xlabel('Sequence Length (bp)')
            plt.ylabel('Number of Sequences')
            plt.title(f'Exact Length Value Counts\n(n={len(sequence_lengths)}, μ={length_stats["mean"]:.1f})')
            plt.xticks(all_lengths, rotation=45 if len(all_lengths) > 15 else 0)
            
            # Add legend for colors
            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor='red', label=f'Min Length ({min_length})'),
                             Patch(facecolor='blue', label=f'Max Length ({max_length})')]
            plt.legend(handles=legend_elements, loc='upper left')
            
            # Plot 2: Sequence Component Length Breakdown
            plt.subplot(1, 2, 2)
            
            # Component lengths (user-defined, fixed values)
            primer_5_len = len(PRIMER_SITE_5_PRIME)
            bsai_5_len = len(BSAI_SITE_5_PRIME)
            bsai_3_len = len(BSAI_SITE_3_PRIME)
            primer_3_len = len(PRIMER_SITE_3_PRIME)
            
            # Calculate original sequence lengths from available data
            # Original length = Final length - primers - BsaI sites
            primer_5_len = len(ready_to_order_final['Primer_5_Prime'].iloc[0])
            bsai_5_len = len(ready_to_order_final['BsaI_5_Prime'].iloc[0])
            bsai_3_len = len(ready_to_order_final['BsaI_3_Prime'].iloc[0])
            primer_3_len = len(ready_to_order_final['Primer_3_Prime'].iloc[0])
            
            # Calculate original sequence lengths
            total_additions = primer_5_len + bsai_5_len + bsai_3_len + primer_3_len
            original_lengths = ready_to_order_final['Final_Sequence_Length'] - total_additions
            original_mean = original_lengths.mean()
            original_min = original_lengths.min()
            original_max = original_lengths.max()
            
            # Create component breakdown
            components = ['5\' Primer', '5\' BsaI', 'Original Sequence', '3\' BsaI', '3\' Primer']
            component_lengths = [primer_5_len, bsai_5_len, original_mean, bsai_3_len, primer_3_len]
            
            # Colors for each component
            colors = ['lightcoral', 'lightgreen', 'lightblue', 'lightgreen', 'lightcoral']
            
            bars = plt.bar(components, component_lengths, color=colors, alpha=0.7, edgecolor='black')
            
            # Add error bar only for Original Sequence (showing min to max range)
            original_idx = 2  # Index of "Original Sequence"
            plt.errorbar(original_idx, original_mean, 
                        yerr=[[original_mean - original_min], [original_max - original_mean]], 
                        fmt='none', color='black', capsize=5, capthick=2)
            
            # Add value labels on bars
            for i, (comp, length) in enumerate(zip(components, component_lengths)):
                if comp == 'Original Sequence':
                    # Show range and mean for original sequence
                    label = f'{original_min}-{original_max}\nμ={original_mean:.0f}'
                else:
                    label = f'{length:.0f}'
                plt.text(i, length + max(component_lengths) * 0.02, label, 
                        ha='center', va='bottom', fontweight='bold', fontsize=9)
            
            plt.ylabel('Length (bp)')
            plt.title('Sequence Component Length Breakdown')
            plt.xticks(rotation=45)
            
            plt.tight_layout(rect=[0, 0, 1, 0.88])  # Leave space for suptitle
            
            print("Saving Page 5 to PDF...")
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            print("Page 5 completed.")
        else:
            print("Ready-to-order data not available - skipping Page 5")
        
    print(f"PDF report created: {pdf_filename}")
    print("Individual plots also saved as PNG files")
    print("Visualizations completed")
else:
    print("No variants were selected - PDF report will not be generated")

# Save component usage statistics
if 'constructs' in locals() and len(constructs) > 0:
    component_stats = pd.DataFrame({
        'Component_Type': ['Promoters', 'RBS'],
        'Used_in_Library': [constructs['promoter_name'].nunique() if 'promoter_name' in constructs.columns else 0, 
                           constructs['rbs_name'].nunique() if 'rbs_name' in constructs.columns else 0],
        'Total_Constructs': [len(constructs), len(constructs)]
    })
    
    component_stats_file = os.path.join(output_dir, "07_component_usage.csv")
    component_stats.to_csv(component_stats_file, index=False)
    print(f"Component usage statistics saved to: {component_stats_file}")

# Print final summary
print("\n" + "="*60)
print("FINAL ANALYSIS SUMMARY")
print("="*60)
if 'constructs' in locals():
    print(f"Total constructs in library: {len(constructs):,}")
if 'final_selection' in locals():
    print(f"Selected diverse variants: {len(final_selection):,}")
if 'output_data_with_sites' in locals():
    print(f"Variants with BsaI sites: {len(output_data_with_sites):,}")
if 'ready_to_order_final' in locals():
    print(f"Ready-to-order variants: {len(ready_to_order_final):,}")

print(f"\nAll results saved to: {output_dir}")

# Save a copy of this script to the output directory for reproducibility
try:
    script_path = __file__  # Gets the path of the current script
    script_filename = os.path.basename(script_path)
    destination_path = os.path.join(output_dir, f"script_used_{script_filename}")
    shutil.copy2(script_path, destination_path)
    print(f"Script saved for reproducibility: {script_filename}")
except Exception as e:
    print(f"Warning: Could not save script copy: {e}")

print("Analysis completed successfully!")

# Analysis Documentation and Future Improvements
# 
# ### Analysis Summary:
# This script performs comprehensive analysis of combinatorial genetic libraries with a focus on:
# 1. **Sequence Diversity**: K-mer analysis and clustering to understand sequence space
# 2. **Expression Analysis**: Intensity-based binning for performance characterization  
# 3. **Intelligent Selection**: Diversity-based selection within intensity bins
# 4. **Practical Outputs**: BsaI-site addition and primer design for cloning
# 
# ### Key Features:
# - Automatic sequence reconstruction from component libraries
# - Comprehensive quality filtering (length, invalid sequences, restriction sites)
# - PCA-based sequence diversity analysis
# - Intensity-based binning with logarithmic or linear scales
# - K-mer distance calculation for sequence diversity
# - Intelligent variant selection balancing performance and diversity
# - Restriction site analysis for cloning compatibility
# - Ready-to-order sequence preparation with primers
# - Comprehensive visualization and reporting
# 
# ### Selection Algorithm:
# The selection process uses a multi-step approach:
# 1. **Intensity Binning**: Sequences are grouped by expression level
# 2. **Diversity Filtering**: Within each bin, sequences are selected for maximum k-mer diversity
# 3. **Restriction Site Filtering**: Sequences with problematic restriction sites are removed
# 4. **Junction Analysis**: Additional restriction sites created by concatenation are checked
# 
# ### Future Improvements:
# 1. **Machine Learning Integration**: Predict sequence properties from k-mer features
# 2. **Advanced Clustering**: Use more sophisticated clustering algorithms (HDBSCAN, spectral clustering)
# 3. **Iterative Refinement**: Use experimental results to refine selection criteria
# 4. **Alternative Diversity Metrics**: Consider other sequence features (GC content, secondary structure)
# 
# ### Usage Notes:
# - Requires 'prot' column in sd03.xlsx for fluorescence values
# - Requires 'Sequence' columns in sd01.xlsx (promoters) and sd02.xlsx (RBS)
# - Automatically handles missing or malformed data
# - K-mer diversity ensures selected sequences are maximally different within each intensity range
# 
# All analysis results, visualizations, and a copy of this script have been saved to a timestamped directory for complete reproducibility. 
