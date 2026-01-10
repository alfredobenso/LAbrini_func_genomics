# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 14:19:13 2026

@author: Angela Re
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.metrics import confusion_matrix, classification_report, roc_curve, auc
from sklearn.pipeline import Pipeline
import joblib
import gc  # For garbage collection
import warnings
import os


# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Set default plot style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12

# Create output directory for figures if it doesn't exist
output_dir = "output_figures"
os.makedirs(output_dir, exist_ok=True)

# Define helper functions for the analysis
def format_p_value(p_value, threshold=0.0001):
    """Helper function to consistently format p-values"""
    if p_value < threshold:
        return f"<{threshold:.4f}"
    else:
        return f"{p_value:.4f}"

def clr_transform(comp):
    """Apply centered log-ratio (CLR) transformation to a composition"""
    # Handle zeros with small pseudocount
    comp = np.array([max(val, 1e-5) for val in comp])
    # Apply CLR transform
    return np.log(comp) - np.mean(np.log(comp))

def clr_transform_df(df):
    """Apply CLR transformation to each row of a dataframe"""
    return df.apply(lambda x: np.log(x + 1e-5) - np.mean(np.log(x + 1e-5)), axis=1)


# Data loading and exploration

# Load data files
try:
    gas_df = pd.read_csv('../data/Polysome_seq/gas_data.txt', sep='\t', index_col=0)
    fructose_df = pd.read_csv('../data/Polysome_seq/fructose_data.txt', sep='\t', index_col=0)

    print(f"Loaded {gas_df.shape[0]} genes from gas dataset")
    print(f"Loaded {fructose_df.shape[0]} genes from fructose dataset")
    print(f"Columns in datasets: {', '.join(gas_df.columns)}")
except Exception as e:
    print(f"Error loading files: {e}")
    print("Please ensure the data files are in the ./data/ directory.")
    
    
# Check the expected column patterns    
    
# Expected column patterns
expected_columns = []
for rep in [1, 2, 3]:
    for phase in ['A', 'B', 'C', 'D']:
        expected_columns.append(f"{phase}_{rep}")

# Check if all expected columns are present
missing_gas = [col for col in expected_columns if col not in gas_df.columns]
missing_fructose = [col for col in expected_columns if col not in fructose_df.columns]

if missing_gas:
    print(f"Warning: Missing columns in gas file: {missing_gas}")
if missing_fructose:
    print(f"Warning: Missing columns in fructose file: {missing_fructose}")
if not missing_gas and not missing_fructose:
    print("All expected columns are present in both datasets.")    
    
    
    
# Display the first few rows of the gas dataset
print("Gas dataset preview:")
gas_df.head()

# Display the first few rows of the fructose dataset
print("Fructose dataset preview:")
fructose_df.head()

# Check sum across phases for each replicate
for rep in [1, 2, 3]:
    gas_sum = gas_df[[f"{phase}_{rep}" for phase in ['A', 'B', 'C', 'D']]].sum(axis=1)
    fructose_sum = fructose_df[[f"{phase}_{rep}" for phase in ['A', 'B', 'C', 'D']]].sum(axis=1)

    print(f"Gas replicate {rep} - Sum statistics: mean={gas_sum.mean():.4f}, min={gas_sum.min():.4f}, max={gas_sum.max():.4f}")
    print(f"Fructose replicate {rep} - Sum statistics: mean={fructose_sum.mean():.4f}, min={fructose_sum.min():.4f}, max={fructose_sum.max():.4f}")
    print()
        
    
# Replicate consistency analysis

def check_replicates(df, condition_name):
    """Analyze replicate consistency by using methods suited for compositional data."""
    print(f"\n----- {condition_name} Replicate Analysis -----")

    phases = ['A', 'B', 'C', 'D']

    # Create profile matrices for each replicate
    rep_matrices = {}
    for rep in [1, 2, 3]:
        rep_matrices[rep] = pd.DataFrame(index=df.index)
        for phase in phases:
            rep_matrices[rep][phase] = df[f"{phase}_{rep}"]

    # 1. CORRELATION ANALYSIS - using Spearman (rank-based)
    print("\n1. Spearman Correlation Analysis between Replicates:")

    # Calculate correlation for each pair of replicates
    corr_results = {}
    for phase in phases:
        corr_1_2 = stats.spearmanr(rep_matrices[1][phase], rep_matrices[2][phase])[0]
        corr_1_3 = stats.spearmanr(rep_matrices[1][phase], rep_matrices[3][phase])[0]
        corr_2_3 = stats.spearmanr(rep_matrices[2][phase], rep_matrices[3][phase])[0]

        corr_results[phase] = {
            "1_vs_2": corr_1_2,
            "1_vs_3": corr_1_3,
            "2_vs_3": corr_2_3
        }

        print(f"\nPhase {phase} (higher is better, >0.8 is good):")
        print(f"  Rep1 vs Rep2: {corr_1_2:.4f}")
        print(f"  Rep1 vs Rep3: {corr_1_3:.4f}")
        print(f"  Rep2 vs Rep3: {corr_2_3:.4f}")

    # Identify potential outlier based on correlation
    corr_outlier_votes = {1: 0, 2: 0, 3: 0}

    for phase in phases:
        results = corr_results[phase]
        # If Rep1 correlates poorly with both Rep2 and Rep3
        if results['1_vs_2'] < 0.7 and results['1_vs_3'] < 0.7 and results['2_vs_3'] > 0.7:
            corr_outlier_votes[1] += 1
        # If Rep2 correlates poorly with both Rep1 and Rep3
        elif results['1_vs_2'] < 0.7 and results['2_vs_3'] < 0.7 and results['1_vs_3'] > 0.7:
            corr_outlier_votes[2] += 1
        # If Rep3 correlates poorly with both Rep1 and Rep2
        elif results['1_vs_3'] < 0.7 and results['2_vs_3'] < 0.7 and results['1_vs_2'] > 0.7:
            corr_outlier_votes[3] += 1

    corr_outlier = max(corr_outlier_votes.items(), key=lambda x: x[1])[0] if max(corr_outlier_votes.values()) > 0 else None

    if corr_outlier:
        print(f"\nBased on correlation analysis, Replicate {corr_outlier} appears problematic")
    else:
        print("\nAll replicates show good correlation consistency")

    # 2. COMPOSITIONAL DISTANCE ANALYSIS
    print("\n2. Compositional Distance Analysis:")

    # Calculate mean composition for each replicate
    mean_comps = {}
    for rep in [1, 2, 3]:
        # Get the gene-wise mean of each phase
        mean_comps[rep] = [rep_matrices[rep][phase].mean() for phase in phases]

    # Calculate Aitchison distances between mean compositions
    aitchison_distances = {
        "1_vs_2": np.sqrt(np.sum((clr_transform(mean_comps[1]) - clr_transform(mean_comps[2])) ** 2)),
        "1_vs_3": np.sqrt(np.sum((clr_transform(mean_comps[1]) - clr_transform(mean_comps[3])) ** 2)),
        "2_vs_3": np.sqrt(np.sum((clr_transform(mean_comps[2]) - clr_transform(mean_comps[3])) ** 2))
    }

    print("Aitchison distances between replicate mean compositions (lower is better):")
    print(f"  Rep1 vs Rep2: {aitchison_distances['1_vs_2']:.4f}")
    print(f"  Rep1 vs Rep3: {aitchison_distances['1_vs_3']:.4f}")
    print(f"  Rep2 vs Rep3: {aitchison_distances['2_vs_3']:.4f}")

    # Identify potential outlier based on Aitchison distance
    rep_avg_dists = {
        1: (aitchison_distances['1_vs_2'] + aitchison_distances['1_vs_3']) / 2,
        2: (aitchison_distances['1_vs_2'] + aitchison_distances['2_vs_3']) / 2,
        3: (aitchison_distances['1_vs_3'] + aitchison_distances['2_vs_3']) / 2
    }

    aitchison_outlier = max(rep_avg_dists.items(), key=lambda x: x[1])[0]

    threshold = 0.5  # Threshold for concerning difference
    if max(rep_avg_dists.values()) > threshold:
        print(f"Replicate {aitchison_outlier} shows highest compositional distance from others")
    else:
        aitchison_outlier = None
        print("All replicates have similar compositional profiles")

    # 3. PRINCIPAL COMPONENT ANALYSIS
    print("\n3. Principal Component Analysis:")

    # Create concatenated dataset with all gene-replicate combinations
    all_data = []
    for rep in [1, 2, 3]:
        for gene in df.index:
            row_data = {
                'Gene': gene,
                'Replicate': f'Rep {rep}'
            }
            # Add each phase value
            for phase in phases:
                row_data[phase] = df.loc[gene, f"{phase}_{rep}"]

            all_data.append(row_data)

    # Create DataFrame from collected data
    pca_data = pd.DataFrame(all_data)

    # Apply CLR transformation to each gene's composition
    for i, row in pca_data.iterrows():
        phase_values = [row[phase] for phase in phases]
        try:
            clr_values = clr_transform(phase_values)
            for j, phase in enumerate(phases):
                pca_data.at[i, phase] = clr_values[j]
        except Exception as e:
            # Just leave original values
            pass

    # Extract features for PCA
    X = pca_data[phases].values

    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(X)

    # Add PCA results back to DataFrame
    pca_data['PC1'] = pca_result[:, 0]
    pca_data['PC2'] = pca_result[:, 1]

    # Plot PCA results
    plt.figure(figsize=(10, 8))
    colors = ['blue', 'green', 'red']
    markers = ['o', 's', '^']

    # Calculate centroids for each replicate
    centroids = {}
    for i, rep in enumerate(['Rep 1', 'Rep 2', 'Rep 3']):
        subset = pca_data[pca_data['Replicate'] == rep]
        plt.scatter(subset['PC1'], subset['PC2'],
                    alpha=0.5, color=colors[i], marker=markers[i], label=rep)

        # Calculate and plot centroid
        centroid = (np.mean(subset['PC1']), np.mean(subset['PC2']))
        centroids[rep] = centroid
        plt.scatter(centroid[0], centroid[1], color=colors[i],
                    marker='*', s=300, edgecolor='black')

    plt.title(f'{condition_name} - PCA of Gene Compositions by Replicate')
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
    plt.legend()
    plt.grid(alpha=0.3)
    # Save PCA replicate plot
    plt.savefig(os.path.join(output_dir, f"{condition_name.lower()}_replicate_pca.png"))
    plt.show()

    # Calculate distances between centroids
    centroid_distances = {
        '1_vs_2': np.sqrt((centroids['Rep 1'][0] - centroids['Rep 2'][0])**2 +
                         (centroids['Rep 1'][1] - centroids['Rep 2'][1])**2),
        '1_vs_3': np.sqrt((centroids['Rep 1'][0] - centroids['Rep 3'][0])**2 +
                         (centroids['Rep 1'][1] - centroids['Rep 3'][1])**2),
        '2_vs_3': np.sqrt((centroids['Rep 2'][0] - centroids['Rep 3'][0])**2 +
                         (centroids['Rep 2'][1] - centroids['Rep 3'][1])**2)
    }

    print("PCA centroid distances between replicates:")
    print(f"  Rep1 vs Rep2: {centroid_distances['1_vs_2']:.4f}")
    print(f"  Rep1 vs Rep3: {centroid_distances['1_vs_3']:.4f}")
    print(f"  Rep2 vs Rep3: {centroid_distances['2_vs_3']:.4f}")

    # Identify potential outlier based on PCA
    pca_avg_dists = {
        1: (centroid_distances['1_vs_2'] + centroid_distances['1_vs_3']) / 2,
        2: (centroid_distances['1_vs_2'] + centroid_distances['2_vs_3']) / 2,
        3: (centroid_distances['1_vs_3'] + centroid_distances['2_vs_3']) / 2
    }

    pca_outlier = max(pca_avg_dists.items(), key=lambda x: x[1])[0]

    if max(pca_avg_dists.values()) > 1.0:  # Threshold for concerning difference
        print(f"PCA shows Replicate {pca_outlier} is most distant from others")
    else:
        pca_outlier = None
        print("PCA shows all replicates cluster together well")
        
    # 4. VISUALIZATION OF REPLICATES
    # Stacked bar chart showing composition differences
    plt.figure(figsize=(12, 6))
    positions = np.arange(3)
    bottom = np.zeros(3)

    for phase in phases:
        heights = [rep_matrices[rep][phase].mean() for rep in [1, 2, 3]]
        plt.bar(positions, heights, bottom=bottom, label=f'Phase {phase}')
        bottom += heights

    plt.xticks(positions, ['Rep 1', 'Rep 2', 'Rep 3'])
    plt.ylabel('Proportion')
    plt.title(f'{condition_name} - Mean Composition by Replicate')
    plt.legend()
    # Save mean composition bar plot
    plt.savefig(os.path.join(output_dir, f"{condition_name.lower()}_replicate_composition.png"))
    plt.show()
    
    # 5. FINAL REPLICATE ASSESSMENT
    print("\n=== FINAL REPLICATE ASSESSMENT ===")

    # Combine evidence from different methods
    outlier_votes = {1: 0, 2: 0, 3: 0}

    if corr_outlier:
        outlier_votes[corr_outlier] += 1
        print(f"- Correlation analysis suggests Replicate {corr_outlier} may be problematic")

    if aitchison_outlier:
        outlier_votes[aitchison_outlier] += 1
        print(f"- Compositional distance analysis suggests Replicate {aitchison_outlier} may be problematic")

    if pca_outlier:
        outlier_votes[pca_outlier] += 1
        print(f"- PCA analysis suggests Replicate {pca_outlier} may be problematic")

    # Final decision - only exclude if it has enough votes
    if sum(outlier_votes.values()) >= 2:  # At least 2 methods agree
        problematic_rep = max(outlier_votes.items(), key=lambda x: x[1])[0]
        print(f"\nRECOMMENDATION: Exclude Replicate {problematic_rep} (most problematic)")

        # Highlight problematic replicate in composition chart
        plt.figure(figsize=(12, 6))
        positions = np.arange(3)
        bottom = np.zeros(3)

        for phase in phases:
            heights = [rep_matrices[rep][phase].mean() for rep in [1, 2, 3]]
            plt.bar(positions, heights, bottom=bottom, label=f'Phase {phase}')
            bottom += heights

        plt.xticks(positions, ['Rep 1', 'Rep 2', 'Rep 3'])
        plt.ylabel('Proportion')
        plt.title(f'{condition_name} - Mean Composition by Replicate')
        plt.legend()

        # Highlight problematic replicate
        plt.axvline(x=problematic_rep-1, color='red', linestyle='--', linewidth=2)
        plt.text(problematic_rep-1, 1.05, 'Problematic', color='red', ha='center')
        # Save highlighted bar plot
        plt.savefig(os.path.join(output_dir, f"{condition_name.lower()}_replicate_composition_problematic.png"))
        plt.show()

        return problematic_rep
    else:
        print("\nRECOMMENDATION: Keep all three replicates - no consistent evidence of problems")
        return None
    
    
    print("=== GAS CONDITION REPLICATE ANALYSIS ===")
    gas_issue = check_replicates(gas_df, "Gas")

    print("=== FRUCTOSE CONDITION REPLICATE ANALYSIS ===")
    fructose_issue = check_replicates(fructose_df, "Fructose")
    

    # Statistical difference analysis between growth conditions
    # Compare the gas and fructose conditions to identify statistical differences 
    # in their polysomal profiles.
    # We calculate means across valid replicates (excluding any problematic
    #ones identified above).

    def calculate_means(df, excluded_rep=None):
        """Calculate mean values across valid replicates, excluding any problematic one"""
        means = pd.DataFrame(index=df.index)
        for phase in ['A', 'B', 'C', 'D']:
            if excluded_rep:
                # Calculate means excluding the problematic replicate
                valid_reps = [i for i in [1, 2, 3] if i != excluded_rep]
                means[phase] = df[[f"{phase}_{r}" for r in valid_reps]].mean(axis=1)
            else:
                # Calculate means with all replicates
                means[phase] = df[[f"{phase}_1", f"{phase}_2", f"{phase}_3"]].mean(axis=1)
        return means

    # Calculate means using only valid replicates
    gas_means = calculate_means(gas_df, gas_issue)
    fructose_means = calculate_means(fructose_df, fructose_issue)

    # Ensure both datasets have the same genes
    common_genes = list(set(gas_means.index).intersection(set(fructose_means.index)))
    
    
    
    # Three complementary statistical methods are used to define differentially translated genes:
    #1. Bootstrap confidence intervals: it identifies significantly different phases 
    #for each gene by defining confidence intervals robust to outliers.
    #2. Variance-weighted Aitchison distance: it weights differences by their variance 
    #to prioritize consistent changes by downweighting phases with high variability.
    #3. Absolute translation difference: it focuses on translation differences by
    # borroings information across genes to improve variance estimates. 


# Setup and initialization for robust differential gene analysis
import time
from datetime import datetime, timedelta
from tqdm.notebook import tqdm  # For progress bars in Jupyter

# Start overall timing
overall_start_time = time.time()

print("\n\n====================================================================")
print("ROBUST DIFFERENTIAL GENE IDENTIFICATION")
print("====================================================================")


# Create containers for our results
diff_genes_methods = {}
top_n = 100  # Number of top genes to display from each method
execution_times = {}  # Track execution time for each method

# Extract valid replicates information for later use

valid_gas_reps = [i for i in [1, 2, 3] if i != gas_issue]
valid_fructose_reps = [i for i in [1, 2, 3] if i != fructose_issue]

print(f"Using {len(valid_gas_reps)} valid gas replicates: {valid_gas_reps}")
print(f"Using {len(valid_fructose_reps)} valid fructose replicates: {valid_fructose_reps}")
print(f"Analysis will examine {len(common_genes)} common genes between datasets")


# METHOD 1: Bootstrap confidence intervals for phase differences
print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Starting METHOD 1: Bootstrap confidence intervals for phase differences")
method1_start = time.time()

# Initialize result storage
bootstrap_results = {}
significant_genes = []

# Number of bootstrap iterations
n_boot = 1000
print(f"  Running {n_boot} bootstrap iterations for each gene and phase...")

# Implemented for phase-specific differences
for phase in ['A', 'B', 'C', 'D']:
    print(f"  Analyzing phase {phase}...")

    for gene in tqdm(common_genes, desc=f"Bootstrap phase {phase}"):
        # Get original values for this gene and phase
        gas_vals = []
        for rep in valid_gas_reps:
            gas_vals.append(gas_df.loc[gene, f"{phase}_{rep}"])

        fructose_vals = []
        for rep in valid_fructose_reps:
            fructose_vals.append(fructose_df.loc[gene, f"{phase}_{rep}"])

        # Calculate observed difference
        observed_diff = np.mean(gas_vals) - np.mean(fructose_vals)

        # Bootstrap sampling
        boot_diffs = []

        for _ in range(n_boot):
            # Resample with replacement
            if len(gas_vals) > 1:
                boot_gas = np.random.choice(gas_vals, size=len(gas_vals), replace=True)
            else:
                boot_gas = gas_vals

            if len(fructose_vals) > 1:
                boot_fructose = np.random.choice(fructose_vals, size=len(fructose_vals), replace=True)
            else:
                boot_fructose = fructose_vals

            # Calculate bootstrap difference
            boot_diff = np.mean(boot_gas) - np.mean(boot_fructose)
            boot_diffs.append(boot_diff)

        # Calculate confidence interval
        boot_diffs = np.array(boot_diffs)
        ci_lower = np.percentile(boot_diffs, 2.5)
        ci_upper = np.percentile(boot_diffs, 97.5)

        # Check if 0 is in the confidence interval
        significant = (ci_lower > 0 and observed_diff > 0) or (ci_upper < 0 and observed_diff < 0)

        # Store results
        if gene not in bootstrap_results:
            bootstrap_results[gene] = {'significant_phases': [], 'max_diff': 0}

        if significant:
            bootstrap_results[gene]['significant_phases'].append(phase)

        # Track maximum absolute difference for this gene
        abs_diff = abs(observed_diff)
        if abs_diff > bootstrap_results[gene]['max_diff']:
            bootstrap_results[gene]['max_diff'] = abs_diff

# Count genes with at least one significant phase difference
print("  Summarizing bootstrap results...")
for gene, result in bootstrap_results.items():
    if len(result['significant_phases']) > 0:
        significant_genes.append(gene)

print(f"Found {len(significant_genes)} genes with at least one significant phase difference")

# Sort genes by maximum phase difference and number of significant phases
bootstrap_df = pd.DataFrame([
    {'gene': gene,
     'num_sig_phases': len(result['significant_phases']),
     'max_diff': result['max_diff']}
    for gene, result in bootstrap_results.items()
])

# Sort first by number of significant phases, then by max difference
bootstrap_df = bootstrap_df.sort_values(['num_sig_phases', 'max_diff'],
                                        ascending=[False, False])

diff_genes_methods['bootstrap'] = list(bootstrap_df.head(top_n)['gene'])

# Calculate execution time
method1_time = time.time() - method1_start
execution_times['Method 1'] = method1_time
print(f"Top 5 genes: {', '.join(diff_genes_methods['bootstrap'][:5])}")
print(f"[{datetime.now().strftime('%H:%M:%S')}] Method 1 completed in {timedelta(seconds=int(method1_time))}")

# METHOD 2: Variance-weighted Aitchison distance
print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Starting METHOD 2: Variance-weighted Aitchison distance")
method2_start = time.time()

# First, calculate CLR for each replicate
print("  Calculating CLR for each replicate...")
gas_clr_byRep = {}
fructose_clr_byRep = {}

for rep in valid_gas_reps:
    # Extract data for this replicate
    rep_data = pd.DataFrame(index=gas_df.index)
    for phase in ['A', 'B', 'C', 'D']:
        rep_data[phase] = gas_df[f"{phase}_{rep}"]
    # Apply CLR transformation
    gas_clr_byRep[rep] = clr_transform_df(rep_data)

for rep in valid_fructose_reps:
    # Extract data for this replicate
    rep_data = pd.DataFrame(index=fructose_df.index)
    for phase in ['A', 'B', 'C', 'D']:
        rep_data[phase] = fructose_df[f"{phase}_{rep}"]
    # Apply CLR transformation
    fructose_clr_byRep[rep] = clr_transform_df(rep_data)

# Calculate weighted Aitchison distances
print("  Calculating variance-weighted distances...")
weighted_distances = {}

for gene in tqdm(common_genes, desc="Weighted distances"):
    if gene not in weighted_distances:
        weighted_distances[gene] = 0

    # Get CLR values for each phase in each replicate
    gas_clrs = {}
    fructose_clrs = {}

    for phase in ['A', 'B', 'C', 'D']:
        gas_clrs[phase] = [gas_clr_byRep[rep].loc[gene, phase] for rep in valid_gas_reps]
        fructose_clrs[phase] = [fructose_clr_byRep[rep].loc[gene, phase] for rep in valid_fructose_reps]

    # Calculate mean and variance for each phase
    gas_clr_means = {phase: np.mean(vals) for phase, vals in gas_clrs.items()}
    fructose_clr_means = {phase: np.mean(vals) for phase, vals in fructose_clrs.items()}

    # Calculate variances, handling single replicate case
    gas_clr_vars = {}
    fructose_clr_vars = {}

    for phase in ['A', 'B', 'C', 'D']:
        if len(gas_clrs[phase]) > 1:
            gas_clr_vars[phase] = np.var(gas_clrs[phase], ddof=1)
        else:
            gas_clr_vars[phase] = 1.0  # Default variance for single replicate

        if len(fructose_clrs[phase]) > 1:
            fructose_clr_vars[phase] = np.var(fructose_clrs[phase], ddof=1)
        else:
            fructose_clr_vars[phase] = 1.0  # Default variance for single replicate

    # Calculate weighted distance
    weighted_dist = 0

    for phase in ['A', 'B', 'C', 'D']:
        # Difference between means
        diff = gas_clr_means[phase] - fructose_clr_means[phase]

        # Pooled variance
        pooled_var = (gas_clr_vars[phase] / len(gas_clrs[phase]) +
                      fructose_clr_vars[phase] / len(fructose_clrs[phase]))

        # Add to weighted distance - similar to t-statistic formula
        if pooled_var > 0:
            weighted_dist += (diff ** 2) / pooled_var

    weighted_distances[gene] = weighted_dist

# Sort genes by weighted distance
weighted_results = pd.DataFrame({
    'gene': list(weighted_distances.keys()),
    'weighted_distance': list(weighted_distances.values())
})

weighted_results = weighted_results.sort_values('weighted_distance', ascending=False)
diff_genes_methods['weighted_aitchison'] = list(weighted_results.head(top_n)['gene'])

# Calculate execution time
method2_time = time.time() - method2_start
execution_times['Method 2'] = method2_time
print(f"Top 5 genes: {', '.join(diff_genes_methods['weighted_aitchison'][:5])}")
print(f"[{datetime.now().strftime('%H:%M:%S')}] Method 2 completed in {timedelta(seconds=int(method2_time))}")


# METHOD 3: Absolute Translation Status Difference 
print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Starting METHOD 3: Absolute Translation Status Difference")
method3_start = time.time()

print("  Calculating translation status from raw data...")
# Calculate translated fraction (B+C+D) for each gene from scratch
gas_translated_new = pd.DataFrame(index=gas_means.index)
gas_translated_new['non_translated'] = gas_means['A']
gas_translated_new['translated'] = gas_means['B'] + gas_means['C'] + gas_means['D']
gas_translated_new['translation_ratio'] = gas_translated_new['translated'] / (gas_translated_new['non_translated'] + 1e-10)

fructose_translated_new = pd.DataFrame(index=fructose_means.index)
fructose_translated_new['non_translated'] = fructose_means['A']
fructose_translated_new['translated'] = fructose_means['B'] + fructose_means['C'] + fructose_means['D']
fructose_translated_new['translation_ratio'] = fructose_translated_new['translated'] / (fructose_translated_new['non_translated'] + 1e-10)

# Find common genes
common_genes_ts = sorted(list(set(gas_translated_new.index) & set(fructose_translated_new.index)))
print(f"  Found {len(common_genes_ts)} common genes for translation status comparison")

# Calculate absolute difference in translation status
translation_diff_df = pd.DataFrame(index=common_genes_ts)
translation_diff_df['gas_status'] = gas_translated_new.loc[common_genes_ts, 'translation_ratio']
translation_diff_df['fructose_status'] = fructose_translated_new.loc[common_genes_ts, 'translation_ratio']
translation_diff_df['diff'] = translation_diff_df['gas_status'] - translation_diff_df['fructose_status']
translation_diff_df['abs_diff'] = np.abs(translation_diff_df['diff'])

# Sort by absolute difference
translation_diff_df = translation_diff_df.sort_values('abs_diff', ascending=False)

# Select top genes
diff_genes_methods['abs_translation_diff'] = list(translation_diff_df.head(top_n).index)

# Calculate execution time
method3_time = time.time() - method3_start
execution_times['Method 3'] = method3_time

print(f"  Largest absolute translation ratio difference: {translation_diff_df['abs_diff'].max():.4f}")
print(f"  Mean absolute translation ratio difference: {translation_diff_df['abs_diff'].mean():.4f}")
print(f"Top 5 genes by translation ratio difference: {', '.join(diff_genes_methods['abs_translation_diff'][:5])}")
print(f"[{datetime.now().strftime('%H:%M:%S')}] Method 3 completed in {timedelta(seconds=int(method3_time))}")



# Save results

# Method abbreviations
method_abbrevs = {
    'abs_translation_diff': 'ab',
    'bootstrap': 'BS',
    'weighted_aitchison': 'WA'
}


# Export results from each individual method

for method_name, genes in diff_genes_methods.items():
    if method_name == 'abs_translation_diff':
        translation_diff_df.to_csv(f"{method_name}_results.csv")   
    elif method_name == 'bootstrap':
        bootstrap_df.to_csv(f"{method_name}_results.csv")
    elif method_name == 'weighted_aitchison':
        weighted_results.to_csv(f"{method_name}_results.csv")
    print(f"  - Exported {method_name} results to '{method_name}_results.csv'")
    
    
    