import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats
import numpy as np

adata = ad.read_h5ad("results/finalized/merged.h5ad")

# Check observation columns
print("Observation columns:")
print(adata.obs.columns.tolist())

# Check variable columns  
print("\nVariable columns:")
print(adata.var.columns.tolist())

# Check obsm keys (embeddings)
print("\nObsm keys:")
print(list(adata.obsm.keys()))

# Create filename-safe version of celltypist labels
# def clean_filename(name):
#     name = str(name)
#     # Replace forward slash and other problematic characters
#     name = name.replace('/', '_')
#     name = name.replace('+', 'pos')
#     name = name.replace('-', 'neg')
#     name = re.sub(r'[\\:*?"<>|]', '_', name)
#     return name

# adata.obs['celltypist_clean'] = adata.obs['celltypist:Human_PF_Lung'].apply(clean_filename)

# # Show the mapping
# print("Label mapping:")
# mapping = dict(zip(adata.obs['celltypist:Human_PF_Lung'].unique(), 
#                   adata.obs['celltypist_clean'].unique()))
# for orig, clean in mapping.items():
#     print(f"'{orig}' -> '{clean}'")

# adata.write_h5ad("results/finalized/merged_cleaned.h5ad")

# Create disease status grouping
def get_disease_status(batch):
    if 'NOR' in batch:
        return 'Normal'
    elif 'SSCLOW' in batch or 'SSCUP' in batch:
        return 'Disease'
    else:
        return 'Other'

adata.obs['disease_status'] = adata.obs['batch'].apply(get_disease_status)

# Create combined grouping for dodged plot
adata.obs['celltype_disease'] = adata.obs['celltypist:Human_PF_Lung'].astype(str) + '_' + adata.obs['disease_status'].astype(str)

# Create dodged box plots
fig, axes = plt.subplots(1, 2, figsize=(30, 12))

# CEACAM6 - filter cell types with at least 50 cells
df_ceacam6 = pd.DataFrame({
    'expression': adata[:, 'CEACAM6'].X.toarray().flatten(),
    'celltype': adata.obs['celltypist:Human_PF_Lung'],
    'disease_status': adata.obs['disease_status']
})
df_ceacam6 = df_ceacam6[df_ceacam6['disease_status'].isin(['Normal', 'Disease'])]

# Filter cell types with sufficient cells
celltype_counts = df_ceacam6['celltype'].value_counts()
print(f"Cell type counts:\n{celltype_counts}")
valid_celltypes = celltype_counts[celltype_counts >= 50].index
df_ceacam6 = df_ceacam6[df_ceacam6['celltype'].isin(valid_celltypes)]

# Option 1: Use violinplot with a wider width
sns.violinplot(data=df_ceacam6, x='celltype', y='expression', hue='disease_status', 
               ax=axes[0], split=False, inner='quart', width=0.8, linewidth=1.5)
axes[0].set_title('CEACAM6 Expression: Normal vs Disease', fontsize=16)
axes[0].tick_params(axis='x', rotation=45, labelsize=12)
axes[0].tick_params(axis='y', labelsize=12)

# Option 2: Use boxplot instead for sparse data
sns.boxplot(data=df_ceacam6, x='celltype', y='expression', hue='disease_status', 
            ax=axes[0], width=0.8, linewidth=1.5)
axes[0].set_title('CEACAM6 Expression: Normal vs Disease', fontsize=16)
axes[0].tick_params(axis='x', rotation=45, labelsize=12)
axes[0].tick_params(axis='y', labelsize=12)

# CEACAM5 - same filtering
df_ceacam5 = pd.DataFrame({
    'expression': adata[:, 'CEACAM5'].X.toarray().flatten(),
    'celltype': adata.obs['celltypist:Human_PF_Lung'],
    'disease_status': adata.obs['disease_status']
})
df_ceacam5 = df_ceacam5[df_ceacam5['disease_status'].isin(['Normal', 'Disease'])]
df_ceacam5 = df_ceacam5[df_ceacam5['celltype'].isin(valid_celltypes)]

sns.boxplot(data=df_ceacam5, x='celltype', y='expression', hue='disease_status',
            ax=axes[1], width=0.8, linewidth=1.5)
axes[1].set_title('CEACAM5 Expression: Normal vs Disease', fontsize=16)
axes[1].tick_params(axis='x', rotation=45, labelsize=12)
axes[1].tick_params(axis='y', labelsize=12)

plt.tight_layout()
plt.savefig('CEACAM_dodged_boxplots.png', dpi=300, bbox_inches='tight')
plt.show()

# CEACAM6 separate plot with better label handling
fig, ax = plt.subplots(1, 1, figsize=(25, 10))
sns.boxplot(data=df_ceacam6, x='celltype', y='expression', hue='disease_status', 
            ax=ax, width=0.8, linewidth=1.5)
ax.set_title('CEACAM6 Expression: Normal vs Disease', fontsize=16)
ax.tick_params(axis='x', rotation=90, labelsize=10)  # Vertical rotation
ax.tick_params(axis='y', labelsize=12)
ax.set_xlabel('Cell Type', fontsize=14)
ax.set_ylabel('Expression', fontsize=14)
plt.tight_layout()
plt.subplots_adjust(bottom=0.25)  # Add more space for labels
plt.savefig('CEACAM6_boxplot.png', dpi=300, bbox_inches='tight')
plt.show()

# CEACAM5 separate plot with better label handling
fig, ax = plt.subplots(1, 1, figsize=(25, 10))
sns.boxplot(data=df_ceacam5, x='celltype', y='expression', hue='disease_status',
            ax=ax, width=0.8, linewidth=1.5)
ax.set_title('CEACAM5 Expression: Normal vs Disease', fontsize=16)
ax.tick_params(axis='x', rotation=90, labelsize=10)  # Vertical rotation
ax.tick_params(axis='y', labelsize=12)
ax.set_xlabel('Cell Type', fontsize=14)
ax.set_ylabel('Expression', fontsize=14)
plt.tight_layout()
plt.subplots_adjust(bottom=0.25)  # Add more space for labels
plt.savefig('CEACAM5_boxplot.png', dpi=300, bbox_inches='tight')
plt.show()

def add_stat_annotation(ax, data, x_col, y_col, hue_col, test='mannwhitneyu'):
    """Add statistical significance annotations to cell type labels"""
    cell_types = data[x_col].unique()
    hue_levels = data[hue_col].unique()
    
    if len(hue_levels) != 2:
        return
    
    # Create mapping of cell type to significance
    sig_mapping = {}
    
    for cell_type in cell_types:
        # Get data for each group
        group1 = data[(data[x_col] == cell_type) & (data[hue_col] == hue_levels[0])][y_col]
        group2 = data[(data[x_col] == cell_type) & (data[hue_col] == hue_levels[1])][y_col]
        
        if len(group1) < 3 or len(group2) < 3:
            sig_mapping[cell_type] = ''
            continue
            
        # Perform statistical test
        if test == 'mannwhitneyu':
            stat, p_value = stats.mannwhitneyu(group1, group2, alternative='two-sided')
        elif test == 'ttest':
            stat, p_value = stats.ttest_ind(group1, group2)
        
        # Determine significance level
        if p_value < 0.001:
            sig_mapping[cell_type] = ' ***'
        elif p_value < 0.01:
            sig_mapping[cell_type] = ' **'
        elif p_value < 0.05:
            sig_mapping[cell_type] = ' *'
        else:
            sig_mapping[cell_type] = ''
    
    # Get current labels and add significance stars
    current_labels = [tick.get_text() for tick in ax.get_xticklabels()]
    new_labels = []
    
    for label_text in current_labels:
        star = sig_mapping.get(label_text, '')
        new_labels.append(f"{label_text}{star}")
    
    # Set the new labels with stars
    ax.set_xticklabels(new_labels, rotation=90, fontsize=10)

# CEACAM6 separate plot with statistics
fig, ax = plt.subplots(1, 1, figsize=(25, 10))
sns.boxplot(data=df_ceacam6, x='celltype', y='expression', hue='disease_status', 
            ax=ax, width=0.8, linewidth=1.5)

# Add statistical annotations to x-axis labels
add_stat_annotation(ax, df_ceacam6, 'celltype', 'expression', 'disease_status')

ax.set_title('CEACAM6 Expression: Normal vs Disease', fontsize=16)
ax.tick_params(axis='y', labelsize=12)
ax.set_xlabel('Cell Type', fontsize=14)
ax.set_ylabel('Expression', fontsize=14)
plt.tight_layout()
plt.subplots_adjust(bottom=0.25)
plt.savefig('CEACAM6_boxplot_with_stats.png', dpi=300, bbox_inches='tight')
plt.show()

# CEACAM5 separate plot with statistics
fig, ax = plt.subplots(1, 1, figsize=(25, 10))
sns.boxplot(data=df_ceacam5, x='celltype', y='expression', hue='disease_status',
            ax=ax, width=0.8, linewidth=1.5)

# Add statistical annotations to x-axis labels
add_stat_annotation(ax, df_ceacam5, 'celltype', 'expression', 'disease_status')

ax.set_title('CEACAM5 Expression: Normal vs Disease', fontsize=16)
ax.tick_params(axis='x', rotation=90, labelsize=10)
ax.tick_params(axis='y', labelsize=12)
ax.set_xlabel('Cell Type', fontsize=14)
ax.set_ylabel('Expression', fontsize=14)
plt.tight_layout()
plt.subplots_adjust(bottom=0.25)
plt.savefig('CEACAM5_boxplot_with_stats.png', dpi=300, bbox_inches='tight')
plt.show()
