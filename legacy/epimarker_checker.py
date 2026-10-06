# AF22 Epigenomic Design & Enrichment Analysis Script
# 1. Download AF22 wildtype ChIP-seq bigWig files (GSE108506)
# 2. Retrieve TSS and ±30kb region for genes
# 3. saRNA search and position mapping
# 4. Interactive Plotly visualization of histone tracks with improved zooming
# 5. Enrichment analysis of ChIP signal at saRNA sites (POS vs NON controls)
# 6. Random Forest classification of control label from epigenomic features

import os
import sys
import subprocess
import requests
import pandas as pd
import pyBigWig
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from scipy import stats

# Style and visualization configuration
plt.style.use('ggplot')
sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2

# --- Configuration ---
def setup_environment():
    # Set the project directory explicitly
    script_dir = r"C:\Users\joep-\Python scripts\20250319 From PC home"
    
    # Change the working directory
    os.chdir(script_dir)
    print(f"Working directory set to: {script_dir}")
    
    # Output folder for plots
    today = datetime.now().strftime("%Y%m%d")
    folder_name = f"{today}_saRNA_designs_with_epigenomic_marker_scanning"
    output_dir = os.path.join(script_dir, folder_name)
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output plots directory: {output_dir}")
    
    # Data folder for downloaded BigWigs
    data_dir = os.path.join(script_dir, 'af22_data')
    os.makedirs(data_dir, exist_ok=True)
    
    return script_dir, output_dir, data_dir

# GEO accession mapping
datasets = {
    'H3K4me1': ('GSM2902405','GSM2902405_No0027_AF22_H3K4me1_bs50.bw'),
    'H3K4me3': ('GSM2902406','GSM2902406_No0028_AF22_H3K4me3_bs50.bw'),
    'H3K9me3': ('GSM2902407','GSM2902407_No0030_AF22_H3K9me3_bs50.bw'),
    'H3K27me3':('GSM2902408','GSM2902408_No0031_AF22_H3K27me3_bs50.bw'),
    'H3K27ac':('GSM2902409','GSM2902409_No0032_AF22_H3K27ac_bs50.bw'),
    'H3K36me3':('GSM2902410','GSM2902410_No0034_AF22_H3K36me3_bs50.bw')
}

# --- Data Retrieval Functions ---

def download_geo(acc, fname, data_dir):
    """Download BigWig files from GEO if not already present"""
    prefix = acc[:7]+'nnn'
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{prefix}/{acc}/suppl/{fname}"
    out = os.path.join(data_dir, fname)
    
    if not os.path.exists(out):
        try:
            r = requests.get(url, stream=True)
            r.raise_for_status()
            with open(out, 'wb') as f:
                for chunk in r.iter_content(1024*1024): 
                    f.write(chunk)
            print(f"Downloaded {fname}")
        except Exception as e:
            print(f"Error downloading {fname}: {str(e)}")
    else:
        print(f"File exists: {fname}")
    
    return out

def load_bigwig(path, chrom, start, end):
    """Load data from bigWig file for the specified region"""
    chrn = f"chr{chrom}" if not chrom.startswith('chr') else chrom
    try:
        with pyBigWig.open(path) as bw:
            vals = bw.values(chrn, start, end)
        return np.nan_to_num(vals)
    except Exception as e:
        print(f"Error loading bigWig {path} for {chrn}:{start}-{end}: {str(e)}")
        return np.zeros(end - start)

def get_gene_and_tss(gene):
    """Retrieve chromosome, strand, and TSS position for a gene"""
    try:
        url = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{gene}?content-type=application/json"
        r = requests.get(url)
        r.raise_for_status()
        x = r.json()
        
        if not x:
            raise ValueError(f"Gene {gene} not found")
            
        eid = next((e['id'] for e in x if e['id'].startswith('ENSG')), x[0]['id'])
        r2 = requests.get(
            f"https://rest.ensembl.org/lookup/id/{eid}?expand=1",
            headers={"Content-Type":"application/json"}
        )
        r2.raise_for_status()
        g = r2.json()
        
        chrom, strand = g['seq_region_name'], g['strand']
        txs = g.get('Transcript', [])
        ct = next((t for t in txs if t.get('is_canonical')), txs[0])
        TSS = ct['start'] if strand == 1 else ct['end']
        
        return chrom, strand, TSS
    except Exception as e:
        print(f"Error retrieving gene {gene}: {str(e)}")
        raise

def get_region(chrom, strand, TSS, flank=30000):
    """Retrieve genomic sequence around TSS"""
    try:
        s, e = max(1, TSS-flank), TSS+flank
        es = 1 if strand == 1 else -1
        r = requests.get(
            f"https://rest.ensembl.org/sequence/region/human/{chrom}:{s}..{e}:{es}?content-type=text/plain"
        )
        r.raise_for_status()
        return r.text.strip().upper(), s
    except Exception as e:
        print(f"Error retrieving region {chrom}:{TSS-flank}-{TSS+flank}: {str(e)}")
        raise

# --- saRNA Design Functions ---

def rna_to_dna(r):
    """Convert RNA sequence to DNA sequence"""
    return r.replace('U', 'T')

def validate_dna(dna):
    """Validate DNA sequence for GC content criteria"""
    if len(dna) != 19:
        return False, 0
    gc = dna.count('G') + dna.count('C')
    pct = gc / 19 * 100
    return (40 <= pct <= 60), pct

def compute_score(gc, rel, ideal=-500):
    """Compute score for saRNA based on GC content and position relative to TSS"""
    gscore = 100 - abs(gc - 50) * 2
    ddev = abs(rel - ideal)
    dscore = 100 - (ddev / abs(ideal) * 50)
    return round(0.7 * gscore + 0.3 * dscore, 2)

def find_sarna(srna, seq, start, TSS):
    """Find saRNA target position and compute its score"""
    dna = rna_to_dna(srna.strip().upper())
    idx = seq.find(dna)
    
    if idx < 0:
        return None, None
        
    abspos = start + idx
    rel = abspos - TSS
    ok, gc = validate_dna(dna)
    
    return rel, compute_score(gc, rel) if ok else (rel, 0)

# --- Visualization Functions ---

def ensure_plotly():
    """Ensure Plotly and related packages are installed"""
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        
        # Check if kaleido is installed for static image export
        try:
            import kaleido
        except ImportError:
            print("Installing kaleido for static image export...")
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'kaleido'])
    except ModuleNotFoundError:
        print("Installing plotly for visualization...")
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'plotly'])
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'kaleido'])
    
    return go, make_subplots

def visualize_tracks(chips, chrom, TSS, flank, sdata, slabels, gene, output_dir, data_dir):
    """
    Create interactive visualization of histone modification tracks with auto-zooming
    to focus on saRNA regions (±300bp around TSS and furthest saRNA)
    """
    go, make_subplots = ensure_plotly()
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    
    # Determine zoom region based on saRNA positions
    valid_pos = [pos for pos, _ in sdata if pos is not None]
    
    # Default to ±5kb if no valid saRNAs
    if not valid_pos:
        zoom_start, zoom_end = -5000, 5000
    else:
        # Find min/max positions and add margin
        min_pos = min(min(valid_pos), 0)  # Include TSS
        max_pos = max(max(valid_pos), 0)  # Include TSS
        zoom_start = min_pos - 300
        zoom_end = max_pos + 300
        
        # Ensure reasonable range (at least ±500bp)
        if zoom_end - zoom_start < 1000:
            center = (zoom_start + zoom_end) / 2
            zoom_start = center - 500
            zoom_end = center + 500
    
    x = np.arange(TSS-flank, TSS+flank) - TSS
    marks = list(chips.keys())
    n = len(marks)
    
    # Create subplots with titles on y-axis
    fig = make_subplots(
        rows=n, 
        cols=1, 
        shared_xaxes=True, 
        vertical_spacing=0.02,
        subplot_titles=[""] * n  # Empty titles, we'll add them on y-axis
    )
    
    # Enhanced color palette
    colpal = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#A65628']
    
    # Add tracks
    for i, mark in enumerate(marks, 1):
        path = os.path.join(data_dir, chips[mark][1])
        vals = load_bigwig(path, chrom, TSS-flank, TSS+flank)
        
        # Add main track line
        fig.add_trace(
            go.Scatter(
                x=x, 
                y=vals,
                mode='lines',
                line=dict(
                    color=colpal[i-1],
                    width=2,
                    shape='spline',
                    smoothing=1.2
                ),
                name=mark,
                fill='tozeroy',
                fillcolor=f'rgba({int(colpal[i-1][1:3], 16)}, {int(colpal[i-1][3:5], 16)}, {int(colpal[i-1][5:7], 16)}, 0.15)'
            ), 
            row=i, 
            col=1
        )
        
        # Add histone mark label on y-axis
        fig.update_yaxes(
            title_text=mark,
            title_font=dict(
                size=14, 
                color=colpal[i-1], 
                family="Arial, sans-serif"
            ),
            row=i, 
            col=1
        )
    
    # Add vertical lines for saRNA positions and TSS
    shapes = []
    
    # TSS line (red)
    shapes.append(dict(
        type='line',
        x0=0, x1=0,
        y0=0, y1=1,
        xref='x', 
        yref='paper',
        line=dict(
            color='#D62728',
            width=2
        )
    ))
    
    # saRNA position lines
    for (pos, score), lbl in zip(sdata, slabels):
        if pos is not None:
            shapes.append(dict(
                type='line',
                x0=pos, x1=pos,
                y0=0, y1=1,
                xref='x', 
                yref='paper',
                line=dict(
                    color='rgba(0,0,0,0.7)',
                    width=1.2,
                    dash='dash'
                )
            ))
            
            # Add a hover marker for each saRNA
            for i in range(1, n+1):
                fig.add_trace(
                    go.Scatter(
                        x=[pos],
                        y=[0],  # Will be adjusted based on track values
                        mode='markers',
                        marker=dict(
                            symbol='diamond',
                            size=8,
                            color='rgba(0,0,0,0.7)'
                        ),
                        name=f"{lbl} ({pos}bp)",
                        showlegend=(i == 1)  # Only show in legend once
                    ), 
                    row=i, 
                    col=1
                )
    
    # Update layout with cleaner styling
    fig.update_layout(
        shapes=shapes,
        template='plotly_white',
        plot_bgcolor='white',
        paper_bgcolor='white',
        height=200*n,
        title=dict(
            text=f"<b>{gene} Epigenome Profile</b>",
            font=dict(family="Arial, sans-serif", size=22, color="#333"),
            x=0.5
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=-0.15,
            xanchor="center",
            x=0.5,
            bgcolor='rgba(255,255,255,0.9)',
            bordercolor='#ddd',
            borderwidth=1
        ),
        margin=dict(l=120, r=50, t=80, b=80)
    )
    
    # Update x-axis with better labeling
    fig.update_xaxes(
        title_text='Position Relative to TSS (bp)',
        tickformat=',',
        showgrid=True,
        gridwidth=1,
        gridcolor='rgba(220,220,220,0.4)',
        zeroline=True,
        zerolinewidth=2,
        zerolinecolor='rgba(220,0,0,0.6)',
        # Add distance markers
        tickvals=[-30000, -20000, -10000, 0, 10000, 20000, 30000],
        ticktext=['-30kb', '-20kb', '-10kb', 'TSS', '+10kb', '+20kb', '+30kb'],
        range=[zoom_start, zoom_end]  # Set initial zoom range
    )
    
    # Update y-axes
    fig.update_yaxes(
        showgrid=True,
        gridwidth=1,
        gridcolor='rgba(220,220,220,0.4)',
        zeroline=False
    )
    
    # Add annotation for TSS
    fig.add_annotation(
        x=0, 
        y=1.05,
        text="TSS",
        showarrow=False,
        xref="x", 
        yref="paper",
        font=dict(family="Arial, sans-serif", size=12, color="#D62728"),
        bgcolor="rgba(255,255,255,0.8)",
        bordercolor="#D62728",
        borderwidth=1,
        borderpad=3,
        row=1, 
        col=1
    )
    
    # Save plots in multiple formats
    html_path = os.path.join(output_dir, f"{gene}_tracks.html")
    png_path = os.path.join(output_dir, f"{gene}_tracks.png")
    pdf_path = os.path.join(output_dir, f"{gene}_tracks.pdf")
    svg_path = os.path.join(output_dir, f"{gene}_tracks.svg")
    
    config = {
        'toImageButtonOptions': {
            'format': 'png',
            'filename': f'{gene}_tracks',
            'scale': 2
        },
        'displayModeBar': True,
        'responsive': True
    }
    
    # Save all formats
    fig.write_html(html_path, include_plotlyjs='cdn', config=config)
    fig.write_image(png_path, width=1200, height=200*n, scale=2)
    fig.write_image(pdf_path, width=1200, height=200*n, scale=2)
    fig.write_image(svg_path, width=1200, height=200*n, scale=2)
    
    print(f"Saved visualization for {gene} in multiple formats")
    return fig

# --- Analysis Functions ---

def enrichment_analysis(chips, chrom, TSS, flank, seqs, labels, regseq, regstart, data_dir):
    """Analyze enrichment of histone marks at saRNA positions"""
    # Find positions for all sequences
    rels = [find_sarna(s, regseq, regstart, TSS)[0] for s in seqs]
    df = pd.DataFrame({'seq': seqs, 'lab': labels, 'pos': rels})
    
    # Debug info
    print(f"\nEnrichment analysis for {chrom}:{TSS}")
    print(f"Found {len(df.dropna())} valid positions out of {len(df)} sequences")
    pos_count = df[df['lab'] == 'POS_CTRL'].pos.dropna().shape[0]
    non_count = df[df['lab'] == 'NON_CTRL'].pos.dropna().shape[0]
    print(f"POS_CTRL: {pos_count} valid positions")
    print(f"NON_CTRL: {non_count} valid positions")
    
    # Skip analysis if insufficient data
    if pos_count == 0 or non_count == 0:
        print("Skipping analysis: missing either positive or negative controls")
        return pd.DataFrame()
    
    # Calculate enrichment for each histone mark
    results = []
    for mark, (_, fname) in chips.items():
        path = os.path.join(data_dir, fname)
        
        # Extract signal values at positive and negative control positions
        posvals = [
            load_bigwig(path, chrom, TSS+int(p), TSS+int(p)+1)[0] 
            for p in df[df['lab'] == 'POS_CTRL'].pos.dropna()
        ]
        
        nonvals = [
            load_bigwig(path, chrom, TSS+int(n), TSS+int(n)+1)[0] 
            for n in df[df['lab'] == 'NON_CTRL'].pos.dropna()
        ]
        
        # Calculate means
        mpos = np.nanmean(posvals) if posvals else np.nan
        mnon = np.nanmean(nonvals) if nonvals else np.nan
        
        # Calculate statistics
        if not (np.isnan(mpos) or np.isnan(mnon)) and len(posvals) > 0 and len(nonvals) > 0:
            try:
                t_stat, p_val = stats.ttest_ind(posvals, nonvals, equal_var=False, nan_policy='omit')
            except Exception:
                t_stat, p_val = np.nan, np.nan
                
            results.append({
                'mark': mark,
                'diff': mpos - mnon,
                'pos_mean': mpos, 
                'non_mean': mnon,
                'pos_count': len(posvals),
                'non_count': len(nonvals),
                'p_value': p_val,
                't_statistic': t_stat
            })
    
    return pd.DataFrame(results)

def plot_global_differences(global_df, output_dir):
    """Create publication-quality plot of global differences between positive and negative controls"""
    if global_df.empty or global_df['diff'].isna().all():
        print("Cannot create plot: No valid difference values")
        return
    
    # Sort by difference magnitude
    sorted_df = global_df.sort_values('diff', ascending=False).reset_index(drop=True)
    
    # Set up the figure
    plt.figure(figsize=(10, 6), dpi=300)
    ax = plt.gca()
    
    # Define colors based on significance
    colors = [
        '#2C7BB6' if p < 0.05 else '#D7D7D7' 
        for p in sorted_df['p_value']
    ]
    
    # Create bars without borders
    bars = ax.bar(
        sorted_df['mark'], 
        sorted_df['diff'], 
        color=colors, 
        width=0.7,
        edgecolor='none'
    )
    
    # Add significance markers
    for i, p in enumerate(sorted_df['p_value']):
        if np.isnan(p):
            continue
            
        y_pos = sorted_df['diff'].iloc[i]
        offset = max(abs(sorted_df['diff'].max()), abs(sorted_df['diff'].min())) * 0.05
        offset = offset if y_pos >= 0 else -offset
        
        if p < 0.001:
            plt.text(i, y_pos + offset, '***', ha='center', va='center', fontsize=14, fontweight='bold')
        elif p < 0.01:
            plt.text(i, y_pos + offset, '**', ha='center', va='center', fontsize=14, fontweight='bold')
        elif p < 0.05:
            plt.text(i, y_pos + offset, '*', ha='center', va='center', fontsize=14, fontweight='bold')
    
    # Clean up the plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    
    # Add a zero line
    plt.axhline(y=0, color='black', linestyle='-', alpha=0.3, linewidth=1.5)
    
    # Improve formatting
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('Mean ChIP Signal Difference\n(POS_CTRL - NON_CTRL)', fontsize=12, fontweight='bold')
    plt.xlabel('Histone Modification', fontsize=12, fontweight='bold')
    plt.title('Differential Epigenetic Marks at saRNA Target Sites', fontsize=14, fontweight='bold')
    plt.grid(axis='y', linestyle='--', alpha=0.3)
    
    # Set reasonable y-limits
    max_abs = max(abs(sorted_df['diff'].min()), abs(sorted_df['diff'].max()))
    if max_abs > 0:
        y_limit = max_abs * 1.2
        plt.ylim(
            -y_limit if sorted_df['diff'].min() < 0 else 0, 
            y_limit if sorted_df['diff'].max() > 0 else 0
        )
    
    # Add a legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#2C7BB6', edgecolor='none', label='p < 0.05'),
        Patch(facecolor='#D7D7D7', edgecolor='none', label='p ≥ 0.05')
    ]
    ax.legend(handles=legend_elements, loc='best', frameon=True, framealpha=0.9)
    
    plt.tight_layout()
    
    # Save in multiple formats
    plt.savefig(os.path.join(output_dir, 'global_differences.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'global_differences.pdf'), format='pdf', bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'global_differences.svg'), format='svg', bbox_inches='tight')
    
    print("Saved global differences plot to output directory")

def feature_importance_plot(clf, feature_names, output_dir):
    """Create a publication-quality feature importance plot"""
    # Get feature importances and create a Series
    imp = pd.Series(clf.feature_importances_, index=feature_names).sort_values(ascending=False)
    
    plt.figure(figsize=(8, 5), dpi=300)
    ax = plt.gca()
    
    # Create bars with no border
    bars = ax.bar(
        range(len(imp)), 
        imp.values, 
        color=plt.cm.viridis(np.linspace(0.1, 0.9, len(imp))), 
        width=0.7, 
        edgecolor='none'
    )
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width()/2., 
            height + 0.01,
            f'{height:.3f}', 
            ha='center', 
            va='bottom', 
            fontsize=9, 
            fontweight='bold'
        )
    
    # Clean up plot styling
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    
    # Formatting
    plt.xticks(range(len(imp)), imp.index, rotation=45, ha='right', fontweight='bold')
    plt.grid(axis='y', linestyle='--', alpha=0.2)
    plt.ylabel('Relative Importance', fontsize=12, fontweight='bold')
    plt.xlabel('Histone Modification', fontsize=12, fontweight='bold')
    plt.title('Feature Importance in saRNA Target Classification', fontsize=14, fontweight='bold')
    
    # Add explanation
    plt.figtext(
        0.5, 0.01, 
        "Higher values indicate greater importance for predicting functional saRNA sites", 
        ha='center', 
        fontsize=10, 
        fontstyle='italic'
    )
    
    # Ensure good scaling
    max_y = max(imp.values) * 1.15
    plt.ylim(0, max_y)
    
    plt.tight_layout()
    
    # Save high-resolution versions
    plt.savefig(os.path.join(output_dir, 'feature_importances.png'), dpi=600, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'feature_importances.pdf'), format='pdf', bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'feature_importances.svg'), format='svg', bbox_inches='tight')
    
    print("Saved feature importance plot to output directory")

# --- Main Function ---

def main():
    # Setup directories
    script_dir, output_dir, data_dir = setup_environment()
    
    # Download bigWig files
    for acc, fname in datasets.values():
        download_geo(acc, fname, data_dir)
    
    # Collect data for global analysis
    all_pos_data = {}  # Store positive control data by mark 
    all_non_data = {}  # Store negative control data by mark
    all_features = []  # Store features for ML
    
    # Load design file
    try:
        design_path = os.path.join(script_dir, 'Designs.xlsx')
        if not os.path.exists(design_path):
            print(f"Design file not found at {design_path}")
            print("Please ensure 'Designs.xlsx' is in your working directory")
            return
            
        df = pd.read_excel(design_path, engine='openpyxl')
        print(f"Loaded design file with columns: {df.columns.tolist()}")
    except Exception as e:
        print(f"Error loading design file: {str(e)}")
        return
    
    # Process each gene
    for col in df.columns:
        # Extract gene name from column header
        gene = col.split()[0]
        if gene == '#' or not gene.isalnum():
            continue
            
        # Determine if column contains positive or negative controls
        is_positive = 'positive control' in str(col).lower()
        
        # Extract valid sequences
        seqs = [str(v).upper().strip() for v in df[col].dropna() if len(str(v).strip()) == 19]
        if not seqs:
            continue
            
        # Set labels based on control type
        labels = ['POS_CTRL' if is_positive else 'NON_CTRL'] * len(seqs)
        
        print(f"\nAnalyzing {gene}: {len(seqs)} sequences, type={'POS_CTRL' if is_positive else 'NON_CTRL'}")
        
        try:
            # Get TSS location
            chrom, strand, TSS = get_gene_and_tss(gene)
            print(f"  TSS location: {chrom}:{TSS} (strand {strand})")
            
            # Get genomic sequence
            regseq, regstart = get_region(chrom, strand, TSS)
            
            # Process each sequence
            for seq, label in zip(seqs, labels):
                rel, score = find_sarna(seq, regseq, regstart, TSS)
                if rel is None:
                    continue
                
                # Collect features for ML
                features = {'gene': gene, 'lab': label, 'pos': rel}
                
                # Extract signals for each histone mark
                for mark in datasets:
                    signal = load_bigwig(
                        os.path.join(data_dir, datasets[mark][1]),
                        chrom,
                        TSS + rel,
                        TSS + rel + 1
                    )[0]
                    
                    features[mark] = signal
                    
                    # Store for global analysis by control type
                    if label == 'POS_CTRL':
                        if mark not in all_pos_data:
                            all_pos_data[mark] = []
                        all_pos_data[mark].append(signal)
                    else:
                        if mark not in all_non_data:
                            all_non_data[mark] = []
                        all_non_data[mark].append(signal)
                
                # Add to features collection for ML
                all_features.append(features)
            
            # Run enrichment analysis
            sumdf = enrichment_analysis(datasets, chrom, TSS, 30000, seqs, labels, regseq, regstart, data_dir)
            
            # Create visualization
            vis_data = [find_sarna(s, regseq, regstart, TSS) for s in seqs]
            vis_labels = [f"{gene}_{i+1}" for i in range(len(seqs))]
            visualize_tracks(datasets, chrom, TSS, 30000, vis_data, vis_labels, gene, output_dir, data_dir)
            
        except Exception as e:
            print(f"Error processing {gene}: {str(e)}")
            continue
    
    # Global analysis across all genes
    print("\n=== Global Analysis Across All Genes ===")
    global_summary = []
    
    for mark in datasets.keys():
        pos_vals = all_pos_data.get(mark, [])
        non_vals = all_non_data.get(mark, [])
        
        print(f"Mark: {mark} - POS samples: {len(pos_vals)}, NON samples: {len(non_vals)}")
        
        if pos_vals and non_vals:  # Only if we have both types
            pos_mean = np.nanmean(pos_vals)
            non_mean = np.nanmean(non_vals)
            diff = pos_mean - non_mean
            
            # Calculate statistical significance
            try:
                t_stat, p_val = stats.ttest_ind(pos_vals, non_vals, equal_var=False, nan_policy='omit')
            except Exception as e:
                print(f"Error in t-test for {mark}: {str(e)}")
                t_stat, p_val = np.nan, np.nan
            
            global_summary.append({
                'mark': mark,
                'diff': diff,
                'pos_mean': pos_mean,
                'non_mean': non_mean,
                'pos_count': len(pos_vals),
                'non_count': len(non_vals),
                'p_value': p_val,
                't_statistic': t_stat
            })
    
    # Create global summary DataFrame
    global_df = pd.DataFrame(global_summary)
    
    if not global_df.empty:
        print("\nGlobal ChIP signal differences (POS_CTRL vs NON_CTRL):")
        print(global_df[['mark', 'diff', 'pos_mean', 'non_mean', 'pos_count', 'non_count', 'p_value']].sort_values('p_value'))
        
        # Create visualization of global differences
        plot_global_differences(global_df, output_dir)
    else:
        print("No global differences found - ensure you have both positive and negative controls")
    
    # Machine learning analysis
    feat_df = pd.DataFrame(all_features)
    
    if not feat_df.empty and all(mark in feat_df.columns for mark in datasets.keys()):
        print("\n=== Machine Learning Analysis ===")
        
        # Prepare features and target variable
        X = feat_df[list(datasets.keys())]
        y = feat_df['lab'].map({'POS_CTRL': 1, 'NON_CTRL': 0})
        
        # Print data summary
        print(f"Total samples: {len(feat_df)}")
        print(f"Class distribution: {y.value_counts().to_dict()}")
        
        # Check for class imbalance
        if y.value_counts().min() < 5:
            print("\nWARNING: Few samples per class detected.")
            print("For better ML results, ensure at least 5 samples of each class.")
        
        # Only proceed if we have enough data
        if len(y.unique()) > 1 and y.value_counts().min() >= 2:
            # Install scikit-learn if needed
            try:
                from sklearn.ensemble import RandomForestClassifier
                from sklearn.model_selection import train_test_split, cross_val_score
                from sklearn.metrics import classification_report
            except ModuleNotFoundError:
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'scikit-learn'])
                from sklearn.ensemble import RandomForestClassifier
                from sklearn.model_selection import train_test_split, cross_val_score
                from sklearn.metrics import classification_report
            
            # Split data for training and testing
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=0.2, stratify=y, random_state=42
            )
            
            # Train Random Forest classifier
            clf = RandomForestClassifier(n_estimators=100, random_state=42)
            clf.fit(X_train, y_train)
            
            # Evaluate performance
            preds = clf.predict(X_test)
            print("\nClassification Report:")
            print(classification_report(y_test, preds, zero_division=0))
            
            # If we have enough samples, perform cross-validation
            if len(y) >= 10:
                cv_scores = cross_val_score(clf, X, y, cv=min(5, len(y)), scoring='accuracy')
                print(f"\nCross-validation accuracy: {cv_scores.mean():.2f} ± {cv_scores.std():.2f}")
            
            # Plot feature importance
            feature_importance_plot(clf, X.columns, output_dir)
        else:
            print("Not enough samples per class for classification. Need at least 2 samples of each class.")
    else:
        print("No feature data available for machine learning analysis.")
    
    print("\nAnalysis complete. Check output directory for results.")

if __name__ == "__main__":
    main()