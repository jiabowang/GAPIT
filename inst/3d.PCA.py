import os
import sys
import glob
import argparse
import pandas as pd
import plotly.express as px
import numpy as np

def parse_arguments():
    \"\"\"Parse command line arguments\"\"\"
    parser = argparse.ArgumentParser(description='Generate interactive 3D visualization for GAPIT PCA results')
    parser.add_argument('-i', '--input', required=True, 
                       help='Input directory: contains GAPIT generated PCA files')
    parser.add_argument('-o', '--output', required=True, 
                       help='Output directory: save HTML files')
    parser.add_argument('-c', '--color', 
                       help='Color file: CSV format, first column is sample ID, second column is color')
    return parser.parse_args()

def ensure_dirs(input_dir, output_dir):
    \"\"\"Ensure directories exist\"\"\"
    if not os.path.isdir(input_dir):
        print(f\"Error: Input directory does not exist: {input_dir}\")
        sys.exit(1)
    os.makedirs(output_dir, exist_ok=True)

def load_color_file(color_file_path, pca_df):
    \"\"\"Load color file and map to PCA data\"\"\"
    try:
        # Try to read color file
        color_df = pd.read_csv(color_file_path)
        
        # Check column count
        if color_df.shape[1] >= 2:
            color_df = color_df.iloc[:, :2]
            color_df.columns = [\"taxa\", \"color\"]
        else:
            raise ValueError(\"Color file must have at least 2 columns\")

        # Create dictionary mapping
        color_dict = dict(zip(color_df[\"taxa\"], color_df[\"color\"]))
        
        # Add color column to PCA dataframe
        pca_df[\"color\"] = pca_df[\"taxa\"].map(color_dict)
        
        # Check missing values
        missing_colors = pca_df[\"color\"].isna().sum()
        if missing_colors > 0:
            print(f\"  Warning: {missing_colors} samples have no corresponding color, using default color\")
            pca_df.loc[pca_df[\"color\"].isna(), \"color\"] = \"black\"
        
        return pca_df
        
    except Exception as e:
        print(f\"  Warning: Failed to read color file {color_file_path}: {e}\")
        pca_df[\"color\"] = \"blue\"  # Default color
        return pca_df

def load_gapit_pca(pca_path, eigenvalues_path=None):
    \"\"\"Load GAPIT PCA output file\"\"\"
    try:
        df = pd.read_csv(pca_path)
    except:
        try:
            df = pd.read_csv(pca_path, sep='\\t')
        except:
            df = pd.read_csv(pca_path, sep='\\s+')
    
    if df.shape[1] < 4:
        raise ValueError(f\"Insufficient columns ({df.shape[1]}), at least 4 needed: sample, PC1, PC2, PC3\")
    
    if df.shape[1] == 4:
        df.columns = [\"taxa\", \"PC1\", \"PC2\", \"PC3\"]
    else:
        df = df.iloc[:, :4]
        df.columns = [\"taxa\", \"PC1\", \"PC2\", \"PC3\"]
    
    # If eigenvalue file exists, calculate variance explained
    if eigenvalues_path and os.path.exists(eigenvalues_path):
        try:
            eigenvalues = load_eigenvalues(eigenvalues_path)
            total_variance = np.sum(eigenvalues)
            pc1_var = eigenvalues[0] / total_variance * 100 if len(eigenvalues) > 0 else 0
            pc2_var = eigenvalues[1] / total_variance * 100 if len(eigenvalues) > 1 else 0
            pc3_var = eigenvalues[2] / total_variance * 100 if len(eigenvalues) > 2 else 0
            
            df.attrs['variance_explained'] = {
                'PC1': pc1_var,
                'PC2': pc2_var,
                'PC3': pc3_var
            }
        except Exception as e:
            print(f\"  Warning: Failed to read eigenvalue file {eigenvalues_path}: {e}\")
    
    return df

def load_eigenvalues(eigenvalues_path):
    \"\"\"Load eigenvalue file\"\"\"
    try:
        eigen_df = pd.read_csv(eigenvalues_path, header=0)
        eigenvalues = eigen_df.iloc[:, 0].tolist()
        return eigenvalues
    except Exception as e:
        print(f\"  Warning: Failed to parse eigenvalue file {eigenvalues_path}: {e}\")
        return []

def make_figure(pca_df, title):
    \"\"\"Create 3D scatter plot\"\"\"
    variance_explained = pca_df.attrs.get('variance_explained', {})
    
    x_label = f\"PC1 ({variance_explained.get('PC1', 0):.2f}%)\" if variance_explained else \"PC1\"
    y_label = f\"PC2 ({variance_explained.get('PC2', 0):.2f}%)\" if variance_explained else \"PC2\"
    z_label = f\"PC3 ({variance_explained.get('PC3', 0):.2f}%)\" if variance_explained else \"PC3\"
    
    if \"color\" in pca_df.columns:
        fig = px.scatter_3d(
            data_frame=pca_df,
            x=\"PC1\",
            y=\"PC2\",
            z=\"PC3\",
            hover_name=\"taxa\",
            title=title,
            labels={
                \"PC1\": x_label,
                \"PC2\": y_label,
                \"PC3\": z_label,
                #\"color\": \"Group\"
            },
            template=\"plotly_white\",
            color=\"taxa\",
            color_discrete_map=dict(zip(pca_df[\"taxa\"], pca_df[\"color\"])) 
        )
        fig.update_traces(marker=dict(size=3, opacity=0.7))
    else:
        fig = px.scatter_3d(
            data_frame=pca_df,
            x=\"PC1\",
            y=\"PC2\",
            z=\"PC3\",
            hover_name=\"taxa\",
            title=title,
            labels={
                \"PC1\": x_label,
                \"PC2\": y_label,
                \"PC3\": z_label
            },
            template=\"plotly_white\",
        )
        fig.update_traces(marker=dict(size=3, opacity=0.7, color='blue'))
    
    if variance_explained:
        total_var = sum([variance_explained.get(f'PC{i}', 0) for i in range(1, 4)])
        fig.update_layout(
            title_text=f\"{title}<br><sup>Cumulative variance explained by first 3 PCs: {total_var:.2f}%</sup>\"
        )
    
    fig.update_layout(
        scene_camera=dict(eye=dict(x=1.4, y=1.4, z=1.2)),
        margin=dict(l=0, r=0, t=80, b=0)
    )
    return fig

def export_html(fig, out_path, offline=True):
    \"\"\"Export HTML file\"\"\"
    include = \"inline\" if offline else \"cdn\"
    fig.write_html(out_path, include_plotlyjs=include, full_html=True)

def main():
    \"\"\"Main function\"\"\"
    args = parse_arguments()
    input_dir = args.input
    output_dir = args.output
    color_file = args.color
    
    ensure_dirs(input_dir, output_dir)
    
    # Find GAPIT PCA files
    # Explicitly define PCA and eigenvalue file names
    pca_path = os.path.join(input_dir, \"GAPIT.Genotype.PCA.csv\")
    eigenvalues_path = os.path.join(input_dir, \"GAPIT.Genotype.PCA_eigenvalues.csv\")
    
    if not os.path.exists(pca_path):
        print(f\"Error: PCA file not found: {pca_path}\")
        return
    
    if not os.path.exists(eigenvalues_path):
        print(\"  No eigenvalue file found, using default axis labels\")
        eigenvalues_path = None
    else:
        print(\"  Eigenvalue file found\")
    
    ok, skip, fail = 0, 0, 0
    base = os.path.basename(pca_path)
    stem, ext = os.path.splitext(base)
    out_html = os.path.join(output_dir, f\"{stem}_3DPCA.html\")
    print(f\"Processing: {base}\")

    try:
        df = load_gapit_pca(pca_path, eigenvalues_path)
        print(f\"  Successfully read data: {df.shape[0]} samples, {df.shape[1]} columns\")
        
        if color_file and os.path.exists(color_file):
            df = load_color_file(color_file, df)
            print(\"  Color information loaded\")
        else:
            print(\"  No color file provided or file not found, using default color\")
            df[\"color\"] = \"red\"
            
    except ValueError as ve:
        print(f\"  Skipped: {ve}\")
        skip += 1
    except Exception as e:
        print(f\"  Failed to read: {e}\")
        fail += 1
    else:
        try:
            title = f\"<b>Interactive 3D PCA</b> — {stem}\"
            fig = make_figure(df, title)
            export_html(fig, out_html, offline=True)
            print(f\"  Saved: {out_html}\")
            ok += 1
        except Exception as e:
            print(f\"  Failed to export: {e}\")
            fail += 1

    print(\"Processing complete!\")

if __name__ == \"__main__\":
    main()