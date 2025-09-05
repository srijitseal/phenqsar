import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit.Chem.Draw import MolsToGridImage
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt
from PIL import Image
from io import BytesIO
from rdkit.Chem import Descriptors
import os
import requests
import gzip
import gdown
import tempfile

# Page configuration
st.set_page_config(
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Custom CSS for better styling
st.markdown("""
    <style>
    .main-header {
        text-align: center;
        color: #2E86C1;
        font-size: 2.5rem;
        font-weight: 700;
        margin-bottom: 0.5rem;
    }
    .subtitle {
        text-align: center;
        color: #5D6D7E;
        font-size: 1.2rem;
        margin-bottom: 2rem;
    }
    .description-box {
        background-color: #F8F9FA;
        padding: 1.5rem;
        border-radius: 10px;
        border-left: 5px solid #2E86C1;
        margin-bottom: 2rem;
    }
    .metric-container {
        background-color: #FFFFFF;
        padding: 1rem;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        text-align: center;
    }
    .download-section {
        background-color: #E8F8F5;
        padding: 1.5rem;
        border-radius: 10px;
        margin-top: 2rem;
        border: 1px solid #A3E4D7;
    }
    .stButton > button {
        width: 100%;
        background-color: #28B463;
        color: white;
        border: none;
        border-radius: 8px;
        padding: 0.75rem;
        font-weight: 600;
        transition: background-color 0.3s;
    }
    .stButton > button:hover {
        background-color: #239B56;
    }
    </style>
""", unsafe_allow_html=True)


import streamlit as st, base64, pathlib

logo_path = pathlib.Path("logo.png")
logo_b64 = base64.b64encode(logo_path.read_bytes()).decode()

st.markdown(
    f"""
    <div style="display:flex; flex-direction:column; align-items:center; gap:8px;">
        <img src="data:image/png;base64,{logo_b64}" width="300" />
        <div style="font-size:20px; font-weight:600; text-align:center;">
            The Dating App for Genes – <span style="font-weight:500;">Cell Painting Focus</span>
        </div>
    </div>
    """,
    unsafe_allow_html=True,
)

# Demo toggle at the top
st.markdown("""
<style>
.demo-toggle {
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 12px;
    padding: 12px;
    margin-bottom: 20px;
    font-weight: 600;
}

.toggle-switch {
    position: relative;
    width: 50px;
    height: 26px;
    background-color: #28a745;
    border-radius: 13px;
    cursor: not-allowed;
    transition: background-color 0.3s;
}

.toggle-slider {
    position: absolute;
    top: 2px;
    right: 2px;
    width: 22px;
    height: 22px;
    background-color: white;
    border-radius: 50%;
    transition: transform 0.3s;
    box-shadow: 0 2px 4px rgba(0,0,0,0.2);
}

.toggle-label {
    font-size: 16px;
    color: #2E86C1;
}
</style>

<div class="demo-toggle">
    <span class="toggle-label"> Select version </span>
    <div class="toggle-switch">
        <div class="toggle-slider"></div>
    </div>
    <span class="toggle-label">Using demo app with a mini dataset with 150 features</span>
</div>
""", unsafe_allow_html=True)

# Description box
st.markdown("""
<div class="description-box">
    <h3>About This Tool</h3>
    <p>This application analyzes compounds and Cell Painting data to identify potential drug hits based on phenotypic similarity between compound and gene perturbation. 
    By comparing CRISPR knockout and ORF overexpression profiles, we can find compounds that pheno-mimic or pheno-oppose the genetic perturbations in the Cell Painting assay.</p>
    <p><strong>How it works:</strong> Enter a gene name to discover compounds with similar or opposite Cell Painting profiles, 
    helping identify potential therapeutic candidates or tool compounds.</p>
</div>
""", unsafe_allow_html=True)

def display_molecule_grid(molecule_objects, jcp_ids=None, title="Molecule Grid", mols_per_row=3, sub_img_size=(400, 400)):
    """
    Displays an RDKit molecule grid in Streamlit with JCP IDs as legends.
    """
    if molecule_objects.empty:
        st.warning("No molecules to display.")
        return

    # Convert to list if it's a pandas Series
    if hasattr(molecule_objects, 'tolist'):
        molecule_objects = molecule_objects.tolist()
    
    # Prepare legends with JCP IDs if provided
    legends = None
    if jcp_ids is not None:
        if hasattr(jcp_ids, 'tolist'):
            jcp_ids = jcp_ids.tolist()
        legends = [str(jcp_id) for jcp_id in jcp_ids]

    grid_image = MolsToGridImage(molecule_objects, molsPerRow=mols_per_row, subImgSize=sub_img_size, legends=legends)
    
    buffer = BytesIO()
    grid_image.save(buffer, format="PNG")
    buffer.seek(0)
    
    st.subheader(title)
    st.image(buffer, width='stretch')

def get_drive_download_url(drive_url):
    """Convert Google Drive sharing URL to direct download URL"""
    if "drive.google.com" in drive_url:
        if "/file/d/" in drive_url:
            file_id = drive_url.split("/file/d/")[1].split("/")[0]
            return f"https://drive.google.com/uc?id={file_id}"
        elif "id=" in drive_url:
            file_id = drive_url.split("id=")[1].split("&")[0]
            return f"https://drive.google.com/uc?id={file_id}"
    return drive_url

def download_from_drive(file_id, is_gzipped=False):
    """Download file from Google Drive using gdown"""
    try:
        # Create a temporary file
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_path = temp_file.name
        
        # Download using gdown
        url = f"https://drive.google.com/uc?id={file_id}"
        gdown.download(url, temp_path, quiet=True)
        
        # Read the downloaded file
        if is_gzipped:
            df = pd.read_csv(temp_path, compression='gzip')
        else:
            df = pd.read_csv(temp_path)
        
        # Clean up temp file
        os.unlink(temp_path)
        
        return df
        
    except Exception as e:
        # Clean up temp file if it exists
        if 'temp_path' in locals() and os.path.exists(temp_path):
            os.unlink(temp_path)
        
        st.error(f"Failed to download file {file_id} using gdown: {e}")
        st.info("Please ensure the Google Drive file is set to 'Anyone with the link can view'")
        st.stop()

@st.cache_data
def load_data():
    """Load data from Google Drive"""
    try:
        # Try local files first as fallback
        if (os.path.exists("./data/mini_Median_JCP_JUMPCP_all_source_compounds_orf_crisprs_profiles_wellpos_cc_var_mad_outlier_featselect_sphering_harmony_150features.csv.gz") and 
            os.path.exists("./data/compound.csv.gz") and 
            os.path.exists("./data/knowndrughits.csv")):
            
            st.info("Using local mini demo data files")
            data = pd.read_csv("./data/mini_Median_JCP_JUMPCP_all_source_compounds_orf_crisprs_profiles_wellpos_cc_var_mad_outlier_featselect_sphering_harmony_150features.csv.gz")
            meta = pd.read_csv("./data/compound.csv.gz").dropna()
            drug_hits = pd.read_csv("./data/knowndrughits.csv")
            return data, meta, drug_hits
            
    except Exception as e:
        st.warning(f"Local files not available: {e}")
    
    # Download from Google Drive with progress
    st.info("Downloading data from Google Drive...")
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    # Download compounds first (smaller file)
    status_text.text("Loading 1/3: Compound metadata...")
    progress_bar.progress(0.1)
    meta = download_from_drive("19umFmtkt2Qww5Urj-7hv98j0b6zSv_rY", is_gzipped=True).dropna()
    progress_bar.progress(0.33)
    
    # Download drug hits (smallest file)
    status_text.text("Loading 2/3: Known drug hits...")
    drug_hits = download_from_drive("1wes2i7QMFBsMWpvw2rDddfOwPRfK2Omh", is_gzipped=False)
    progress_bar.progress(0.66)
    
    # Download main data last (largest file)
    status_text.text("Loading 3/3: Mini dataset...")
    data = download_from_drive("1U1HGmpAtCUAEXjWglr_fx9tEPzunrq0e", is_gzipped=True)
    progress_bar.progress(1.0)
    
    status_text.text("✅ All files loaded successfully!")
    progress_bar.empty()
    status_text.empty()
    
    return data, meta, drug_hits

# Initialize session state
if 'data_loaded' not in st.session_state:
    st.session_state.data_loaded = False
    st.session_state.all_genes = []

# Load data only once and cache in session state
if not st.session_state.data_loaded:
    with st.spinner("Loading data..."):
        data, meta, drug_hits = load_data()
        
        data = data.drop(columns="Metadata_JCP2022")
        data = data.groupby(["Name", "Type"]).median().reset_index()

        compounds = data[data["Type"] == "Compound"]
        crisprs = data[data["Type"] == "CRISPR"]
        orfs = data[data["Type"] == "ORF"]
        compound_descriptors = compounds.iloc[:, 2:].values

        # Get unique gene names from CRISPR and ORF data
        crispr_genes = set(crisprs["Name"].unique())
        orf_genes = set(orfs["Name"].unique())
        all_genes = sorted(list(crispr_genes.union(orf_genes)))
        
        # Store in session state
        st.session_state.data = data
        st.session_state.meta = meta  
        st.session_state.drug_hits = drug_hits
        st.session_state.compounds = compounds
        st.session_state.crisprs = crisprs
        st.session_state.orfs = orfs
        st.session_state.compound_descriptors = compound_descriptors
        st.session_state.crispr_genes = crispr_genes
        st.session_state.orf_genes = orf_genes
        st.session_state.all_genes = all_genes
        st.session_state.data_loaded = True

        st.success(f"Data loaded successfully! Found {len(crispr_genes)} CRISPR genes and {len(orf_genes)} ORF genes.")
else:
    # Use cached data from session state
    data = st.session_state.data
    meta = st.session_state.meta
    drug_hits = st.session_state.drug_hits
    compounds = st.session_state.compounds
    crisprs = st.session_state.crisprs
    orfs = st.session_state.orfs
    compound_descriptors = st.session_state.compound_descriptors
    crispr_genes = st.session_state.crispr_genes
    orf_genes = st.session_state.orf_genes
    all_genes = st.session_state.all_genes

# Input section with better styling
st.markdown("### Select Target Gene")

# Create two columns for input options
col1, col2 = st.columns([2, 1])

with col1:
    gene_name = st.selectbox(
        "Gene Symbol",
        options=[""] + all_genes,
        help=f"Select from {len(all_genes)} available genes with CRISPR/ORF data"
    )

with col2:
    st.markdown("**Available Data Types:**")
    if gene_name and gene_name in all_genes:
        available_types = []
        if gene_name in crispr_genes:
            available_types.append("CRISPR")
        if gene_name in orf_genes:
            available_types.append("ORF")
        st.info(f"✓ {', '.join(available_types)}")
    else:
        st.info("Select a gene to see available data")

if not gene_name:
    st.info("Please select a gene name to begin the analysis")
    st.stop()

st.markdown("---")

def calculate_similarity(gene_name, perturbation_type="CRISPR"):
    perturbation_data = crisprs if perturbation_type == "CRISPR" else orfs
    if perturbation_data[perturbation_data["Name"] == gene_name].empty:
        st.error(f"No {perturbation_type} data found for gene: {gene_name}")
        return None
    target = perturbation_data[perturbation_data["Name"] == gene_name].iloc[:, 2:].values
    similarities = cosine_similarity(target, compound_descriptors).flatten()
    compounds_copy = compounds.copy()
    compounds_copy["Similarity"] = similarities
    return compounds_copy.sort_values(by="Similarity", ascending=False)

st.markdown("## Cell Painting Compound Hits")
st.markdown(f"*Compounds with similar or opposite Cell Painting profiles to {gene_name} in the Cell Painting assay*")

cp_data = pd.DataFrame()

crispr_results = calculate_similarity(gene_name, 'CRISPR')
if crispr_results is not None:
    crispr_results = crispr_results[abs(crispr_results["Similarity"])>0.2]
    if not crispr_results.empty:
        top_crispr = crispr_results.head(100).copy()
        top_crispr["Calculation"] = "Most similar to CRISPR"
        bottom_crispr = crispr_results.tail(100).copy()
        bottom_crispr["Calculation"] = "Least similar to CRISPR"
        cp_data = pd.concat([cp_data, top_crispr, bottom_crispr])

orf_results = calculate_similarity(gene_name, 'ORF')
if orf_results is not None:
    orf_results = orf_results[abs(orf_results["Similarity"])>0.2]
    if not orf_results.empty:
        top_orf = orf_results.head(20).copy()
        top_orf["Calculation"] = "Most similar to ORF"
        bottom_orf = orf_results.tail(20).copy()
        bottom_orf["Calculation"] = "Least similar to ORF"
        cp_data = pd.concat([cp_data, top_orf, bottom_orf])

if not cp_data.empty:
    cp_data = cp_data.copy()
    cp_data["gene"] = gene_name

if cp_data.empty:
    st.error(f"No Cell Painting data found for gene: {gene_name}")
else:
    merged_cp_data = pd.merge(cp_data, meta[["Metadata_JCP2022", "Metadata_InChIKey", "Metadata_InChI"]], 
                           left_on="Name", right_on="Metadata_InChI").sort_values(by="Similarity", ascending=False).reset_index(drop=True)
    merged_cp_data["molecule_obj"] = merged_cp_data["Metadata_InChI"].apply(lambda x: Chem.MolFromInchi(x) if x else None)
    merged_cp_data["Metadata_Smiles"] = merged_cp_data["molecule_obj"].apply(lambda x: Chem.MolToSmiles(x) if x else None)
    merged_cp_data["MX"] = merged_cp_data["molecule_obj"].apply(lambda x: Descriptors.MolWt(x) if x else None)
    merged_cp_data = merged_cp_data[merged_cp_data["MX"] < 1000]
    merged_cp_data = merged_cp_data[["Metadata_JCP2022", "Metadata_Smiles", "Metadata_InChIKey", "Calculation", "gene", "Metadata_InChI", "molecule_obj", "Similarity"]]

    st.success(f"Found **{len(merged_cp_data)}** Cell Painting compounds with significant similarity to **{gene_name}**")
    
    # CRISPR Results
    st.markdown("### CRISPR Knockout-based Hits")
    st.markdown("*Compounds that mimic or oppose the effects of knocking out your target gene*")
    crispr_hits = merged_cp_data[merged_cp_data["Calculation"].str.contains("CRISPR", na=False)]
    if not crispr_hits.empty:
        top_crispr_hits = crispr_hits[crispr_hits["Calculation"] == "Most similar to CRISPR"].head(20)
        if not top_crispr_hits.empty:
            display_molecule_grid(
                top_crispr_hits["molecule_obj"],
                jcp_ids=top_crispr_hits["Metadata_JCP2022"],
                title="Top 20 CP Most Similar to CRISPR",
                mols_per_row=3,
                sub_img_size=(400, 400)
            )
            
            st.subheader("Top CRISPR Hit Details")
            display_cols = ["Metadata_JCP2022", "Similarity", "Calculation"]
            st.dataframe(top_crispr_hits[display_cols].head(10))
        
        bottom_crispr_hits = crispr_hits[crispr_hits["Calculation"] == "Least similar to CRISPR"].tail(20)
        if not bottom_crispr_hits.empty:
            st.subheader("CRISPR Dissimilar Hits (Potential Opposite Effect)")
            display_molecule_grid(
                bottom_crispr_hits["molecule_obj"],
                jcp_ids=bottom_crispr_hits["Metadata_JCP2022"],
                title="Bottom 20 CP Least Similar to CRISPR",
                mols_per_row=3,
                sub_img_size=(400, 400)
            )
    
    # ORF Results
    st.markdown("### ORF Overexpression-based Hits")
    st.markdown(f"*Compounds that mimic or oppose the effects of overexpressing of {gene_name} gene*")
    orf_hits = merged_cp_data[merged_cp_data["Calculation"].str.contains("ORF", na=False)]
    if not orf_hits.empty:
        top_orf_hits = orf_hits[orf_hits["Calculation"] == "Most similar to ORF"].head(20)
        if not top_orf_hits.empty:
            display_molecule_grid(
                top_orf_hits["molecule_obj"],
                jcp_ids=top_orf_hits["Metadata_JCP2022"],
                title="Top 20 CP Most Similar to ORF",
                mols_per_row=3,
                sub_img_size=(400, 400)
            )
            
            st.subheader("Top ORF Hit Details")
            display_cols = ["Metadata_JCP2022", "Similarity", "Calculation"]
            st.dataframe(top_orf_hits[display_cols].head(10))
        
        bottom_orf_hits = orf_hits[orf_hits["Calculation"] == "Least similar to ORF"].tail(20)
        if not bottom_orf_hits.empty:
            st.subheader("ORF Dissimilar Hits (Potential Opposite Effect)")
            display_molecule_grid(
                bottom_orf_hits["molecule_obj"],
                jcp_ids=bottom_orf_hits["Metadata_JCP2022"],
                title="Bottom 20 CP Least Similar to ORF",
                mols_per_row=3,
                sub_img_size=(400, 400)
            )
    else:
        st.info(f"No ORF-based hits found for {gene_name}.")

st.markdown("## Known Drug Hits Analysis")
st.markdown(f"*Comparison with known therapeutic compounds targeting {gene_name}*")
hits = drug_hits[drug_hits["Uniprot Code"]==gene_name].SMILES.to_list()

if len(hits) == 0:
    st.write(f"No known drug hits for {gene_name}")
elif cp_data.empty:
    st.write("No Cell Painting data available for comparison with known drug hits")
else:
    query_fps = []
    for query_smiles in hits:
        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol:
            morgan_gen = GetMorganGenerator(radius=2, fpSize=2048)
            query_fp = morgan_gen.GetFingerprint(query_mol)
            query_fps.append(query_fp)

    def calculate_max_similarity(row_inchi):
        mol = Chem.MolFromInchi(row_inchi)
        if mol:
            morgan_gen = GetMorganGenerator(radius=2, fpSize=2048)
            fp = morgan_gen.GetFingerprint(mol)
            similarities = [DataStructs.TanimotoSimilarity(query_fp, fp) for query_fp in query_fps]
            return max(similarities)
        else:
            return None

    merged_cp_data['Tanimoto_Similarity'] = merged_cp_data['Metadata_InChI'].apply(calculate_max_similarity)
    
    filtered_cp_hits = merged_cp_data[merged_cp_data['Tanimoto_Similarity'] > 0.55]

    if filtered_cp_hits.empty:
        st.warning("No highly similar Cell Painting compounds found compared to known drug hits.")
    else:
        st.subheader("Cell Painting Hits with Chemical Similarity > 0.55 to Known Drugs")
        st.write(f"Found {len(filtered_cp_hits)} compounds with high similarity to known drugs")
        st.dataframe(filtered_cp_hits[["Metadata_JCP2022", "Similarity", "Tanimoto_Similarity", "Calculation"]])
        
        if not filtered_cp_hits.empty:
            display_molecule_grid(
                filtered_cp_hits["molecule_obj"].head(20),
                jcp_ids=filtered_cp_hits["Metadata_JCP2022"].head(20),
                title="CP Hits Similar to Known Drugs",
                mols_per_row=3,
                sub_img_size=(400, 400)
            )

    fig, ax = plt.subplots()
    merged_cp_data["Tanimoto_Similarity"].hist(ax=ax, bins=20, color='blue', edgecolor='black')
    ax.set_title("Tanimoto Similarity Distribution - Phenotypic Hits from Cell Painting vs Known Drugs")
    ax.set_xlabel("Tanimoto Similarity")
    ax.set_ylabel("Number of CP Compounds")
    st.pyplot(fig)

st.markdown("## Summary Statistics")
if not cp_data.empty:
    # Statistics in styled containers
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("""
        <div class="metric-container">
            <h3 style="color: #2E86C1; margin: 0;">Total Hits from Cell Painting</h3>
            <h2 style="color: #34495E; margin: 10px 0;">{}</h2>
        </div>
        """.format(len(merged_cp_data)), unsafe_allow_html=True)
    
    with col2:
        high_sim = len(merged_cp_data[merged_cp_data["Similarity"] > 0.5])
        st.markdown("""
        <div class="metric-container">
            <h3 style="color: #2E86C1; margin: 0;">High Phenotypic Similarity (>0.5)</h3>
            <h2 style="color: #34495E; margin: 10px 0;">{}</h2>
        </div>
        """.format(high_sim), unsafe_allow_html=True)
    
    with col3:
        if 'Tanimoto_Similarity' in merged_cp_data.columns:
            drug_like = len(merged_cp_data[merged_cp_data['Tanimoto_Similarity'] > 0.3])
            metric_label = "Drug-like (Tanimoto >0.3)"
            metric_value = drug_like
        else:
            metric_label = "Known Drug Hits"
            metric_value = len(hits)
        
        st.markdown("""
        <div class="metric-container">
            <h3 style="color: #2E86C1; margin: 0;">{}</h3>
            <h2 style="color: #34495E; margin: 10px 0;">{}</h2>
        </div>
        """.format(metric_label, metric_value), unsafe_allow_html=True)

    # Download section with styled container
    st.markdown("""
    <div class="download-section">
        <h3 style="color: #196F3D; margin-top: 0;">Download Results</h3>
        <p>Download your analysis results as a CSV file for further investigation or sharing.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Prepare CSV data for download
    csv_data = merged_cp_data.to_csv(index=False)
    
    # Download button
    st.download_button(
        label="Download Results as CSV",
        data=csv_data,
        file_name=f"PhenQSAR_{gene_name}_hits.csv",
        mime="text/csv",
        help="Click to download all compound hits and their similarity scores"
    )
    
else:
    st.info("No data available to generate summary statistics or download results.")
