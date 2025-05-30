import streamlit as st
from io import StringIO, TextIOWrapper, BytesIO
import gzip
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import zipfile
import plotly.express as px

st.set_page_config(page_title="LongRead OneClickQC", layout="wide")

st.title("📊 LongRead OneClickQC")
st.markdown("Upload Oxford Nanopore or PacBio long-read FASTQ files for QC, filtering, and summary.")

# Upload file
uploaded_file = st.file_uploader("Upload a .fastq or .fastq.gz file", type=["fastq", "fastq.gz"])

# Helper: FASTQ reader
def parse_fastq(file):
    if file.name.endswith(".gz"):
        handle = gzip.open(file, "rt")
    else:
        handle = TextIOWrapper(file)
    return list(SeqIO.parse(handle, "fastq"))

# Helper: Convert to DataFrame
def fastq_to_df(records):
    data = []
    for rec in records:
        quals = rec.letter_annotations.get("phred_quality", [])
        avg_q = sum(quals)/len(quals) if quals else 0
        gc = 100 * (rec.seq.count("G") + rec.seq.count("C")) / len(rec.seq)
        data.append({
            "id": rec.id,
            "length": len(rec.seq),
            "avg_quality": avg_q,
            "gc_content": gc
        })
    return pd.DataFrame(data)

# Helper: N50
def get_n50(lengths):
    sorted_lens = sorted(lengths, reverse=True)
    total_bases = sum(sorted_lens)
    cum = 0
    for l in sorted_lens:
        cum += l
        if cum >= total_bases / 2:
            return l
    return 0

# Helper: Plotting (matplotlib)
def plot_hist(df, column, title, xlabel, log_scale=False):
    fig, ax = plt.subplots()
    sns.histplot(df[column], bins=50, ax=ax)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    if log_scale:
        ax.set_yscale('log')
    return fig

if uploaded_file:
    st.success("✅ File uploaded successfully.")
    records = parse_fastq(uploaded_file)
    df = fastq_to_df(records)

    # Summary stats
    st.subheader("📈 Summary Statistics (Before Filtering)")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Reads", len(df))
        st.metric("Total Bases", f"{df['length'].sum():,}")
    with col2:
        st.metric("Avg Read Length", f"{df['length'].mean():.2f} bp")
        st.metric("Max Length", f"{df['length'].max()} bp")
    with col3:
        st.metric("N50", f"{get_n50(df['length'])} bp")
        st.metric("Median Length", f"{int(df['length'].median())} bp")

    # Filters
    st.subheader("🧹 Filtering Options")
    min_len = st.slider("Minimum Read Length", 50, 10000, 100)
    min_q = st.slider("Minimum Avg Quality", 40, 80, 35)
    filtered_df = df[(df["length"] >= min_len) & (df["avg_quality"] >= min_q)]
    filtered_ids = set(filtered_df["id"])
    st.success(f"✅ Filtered Reads: {len(filtered_df)}")

    # Interactive Plotly Histogram
    st.subheader("📊 Interactive Read Length Distribution")
    fig = px.histogram(df, x="length", nbins=100, title="Read Lengths Distribution (Zoom & Select to Filter)")
    fig.update_layout(dragmode="select")
    st.plotly_chart(fig, use_container_width=True)

    # Static Visualizations
    st.subheader("📉 Static Histograms")
    log_scale = st.checkbox("Log Scale (Y-Axis for Histograms)", value=False)
    col1, col2, col3 = st.columns(3)
    with col1:
        st.pyplot(plot_hist(df, "length", "Read Lengths", "Length", log_scale))
    with col2:
        st.pyplot(plot_hist(df, "avg_quality", "Average Quality", "Quality", log_scale))
    with col3:
        st.pyplot(plot_hist(df, "gc_content", "GC Content", "GC%", log_scale))

    # Downloads
    st.subheader("⬇️ Download Outputs")

    # FASTQ writer
    def write_fastq_bytes(filtered_ids):
        txt_buf = StringIO()
        for rec in records:
            if rec.id in filtered_ids:
                SeqIO.write(rec, txt_buf, "fastq")
        txt_buf.seek(0)
        byte_buf = BytesIO(txt_buf.getvalue().encode('utf-8'))
        byte_buf.seek(0)
        return byte_buf

    filtered_fastq_bytes = write_fastq_bytes(filtered_ids)
    csv_data = filtered_df.to_csv(index=False)

    # ZIP Download
    def create_zip():
        zip_buf = BytesIO()
        with zipfile.ZipFile(zip_buf, "w") as z:
            z.writestr("filtered.fastq", filtered_fastq_bytes.getvalue().decode('utf-8'))
            z.writestr("summary.csv", csv_data)
        zip_buf.seek(0)
        return zip_buf

    st.download_button("Download Filtered FASTQ", filtered_fastq_bytes, "filtered.fastq", "text/plain")
    st.download_button("Download Summary CSV", csv_data, "summary.csv", "text/csv")
    st.download_button("Download All as ZIP", create_zip(), "oneclickqc_outputs.zip", "application/zip")

    