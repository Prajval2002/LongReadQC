import streamlit as st
from io import StringIO, TextIOWrapper, BytesIO
import gzip
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import zipfile
import plotly.express as px

# Set page layout
st.set_page_config(page_title="LongRead OneClickQC", layout="wide")
st.title("📊 LongRead OneClickQC")
st.markdown("Upload Oxford Nanopore or PacBio long-read FASTQ files for QC, filtering, and summary.")

# Upload FASTQ or FASTQ.GZ
uploaded_file = st.file_uploader("Upload a .fastq or .fastq.gz file", type=["fastq", "fastq.gz"])

# Parse FASTQ records
def parse_fastq(file):
    if file.name.endswith(".gz"):
        handle = gzip.open(file, "rt")
    else:
        handle = TextIOWrapper(file)
    for record in SeqIO.parse(handle, "fastq"):
        yield record

# FASTQ to DataFrame
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

# N50 Calculation
def get_n50(lengths):
    sorted_lens = sorted(lengths, reverse=True)
    total = sum(sorted_lens)
    cum_sum = 0
    for l in sorted_lens:
        cum_sum += l
        if cum_sum >= total / 2:
            return l
    return 0

# Matplotlib Plot
def plot_hist(df, column, title, xlabel, log_scale=False):
    fig, ax = plt.subplots()
    sns.histplot(df[column], bins=50, ax=ax, kde=False)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    if log_scale:
        ax.set_yscale("log")
    fig.tight_layout()
    return fig

if uploaded_file:
    st.success("✅ File uploaded successfully. Processing...")

    # Read and store raw records once
    raw_records = list(parse_fastq(uploaded_file))
    df = fastq_to_df(raw_records)

    # Summary before filtering
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

    # Filtering UI
    st.subheader("🧹 Filtering Options")
    min_len = st.slider("Minimum Read Length", 50, 10000, 100)
    min_q = st.slider("Minimum Avg Quality", 0, 40, 10)
    filtered_df = df[(df["length"] >= min_len) & (df["avg_quality"] >= min_q)]
    filtered_ids = set(filtered_df["id"])
    st.success(f"✅ Filtered Reads: {len(filtered_df)}")

    # Interactive Plotly Histogram
    st.subheader("📊 Interactive Read Length Distribution")
    fig = px.histogram(df, x="length", nbins=100, title="Read Lengths Distribution")
    fig.update_layout(dragmode="select", xaxis_title="Length", yaxis_title="Count")
    st.plotly_chart(fig, use_container_width=True)

    # Static Histograms
    st.subheader("📉 Static Histograms")
    log_scale = st.checkbox("Log Scale (Y-Axis)", value=False)
    c1, c2, c3 = st.columns(3)
    with c1:
        st.pyplot(plot_hist(df, "length", "Read Length", "Length", log_scale))
    with c2:
        st.pyplot(plot_hist(df, "avg_quality", "Average Quality", "Quality", log_scale))
    with c3:
        st.pyplot(plot_hist(df, "gc_content", "GC Content", "GC %", log_scale))

    # FASTQ Writer
    def write_fastq_bytes(filtered_ids):
        output = StringIO()
        for rec in raw_records:
            if rec.id in filtered_ids:
                SeqIO.write(rec, output, "fastq")
        output.seek(0)
        return BytesIO(output.getvalue().encode("utf-8"))

    # GZIP writer
    def write_gzip_fastq(filtered_ids):
        buf = BytesIO()
        with gzip.GzipFile(fileobj=buf, mode='wb') as gz:
            for rec in raw_records:
                if rec.id in filtered_ids:
                    SeqIO.write(rec, TextIOWrapper(gz, write_through=True), "fastq")
        buf.seek(0)
        return buf

    # Create Zip
    def create_zip(fastq_bytes, csv_data):
        zip_buf = BytesIO()
        with zipfile.ZipFile(zip_buf, "w") as zf:
            zf.writestr("filtered.fastq", fastq_bytes.getvalue().decode("utf-8"))
            zf.writestr("summary.csv", csv_data)
        zip_buf.seek(0)
        return zip_buf

    st.subheader("⬇️ Download Outputs")
    filtered_fastq = write_fastq_bytes(filtered_ids)
    filtered_csv = filtered_df.to_csv(index=False)
    st.download_button("Download Filtered FASTQ", filtered_fastq, "filtered.fastq", "text/plain")
    st.download_button("Download Summary CSV", filtered_csv, "summary.csv", "text/csv")
    st.download_button("Download All as ZIP", create_zip(filtered_fastq, filtered_csv), "oneclickqc_outputs.zip", "application/zip")
