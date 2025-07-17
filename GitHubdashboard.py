#imports
import pandas as pd
import plotly.express as px
import streamlit as st 
from Counts import getaacounts
from extractgps4 import extract_helix_4mers 
from extractgps3 import extract_helix_trimers
from extractgps2 import extract_helix_dimers
amino_acids = ['H', 'R', 'K', 'I', 'F', 'L', 'W', 'A', 'M',
               'P', 'C', 'N', 'V', 'G', 'S', 'Q', 'Y', 'D', 'E', 'T']
#page configuration shows up on the top 
st.set_page_config(page_title="Dashboard",
                   page_icon=":atom_symbol:",
                   layout="wide"
)
#getting the file
uploaded_file = st.file_uploader("Upload your Excel file", type=["xlsx"])
#set the default data(If you have more computing power then you can show more rows)
df = pd.read_excel(
    io = r'C:\Users\Darsh\Downloads\Dashboard\Large_sheet.xlsx',
    engine='openpyxl',
    sheet_name='Sheet1',
    skiprows=0,
    usecols='A:S',
    nrows=1000,
)
# this is to upload your own xcel sheet( and again if you have more computing power then you can show more rows)
if uploaded_file:
    xls = pd.ExcelFile(uploaded_file, engine='openpyxl')
    sheet_names = xls.sheet_names
    selected_sheet = st.selectbox("Choose a sheet", sheet_names)
    df = pd.read_excel(
        xls,
        sheet_name=selected_sheet,
        usecols='A:S',
        nrows=1000,
    )

    st.success(f"Loaded '{selected_sheet}' successfully")
else:
    st.warning("Remember you can only upload .xlsx files")

#side bar code
all_residues = set("".join(df["UniqueAminoAcids"].dropna()))

# residues filter
st.sidebar.header("Filters")
residues = st.sidebar.multiselect(
    "Select Residues:",
    options=sorted(all_residues),
)

#species filter
Species = st.sidebar.multiselect(
    "Select Species Name:",
    options=df["SpeciesSciName"].unique(),
)
# name search/filter
FileName = st.sidebar.multiselect(
    "Select FileName:",
    options=df["FileName"].unique(),
)
# min chains filter
min_chains = st.sidebar.number_input(
    "Minimum Number of Chains:",
    min_value=int(df["NumOfChains"].min()),
    max_value=int(df["NumOfChains"].max()),
    value=int(df["NumOfChains"].min())
)
df["SequenceLength"] = df["FullSequence"].apply(lambda x: len(x) if pd.notnull(x) else 0)

# Sidebar slider for sequence length
min_length, max_length = st.sidebar.slider(
    "Select Sequence Length Range:",
    min_value=int(df["SequenceLength"].min()),
    max_value=int(df["SequenceLength"].max()),
    value=(int(df["SequenceLength"].min()), int(df["SequenceLength"].max())),
    step=1
)
filtered_df = df[
    (df["NumOfChains"] >= min_chains) &
    (df["SequenceLength"] >= min_length) &
    (df["SequenceLength"] <= max_length) &
    (df["SpeciesSciName"].isin(Species) if Species else True) &
    (df["FileName"].isin(FileName) if FileName else True) &
    (df["UniqueAminoAcids"].apply(lambda x: all(res in x for res in residues)) if residues else True) #&
]
# filtered data output
st.write("Filtered Data:", filtered_df)
st.write(f"Filtered Data Shape: {filtered_df.shape}")
total_acid_counts = getaacounts(filtered_df)
df_counts = pd.DataFrame({
    'Residue': amino_acids,
    'Frequency': list(total_acid_counts.values())
})
# display amino acid frequency throughout the data
fig = px.bar(
    df_counts,
    x='Residue',
    y='Frequency',
    title='Amino Acid Frequency from DataFrame Sequences',
    color_discrete_sequence=['blue']
)
fig.update_layout(
    xaxis_title='Residues',
    yaxis_title='Frequency Through All Data'
)
# Display 4mer, trimer, and dimer groups 
st.plotly_chart(fig, use_container_width=True) 
groups_4 = extract_helix_4mers(filtered_df)
groups_3 = extract_helix_trimers(filtered_df)
groups_2 = extract_helix_dimers(filtered_df)
top_100_4 = sorted(groups_4.items(), key=lambda x: x[1], reverse=True)[:100]
groups, counts = zip(*top_100_4)
fig2 = px.bar(
        x=groups,
        y=counts,
        color_discrete_sequence=['red'],
        title='Occurrence'
)
fig2.update_layout(
    xaxis_title='4-mer Group',
    yaxis_title='Occurrence in Helixes',
    title='Top 100 4-mer Group Frequency in Helixes',
)
st.plotly_chart(fig2, use_container_width=True)
top_100_3 = sorted(groups_3.items(), key=lambda x: x[1], reverse=True)[:100]
groups, counts = zip(*top_100_3)
fig3 = px.bar(
        x=groups,
        y=counts,
        color_discrete_sequence=['orange'],
        title='Occurrence'
)
fig3.update_layout(
    xaxis_title='Trimer Group',
    yaxis_title='Occurrence in Helixes',
    title='Top 100 trimer Group Frequency in Helixes',
)
st.plotly_chart(fig3, use_container_width=True)
top_100_2 = sorted(groups_2.items(), key=lambda x: x[1], reverse=True)[:100]
groups, counts = zip(*top_100_2)
fig4 = px.bar(
        x=groups,
        y=counts,
        color_discrete_sequence=['green'],
        title='Occurrence'
)
fig4.update_layout(
    xaxis_title='4-mer Group',
    yaxis_title='Occurrence in Helixes',
    title='Top 100 dimer Group Frequency in Helixes',
)
st.plotly_chart(fig4, use_container_width=True)
