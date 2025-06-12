import streamlit as st
import pandas as pd
from utils.lipinski import detect_smiles_column, analyze_lipinski
import io

# Configure page
st.set_page_config(page_title="Lipinski's Rule of Five Analyzer", layout="wide")

# Title and description
st.title("Lipinski's Rule of Five Analyzer")
st.markdown("""
This tool analyzes chemical compounds based on Lipinski's Rule of Five using RDKit.
Upload a CSV file containing SMILES strings to analyze drug-likeness.
""")

# File upload
uploaded_file = st.file_uploader("Upload CSV file", type=['csv'])

if uploaded_file is not None:
    try:
        # Read CSV file
        df = pd.read_csv(uploaded_file)
        
        # Detect SMILES column
        try:
            default_smiles_col = detect_smiles_column(df)
            smiles_col = st.selectbox(
                "Select the SMILES column", 
                options=df.columns,
                index=df.columns.get_loc(default_smiles_col)
            )
        except ValueError as e:
            st.error(str(e))
            st.stop()
        except Exception as e:
            st.error(f"An error occurred while detecting SMILES column: {str(e)}")
            st.stop()
        
        # Analyze data
        with st.spinner('Analyzing compounds...'):
            try:
                result_df = analyze_lipinski(df, smiles_col)
                
                # Display statistics
                valid_smiles = result_df['SMILES_Valid'].sum()
                total_compounds = len(result_df)
                pass_count = (result_df['LipinskiResult'] == 'Pass').sum()
                fail_count = (result_df['LipinskiResult'] == 'Fail').sum()
                
                st.success(f"Analysis complete! Processed {valid_smiles} valid SMILES out of {total_compounds} compounds.")
                
                # Show summary stats
                col1, col2, col3 = st.columns(3)
                col1.metric("Total Compounds", total_compounds)
                col2.metric("Valid SMILES", valid_smiles, f"{valid_smiles/total_compounds:.1%}")
                col3.metric("Pass Rate", f"{pass_count}/{valid_smiles}", f"{pass_count/valid_smiles:.1%}" if valid_smiles else "N/A")
                
                # Visualization
                st.subheader("Pass/Fail Distribution")
                if valid_smiles > 0:
                    chart_data = pd.DataFrame({
                        'Result': ['Pass', 'Fail'],
                        'Count': [pass_count, fail_count]
                    })
                    st.bar_chart(chart_data.set_index('Result'))
                else:
                    st.warning("No valid SMILES to display results.")
                
                # Show data table
                st.subheader("Results Table")
                st.dataframe(result_df)
                
                # Download button
                output = io.StringIO()
                result_df.to_csv(output, index=False)
                output.seek(0)
                
                st.download_button(
                    label="Download Results as CSV",
                    data=output,
                    file_name='lipinski_results.csv',
                    mime='text/csv'
                )
                
            except Exception as e:
                st.error(f"An error occurred during analysis: {str(e)}")
                
    except Exception as e:
        st.error(f"Error reading CSV file: {str(e)}")
else:
    st.info("Please upload a CSV file to begin analysis.")
