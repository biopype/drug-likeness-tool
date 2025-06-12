import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from typing import Optional
import re

def detect_smiles_column(df: pd.DataFrame) -> str:
    """
    Detect the SMILES column in a DataFrame by checking common naming variations.
    
    Args:
        df: Input DataFrame containing chemical compounds data
        
    Returns:
        str: Name of the detected SMILES column
        
    Raises:
        ValueError: If no valid SMILES column is found
    """
    # Common SMILES column name variations (case-insensitive)
    possible_names = ['smiles', 'smi', 'smile', 'structure']
    
    for col in df.columns:
        # Remove any non-alphanumeric characters and make lowercase
        col_clean = re.sub(r'[^a-zA-Z0-9]', '', col).lower()
        if col_clean in possible_names:
            return col
    
    # If no column found with exact match, try partial matches
    for col in df.columns:
        col_lower = col.lower()
        if any(name in col_lower for name in possible_names):
            return col
    
    raise ValueError("No valid SMILES column found. Please ensure your CSV contains a column with SMILES strings.")

def analyze_lipinski(df: pd.DataFrame, smiles_col: str) -> pd.DataFrame:
    """
    Analyze compounds using Lipinski's Rule of Five and add computed properties to DataFrame.
    
    Args:
        df: Input DataFrame containing chemical compounds
        smiles_col: Name of the column containing SMILES strings
        
    Returns:
        pd.DataFrame: Original DataFrame with added computed columns
    """
    # Create a copy to avoid modifying the original DataFrame
    result_df = df.copy()
    
    # Initialize new columns
    result_df['MolWt'] = None
    result_df['LogP'] = None
    result_df['NumHDonors'] = None
    result_df['NumHAcceptors'] = None
    result_df['LipinskiViolations'] = 0
    result_df['LipinskiResult'] = "Invalid SMILES"
    result_df['SMILES_Valid'] = False
    
    for idx, row in df.iterrows():
        smiles = row[smiles_col]
        
        if not isinstance(smiles, str) or not smiles.strip():
            continue
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
            
        # Mark as valid SMILES
        result_df.at[idx, 'SMILES_Valid'] = True
        
        # Calculate properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        h_donors = Lipinski.NumHDonors(mol)
        h_acceptors = Lipinski.NumHAcceptors(mol)
        
        # Count violations of Lipinski's rule (max 4 violations)
        violations = 0
        if mw > 500:
            violations += 1
        if logp > 5:
            violations += 1
        if h_donors > 5:
            violations += 1
        if h_acceptors > 10:
            violations += 1
            
        # Update the row with calculated values
        result_df.at[idx, 'MolWt'] = mw
        result_df.at[idx, 'LogP'] = logp
        result_df.at[idx, 'NumHDonors'] = h_donors
        result_df.at[idx, 'NumHAcceptors'] = h_acceptors
        result_df.at[idx, 'LipinskiViolations'] = violations
        result_df.at[idx, 'LipinskiResult'] = "Pass" if violations == 0 else "Fail"
    
    return result_df