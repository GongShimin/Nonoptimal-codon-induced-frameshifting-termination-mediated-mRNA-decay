import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import math

# Define the exponential decay fitting function
def exponential_decay(t, a, k):
    return a * np.exp(-k * t)

# Read data from a file
def read_data(file_path):
    try:
        # Assuming the txt file is tab-separated or space-separated
        return pd.read_csv(file_path, delimiter='\t', index_col=0)
    except FileNotFoundError:
        print(f"File {file_path} not found.")
        return None

# Fit the exponential decay model and calculate the half-life
def fit_data(data):
    times = np.array([0, 4, 8, 12])  # Time points in hours
    genes_half_life = {}
    
    for gene in data.index:
        expression_levels = data.loc[gene].dropna().values
        if len(expression_levels) < 2:  # Ensure there are at least two data points for fitting
            genes_half_life[gene] = float('inf')
            continue
        
        # Fit the exponential decay model
        popt, _ = curve_fit(exponential_decay, times, expression_levels, p0=[1, 0.1])
        a, k = popt
        half_life = math.log(2) / k if k != 0 else float('inf')
        genes_half_life[gene] = half_life

    return genes_half_life

def main():
    file_path = 'data.txt'  # Assume the input file is in txt format
    data = read_data(file_path)
    if data is not None:
        half_lives = fit_data(data)
        with open('half_life_results.txt', 'w') as fw:
            for gene, half_life in half_lives.items():
                fw.write(f"{gene}\t{half_life}\n")

if __name__ == "__main__":
    main()