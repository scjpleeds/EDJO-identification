import numpy as np
import pandas as pd

# Function to efficiently process regions using vectorized operations

def process_regions(data):
    # Assuming `data` is a DataFrame with regional data
    results = data.copy()
    
    # Example vectorized operation: calculating the mean and standard deviation
    results['mean'] = data.mean(axis=1)
    results['std_dev'] = data.std(axis=1)
    
    # More efficiency improvements can be added as needed
    return results

# Example usage
if __name__ == '__main__':
    # Creating a sample DataFrame as an example
    num_regions = 1000
    num_samples = 50
    sample_data = np.random.rand(num_regions, num_samples)
    df = pd.DataFrame(sample_data, columns=[f'Sample_{i}' for i in range(num_samples)])
    
    # Process the sample data
    processed_results = process_regions(df)
    print(processed_results.head())