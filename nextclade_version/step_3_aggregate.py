import pandas as pd
import argparse

def process_data(file_path):
    # Load the TSV file into a DataFrame
    df = pd.read_csv(file_path, sep='\t', header=None)

    # Assign column names based on our previous analysis
    df.columns = ["Sequence", "Day", "Month", "Unknown", "Date", "Country", "EPI"]

    # Extract the year from the 'Date' column
    df['Year'] = df['Date'].str.extract('(\d{4})')

    # Aggregate and count by 'Country' only
    count_by_country = df.groupby('Country').size().reset_index(name='Count').sort_values(by='Count', ascending=False)

    # Aggregate and count by 'Year' only
    count_by_year = df.groupby('Year').size().reset_index(name='Count').sort_values(by='Year', ascending=False)

    # Print the results
    print("Counts by Country (Top 10):")
    print(count_by_country.head(10))
    print("\nCounts by Year:")
    print(count_by_year)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a TSV file to get aggregated counts by country and year.")
    parser.add_argument("file_path", type=str, help="Path to the TSV file.")

    args = parser.parse_args()
    process_data(args.file_path)
