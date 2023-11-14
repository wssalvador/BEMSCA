import math
import sqlite3
import pandas as pd

# -----------------------------------------------------------------------------
#    FUNCTIONS
# -----------------------------------------------------------------------------
# Function for calculating std from sem
def sem_to_std(sem, sample_size):
    std = sem * math.sqrt(sample_size)
    return std


# Function for calculating alpha and beta parameters of beta probability distribution
def alpha_beta(avg, std):
    precision = ((avg * (1 - avg)) / std**2) - 1
    alpha = avg * precision
    beta = (1 - avg) * precision
    return (alpha, beta)


# Function for accessing database data when BEMSCA starts up
def get_database_data():
    # Connect to database
    connection = sqlite3.connect("database.db")

    # Create cursor
    cursor = connection.cursor()

    # Get table names from database
    cursor.execute(f'SELECT Name FROM sqlite_master WHERE TYPE = "table"')
    tables = [name[0] for name in cursor.fetchall()]

    # Extract data for all tables
    table_dict = {}
    for table in tables:
        # Get column names from table
        cursor.execute(f'PRAGMA table_info("{table}")')
        column_names = [info[1] for info in cursor.fetchall()]

        # Get entries from table
        cursor.execute(f'SELECT * FROM "{table}"')
        entries = cursor.fetchall()

        # Insert entries into a pandas data frame to keep data organized
        index = []
        entry_dict = {}
        for i, column_name in enumerate(column_names):
            # Use 'Name' column as index of pandas dataframes
            if column_name == 'Name':
                index = [entry[i] for entry in entries]

            # Use all other columns as contents of pandas dataframes
            else:
                entry_dict[column_name] = [entry[i] for entry in entries]

        # Save data as pandas dataframe to dictionary with all dataframes
        table_dict[table] = pd.DataFrame(entry_dict, index=index)

    return table_dict