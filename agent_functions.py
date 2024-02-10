import pandas as pd
from chembl_webresource_client.new_client import new_client
from mordred import Calculator, descriptors
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from typing import List

from typing_extensions import Annotated

import autogen
from autogen.cache import Cache

config_list = autogen.config_list_from_json(env_or_file="OAI_CONFIG_LIST")

def chat_friendly_df_summary(func):
    # TODO: annnotations are needed even for decorators
    def wrapper(*args: str, **kwargs: str):
        # Call the function and get the DataFrame output
        obj = func(*args, **kwargs)

        if isinstance(obj, pd.DataFrame):
            # Convert DataFrame to dictionary
            result_dict = {
                'columns': obj.columns.tolist(),
                'num_rows': obj.shape[0]
            }
        elif isinstance(obj, pd.Series):
            # Convert Series to dictionary
            result_dict = {
                'index': obj.index.tolist(),
                'values': obj.values.tolist()
            }
        
        return result_dict
    
    return wrapper

def summarize_df(obj: pd.DataFrame, file_name: str=None) -> dict:
    if isinstance(obj, pd.DataFrame):
        # Convert DataFrame to dictionary
        result_dict = {
            'columns': obj.columns.tolist(),
            'num_rows': obj.shape[0]
        }
    elif isinstance(obj, pd.Series):
        # Convert Series to dictionary
        result_dict = {
            'index': obj.index.tolist(),
            'values': obj.values.tolist()
        }
    
    if file_name:
        result_dict['csv_output'] = file_name
    return result_dict


def target_query_results_file(target_name: str) -> str:
    """
    Returns the filename for the target query results.

    Args:
        target_name (str): The name of the target.

    Returns:
        str: The filename for the target query results.
    """
    return f"target_query_results_{target_name}.csv"

#@chat_friendly_df_summary
async def download_protein_results(target_query_string: str) -> dict:
    """ 
    Download the protein results from ChEMBL and save to a csv file
    
    Args:
        target_query_string: str - the query string used to search for protein targets in ChEMBL

    Returns:
        pd.DataFrame - proteins returned from the query to the ChEMBL database
    """
    target = new_client.target
    target_result = target.search(target_query_string)
    
    # Save results to a table
    df = pd.DataFrame(target_result)
    file_name = target_query_results_file(target_query_string)
    df.to_csv(file_name, index=False)

    return summarize_df(df, )

# @chat_friendly_df_summary
async def select_target_from_query_results(chembl_id: str, target_query_string: str) -> dict:
    """ Select a target from the ChEMBL protein results csv file

    Args:
        chembl_id (str): The ChEMBL ID of the chosen target protein
        target_id (str): The query string used to pull results from ChEMBL

    Returns:
        pd.Series: The selected target protein as a pandas Series

    Raises:
        ValueError: If the specified chembl_id is not found in the protein results dataframe
    """
    # Load the protein results from the csv file
    protein_result_df = pd.read_csv(target_query_results_file(target_query_string))

    # returns an error message if the id does not appear in the table  
    if not protein_result_df.target_chembl_id.eq(chembl_id).any():
        raise ValueError(f"Error: {chembl_id} not found in the protein results dataframe")
    else:
        return summarize_df(protein_result_df[protein_result_df['target_chembl_id'] == chembl_id])

#@chat_friendly_df_summary
async def generate_activity_data(chembl_id:str, standard_type:str='IC50') -> pd.DataFrame:
    """
    Generate a pandas DataFrame containing activity data for a given ChEMBL target id and save it to a csv.

    Parameters:
    chembl_id (str): The ChEMBL ID of a Protein Target.
    standard_type (str, optional): The standard type of the activity data. Defaults to IC50. Standard types include: IC50, EC50, Ki, Kd, and Potency.

    Returns:
    pd.DataFrame: The generated DataFrame containing compound activity data for the given protein target.

    Raises:
    ValueError: If no activity data is found for the given target ID.
    ValueError: If activity data is found for the target ID, but not with the specified standard type.
    """
    
    activity = new_client.activity
    activity_result = activity.filter(target_chembl_id=chembl_id)

    if not activity_result:
        raise ValueError(f"No activity data found for {chembl_id}. Please try another target.")

    activity_result = activity_result.filter(standard_type=standard_type)
        
    if not activity_result:
        raise ValueError(f"Activity data was found for this target, but not with the standard type '{standard_type}'. Please try another standard type.")
    # Save results to a dataframe
    df = pd.DataFrame.from_dict(activity_result)

    # drop rows with no smiles
    num_rows = df.shape[0]
    df.dropna(subset=['canonical_smiles'], inplace=True)
    num_rows_dropped = num_rows - df.shape[0]
    print("Dropped ", num_rows_dropped, " rows with no smiles.")

    # save the dataframe to a csv file
    file_name = f'activity_data_{chembl_id}_{standard_type}.csv'
    df.to_csv(file_name, index=False)
    return summarize_df(df, file_name)

# @chat_friendly_df_summary
async def calculate_lipinski_descriptors(activity_file: str) -> pd.DataFrame:
    """
    Calculate Lipinski descriptors for compounds based on their SMILES strings.

    Parameters:
    activity_file (str): Path to the activity file containing SMILES strings.

    Returns:
    pd.DataFrame: DataFrame containing the calculated Lipinski descriptors.
    """

    # Get the smiles strings for the compounds
    activity_data = pd.read_csv(activity_file)

    # get a list of the smiles strings from the activity data dataframe
    smiles_series = list(activity_data['canonical_smiles'])

    # Calculate the Lipinski descriptors for each compound
    lipinski_descriptors = []
    for smiles in smiles_series:
        mol = Chem.MolFromSmiles(smiles)
        lipinski_descriptors.append([
            Descriptors.MolWt(mol),
            Lipinski.NumHDonors(mol),
            Lipinski.NumHAcceptors(mol),
            Descriptors.NumRotatableBonds(mol),
            Descriptors.MolLogP(mol)
        ])

    # Create a DataFrame with the Lipinski descriptors
    lipinski_df = pd.DataFrame(lipinski_descriptors, columns=['Molecular Weight', 'H-Bond Donors', 'H-Bond Acceptors', 'Rotatable Bonds', 'LogP'])

    # save the dataframe to a csv file

    # get the target id and actvity type from the activity file name
    target_id = activity_file.split('_')[2]
    activity_type = activity_file.split('_')[3].split('.')[0]
    file_name = f'{target_id}_{activity_type}_lipinski.csv'
    lipinski_df.to_csv(file_name, index=False)

    # Print the Lipinski descriptors
    return summarize_df(lipinski_df, file_name)


def calculate_mordred_descriptors(smiles: List[str]) -> pd.DataFrame:
    calc = Calculator(descriptors, ignore_3D=True)

    mordred_descriptors = calc.pandas([Chem.MolFromSmiles(s) for s in smiles])

    mordred_descriptors.to_csv('mordred_descriptors.csv', index=False)

    return summarize_df(mordred_descriptors)