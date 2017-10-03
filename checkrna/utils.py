import decimal
import numpy as np
import pandas as pd
import pynmrstar as bmrb
bmrb.CONVERT_DATATYPES = True

# Python 2.x compatibility
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

def read_bmrb(bmrb_id, bmr_path='bmrbs'):
    """
    Reads nmrstar file of a given BMRB entry and returns its content

    Parameters
    ----------
    bmrb_id : int
        bmrb accession number
    bmr_path : string
        path to the folder with nmrstar files

    Returns
    -------
    dataframes : list
        list of DataFrames containing RNA 13C chemical shift data
    ref_sequences : list
       list of sequences used as reference
        '-' is used in a sequence when a reference residue doesn't exist
    bmrb_cont : 
        PyNMRSTAR object representing the nmrstar file
    """
    try:
        bmrb_cont = bmrb.Entry.from_file('{}/bmr{}.str'.format(bmr_path,
                                                               bmrb_id))
    except FileNotFoundError:
        bmrb_cont = bmrb.Entry.from_database(bmrb_id)

    acs = bmrb_cont.get_saveframes_by_category('assigned_chemical_shifts')

    dataframes = []
    ref_sequences = []
    carbons = ["C1'", "C2'", "C3'", "C4'", "C5'",
               "C2", "C4", "C8", "C5", "C6"]

    for d in acs:

        cs = bmrb_cont.get_saveframe_by_name(
            d.name).get_loop_by_category('Atom_chem_shift')

        data = {'id': cs['ID'],
                'entity_assembly_id': cs['Entity_assembly_ID'],
                'seq_id': [int(i) for i in cs['Seq_ID']],
                'comp_id': cs['Comp_ID'],
                'atom_id': cs['Atom_ID'],
                'cs_val': [float(i) for i in cs['Val']]}

        df = pd.DataFrame(data, columns=list(data.keys())[:])

        n_of_seqs = list(sorted(set(df['entity_assembly_id'])))

        for n_seq in n_of_seqs:

            df_seq = df[df['entity_assembly_id'] == n_seq]
            index_seq = df_seq.groupby('seq_id').first()['comp_id']

            if len(index_seq) >= 3:
                seq = df_seq.groupby('seq_id').first()['comp_id'].values
                ref_seq = ""

                if index_seq.index[0] == 1:
                    ref_seq = seq[0]
                else:
                    ref_seq = "-"

                if index_seq.index[1] == 2:
                    ref_seq += seq[1]
                else:
                    ref_seq += "-"

                ref_seq += seq[-1]

                df_seq = df_seq[df_seq['atom_id'].isin(carbons)]

                if not df_seq.empty:
                    ref_sequences.append(ref_seq)
                    dataframes.append(df_seq)

    return dataframes, ref_sequences, bmrb_cont


def extract_refcs(df):
    """
    Extract the CS of the five reference nuclei from df and combines them with
    the reference ranges DataFrame

    Parameters
    ----------    
    df : DataFrame
        DataFrame with all 13C chemical shifts (CS) of a RNA structure

    Returns
    ----------
    merged_df : DataFrame
        DataFrame with 'cs_val' from df merged with the reference ranges
    in_range : Series
        Series with boolean values, True/False if cs_val is inside/outside the
        range for the corresponding 13C chemical shift reference values

    """
    last = df['seq_id'].values[-1]

    reference = {'comp_id': ["G", "G", "C", "C", "C"],
                 'atom_id': ["C8", "C8", "C5", "C1'", "C3'"],
                 'seq_id': [1, 2, last, last, last],
                 'min_val': [138.7, 136.4, 97.4, 92.5, 69.4],
                 'max_val': [139.7, 137.6, 98.8, 93.4, 70.4],
                 'mean_val': [139.20, 137.00, 98.10, 92.95, 69.90]}

    df_ranges = pd.DataFrame.from_dict(reference)

    merged_df = pd.merge(df_ranges, df, on=['atom_id', 'comp_id', 'seq_id'],
                         how='left')

    in_range = (merged_df['cs_val'] >= merged_df['min_val']) & (
                merged_df['cs_val'] <= merged_df['max_val'])

    return merged_df, in_range


def compute_error(dataset, merged_df):
    """
    Searches for a systematic error in 13C chemical shifts based on deviation 
    from expected values and returns the error if found

    Parameters
    ----------    
    dataset : DataFrame
        DataFrame with a subset of 13C chemical shifts values
    merged_df : DataFrame
        DataFrame with 'cs_val' from df merged with the reference ranges

    Returns
    ----------
    error : float
        systematic error computed as the mean of the differences between cs_val
        and mean_val
    """
    deltas = merged_df['cs_val'] - merged_df['mean_val']
    deviation = np.nanstd(deltas)
    error = np.nanmean(deltas)

    if deviation <= 1.:
        return error


def write(results):
    """
    Saves or returns bmrb_cont with corrected 13C chemical shifts

    Parameters
    ----------
    Results : tuple
        bmrb_cont : PyNMRSTAR Entry
            object represeting the nmrstar file content
        error : float
            systematic error found for 13C chemical shifts
        carbons : list
            list of carbon nuclei used for selecting the subset of chemical
            shifts over which correction is computed
        dataset : DataFrame
            DataFrame with a subset of 13C chemical shifts values
        fmt : str
            Available option 'nmrstar','csv' or 'df' (default 'nmrstar')

    Returns
    ----------
    Corrected 13C chemical shifts. 
    3 options are available:
        * 'nmrstar', saves a nmrstar file with corrected 13C chemical shifts
        * 'csv', saves a csv file with corrected 13C chemical shifts
        * 'df', returns a DataFrame with corrected 13C chemical shifts    
    """
    correct_datasets = []
    for res in results:
        bmrb_id, bmrb_cont, error, carbons, dataset, fmt = res

        n_dataset = dataset['entity_assembly_id'].unique()

        if fmt == 'nmrstar':
            tag = ['Entity_ID', 'Atom_ID', 'Seq_ID', 'Comp_ID', 'Val',
                   'Entity_assembly_ID']

            for cs_sf in bmrb_cont.get_saveframes_by_category('assigned_chemical_shifts'):
                shift_loop = cs_sf.get_loop_by_category('atom_chem_shift')
                cs_index = shift_loop.columns.index('Val')

                for rownum, rowdata in enumerate(shift_loop.get_tag(tag)):
                    for carbon in carbons:
                        # If the nuclei is an RNA carbon
                        if rowdata[1] == carbon and rowdata[-1] == n_dataset:
                            error = decimal.Decimal(error)
                            #error = round(decimal.Decimal(error), 3)
                            shift_loop.data[rownum][cs_index] -= error

        elif fmt == 'csv' or fmt == 'df':

            for carbon in carbons:
                c = dataset.ix[(dataset['atom_id'] == carbon),
                               'cs_val'] - error
                dataset.ix[(dataset['atom_id'] == carbon), 'cs_val'] = c

            correct_datasets.append(dataset)

    if fmt == 'nmrstar':
        fb = open('bmr{}_correct.str'.format(bmrb_id), 'w')
        fb.write(str(bmrb_cont))
        fb.close()

    elif fmt == 'csv':
        final_dataset = pd.concat(correct_datasets)[['id',
                                                     'entity_assembly_id',
                                                     'comp_id',
                                                     'seq_id',
                                                     'atom_id',
                                                     'cs_val']]

        final_dataset.to_csv('bmr{}_correct.csv'.format(bmrb_id))

    elif fmt == 'df':
        final_dataset = pd.concat(correct_datasets)[['id',
                                                     'entity_assembly_id',
                                                     'comp_id',
                                                     'seq_id',
                                                     'atom_id',
                                                     'cs_val']]

        return(final_dataset)
