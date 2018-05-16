import warnings
import numpy as np
from .utils import read_bmrb, extract_refcs, compute_error, write
warnings.simplefilter("always")

def checkcarbons(bmrb_id, fmt='nmrstar', cutoff=.5): 
    """
    Checks for systematic errors in the 13C chemical shifts of a bmrb entry
    if fmt = 'nmrstar', saves a nmrstar file with corrected 13C chemical shifts
    if fmt = 'csv', saves a csv file with corrected 13C chemical shifts
    if fmt = 'df', returns a DataFrame with corrected 13C chemical shifts
    Parameters
    ----------    
    bmrb_id  : int or str
        BMRB accession number or path to nmrstar file
    fmt : str
        'nmrstar','csv' or 'df'. default value = 'nmrstar'
    cutoff : float
        cutoff for error deviation. default value = 0.5
    """
    results = []
    dataframes, ref_sequences, bmrb_cont, correct_filename = read_bmrb(bmrb_id)   
    
    if fmt not in ['nmrstar', 'csv', 'df']:
        raise ValueError('The output format {}, is not supported.'.format(fmt))
    
    if type(bmrb_id) == int:
        legend = 'BMRB id'
        
    if type(bmrb_id) == str:
        legend = 'File'
    
    error = None
    error_nb = None
    error_r =None
    
    for idx, dataset in enumerate(dataframes):
        carbons = ["C1'", "C2'", "C3'", "C4'",
                   "C5'", "C2", "C4", "C8", "C5", "C6"]
        ref_seq = ref_sequences[idx]

        if ref_seq != 'GGC':
            #warnings.warn("No 5'-GG/3'-C terminal sequence")
            warnings.warn("The method cannot be applied because this RNA molecule has no 5'-GG/3'-C terminal sequence")

        else:
            merged_df, in_range = extract_refcs(dataset)

            # There are less than 2 reference 13C chemical shifts
            if len(merged_df['cs_val'].dropna()) < 2:
                print('{:} {:} has not enough '
                          'reference nuclei'.format(legend,bmrb_id))

            # There are at least 2 of the 5 reference 13C values
            else:
                # All reference values are inside the expected ranges
                if np.all(in_range):
                    print('Nitrogenous base 13C chemical shifts '
                          'of {:} {:} are correct\n'
                          'Ribose 13C chemical shifts of {:} '
                          '{:} are correct'.format(legend,bmrb_id,legend,bmrb_id))

                # All reference values are outside the expected ranges
                elif not np.any(in_range):
                    error = compute_error(merged_df,cutoff)

                    if error is not None:
                        results.append((bmrb_id, bmrb_cont, error,
                                        carbons, dataset, fmt))
                        print('Nitrogenous base 13C chemical shifts '
                              'of {:} {:} have a systematic error '
                              'of {:.2f}ppm\nRibose 13C chemical shifts '
                              'of {:} {:} have a systematic error '
                              'of {:.2f}ppm'.format(legend,bmrb_id,error,legend,bmrb_id,error))
                    else:
                        print('13C chemical shifts of {:} {:} have'
                              ' non-systematic errors'.format(legend,bmrb_id))

                # A portion of the reference chemical shift
                # are inside the expected ranges, another part is outside
                else:
                    nbase = [cs for cs in in_range[:3]]
                    ribose = [cs for cs in in_range[3:]]

                    # All reference values of nitrogenous bases are inside the
                    # expected ranges
                    if np.all(nbase):
                        nb = 'Nitrogenous base 13C chemical shifts of {:} {:} are correct'.format(legend,bmrb_id,legend)
                        
                    # All reference  values of ribose are outside the
                    # expected ranges
                    if not np.any(ribose):
                        error_r = compute_error(merged_df.iloc[3:],cutoff)

                        if error_r is not None:
                            results.append((bmrb_id, bmrb_cont, error_r,
                                            carbons[:5], dataset, fmt))
                           
                            print('{:}\nRibose 13C chemical shifts of {:} '
                                  '{:} have a systematic error of {:.2f} '
                                  'ppm'.format(nb,legend,bmrb_id, error_r))                  
                    
                    # All reference values of nitrogenous bases are
                    # outside the expected ranges
                    if not np.any(nbase):
                        error_nb = compute_error(merged_df.iloc[:3],cutoff)

                        if error_nb is not None:
                            results.append((bmrb_id, bmrb_cont, error_nb,
                                            carbons[:5], dataset, fmt))
                            nb = 'Nitrogenous base 13C chemical shifts of {:} {:} have a systematic error of {:.2f} ppm\n'.format(legend,bmrb_id,error_nb)
                    
                    # All reference chemical shift values of ribose are
                    # inside the expected ranges
                    if np.all(ribose):
                        print('{:}\nRibose 13C chemical shifts of [.] '
                                '{:} are correct'.format(nb,legend,bmrb_id))
                    
                    if error_r is None and error_nb is None:
                        print('13C chemical shifts of {:} {:} have non-'
                              'systematic errors'.format(legend,bmrb_id))

    if results:
        return write(results,correct_filename,fmt)
