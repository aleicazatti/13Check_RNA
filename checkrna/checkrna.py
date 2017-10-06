import warnings
import numpy as np
from .utils import read_bmrb, extract_refcs, compute_error, write


def checkcarbons(bmrb_id, fmt='nmrstar'):
    """
    Checks for systematic errors in the 13C chemical shifts of a bmrb entry

    if fmt = 'nmrstar', saves a nmrstar file with corrected 13C chemical shifts
    if fmt = 'csv', saves a csv file with corrected 13C chemical shifts
    if fmt = 'df', returns a DataFrame with corrected 13C chemical shifts

    Parameters
    ----------    
    bmrb_id  : int 
        accesion number of the bmrb entry to check for systematic errors

    fmt : str
        'nmrstar','csv' or 'df'. default value = 'nmrstar'

    """
    results = []

    dataframes, ref_sequences, bmrb_cont = read_bmrb(bmrb_id)

    if fmt not in ['nmrstar', 'csv', 'df']:
        raise ValueError('The output format {}, is not supported.'.format(fmt))

    for idx, dataset in enumerate(dataframes):
        carbons = ["C1'", "C2'", "C3'", "C4'",
                   "C5'", "C2", "C4", "C8", "C5", "C6"]
        ref_seq = ref_sequences[idx]

        if ref_seq != 'GGC':
            warnings.warn("No 5'-GG/3'-C terminal sequence")

        else:
            merged_df, in_range = extract_refcs(dataset)

            # There are less than 2 reference 13C chemical shifts
            if len(merged_df['cs_val'].dropna()) < 2:
                print('bmrb id {:d} has not enough '
                      'reference nuclei, {}'.format(bmrb_id, ref_seq))

            # There are at least 2 of the 5 reference 13C values
            else:
                # All reference values are inside the expected ranges
                if np.all(in_range):
                    print('13C chemical shifts of bmrb id'
                          ' {:d} are correct, {}'.format(bmrb_id, ref_seq))

                # All reference values are outside the expected ranges
                elif not np.any(in_range):
                    error = compute_error(merged_df)

                    if error is not None:
                        results.append((bmrb_id, bmrb_cont, error,
                                        carbons, dataset, fmt))
                        print('13C chemical shifts of bmrb id {:d} have a '
                              'systematic error of '
                              '{:.2f} ppm, {}'.format(bmrb_id, error, ref_seq))
                    else:
                        print('13C chemical shifts of bmrb id {:d} have'
                              ' non-systematic errors, {}'.format(bmrb_id,
                                                                  ref_seq))

                # A portion of the reference chemical shift
                # are inside the expected ranges, another part is outside
                else:
                    nbase = [cs for cs in in_range[:3]]
                    ribose = [cs for cs in in_range[3:]]

                    # All reference values of nitrogenous bases are inside the
                    # expected ranges
                    # All reference  values of ribose are outside the
                    # expected ranges
                    if np.all(nbase) and not np.any(ribose):
                        error = compute_error(merged_df.iloc[3:])

                        if error is not None:
                            results.append((bmrb_id, bmrb_cont, error,
                                            carbons[:5], dataset, fmt))
                            print('Ribose 13C chemical shifts of bmrb id '
                                  '{:d} have a systematic error of {:.2f} '
                                  'ppm, {}'.format(bmrb_id, error, ref_seq))

                    # All reference chemical shift values of ribose are
                    # inside the expected ranges
                    # All reference values of nitrogenous bases are
                    # inside the expected ranges
                    elif not np.any(nbase) and np.all(ribose):
                        error = compute_error(merged_df.iloc[:3])

                        if error is not None:
                            results.append((bmrb_id, bmrb_cont, error,
                                            carbons[:5], dataset, fmt))
                            print('Nitrogenous base 13C chemical shifts of '
                                  'bmrb id {:d} have a systematic error of '
                                  '{:.2f} ppm, {}'.format(bmrb_id, error,
                                                          ref_seq))
                    else:
                        print('13C chemical shifts of bmrb id {:d} have non '
                              'systematic errors, {}'.format(bmrb_id,
                                                             ref_seq))

    if results:
        return write(results)
