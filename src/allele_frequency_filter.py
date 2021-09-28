import logging


def filter1(tally_dict, num_bison, num_cattle):
    """
    This function will filter the data based on...
            1. More than 0 homozygous alternative cows -> could indicate common ancestry
            2. More than 0 heterozygous cows
            3. Heterozygous bison frequency less than 30% of total called bison.
    @param tally_dict: Dictionary containing the records tally of genotypic frequencies of the samples.
    @return: True if the tally dict passes all checks, False if it does not pass all checks.
    """
    try:
        # Initial REJECT if less than 30% of cows are called. ***
        assert(tally_dict['total_cows_called'] / num_cattle > 0.3)
        # Fix ZeroDivisionError
        assert(tally_dict["total_cows_called"] != 0)
        # No homozygous alternative cattle samples
        assert(tally_dict['cow_hom_alt'] == 0)
        # No heterozygous cattle samples
        assert(tally_dict['cow_het'] == 0)
        # Only one heterozygous bison sample
        assert(tally_dict['bison_het'] == 1)
        # No homozygous reference bison
        assert(tally_dict['bison_hom_ref'] == 0)
        # Initial REJECT if less than 30% of bison are called. ***
        assert(tally_dict['total_bison_called'] / num_bison > 0.3)
        # Ensure there are homozygous alternative bison present
        assert(tally_dict['bison_hom_alt'] != 0)
        return True
    except AssertionError:
        return False


def sample_tally(hom_refs, hom_alts, hets, samp_ref_dict, num_bison, num_cattle):
    """
    This function counts the amount of homozygous reference, heterozygous
    and homozygous alternative individuals of
    a particular species. (i.e. bison, cow, water buffalo)

    :param hom_refs: List of sample names that are homozygous for the reference allele
    :param hom_alts: List of sample names that are homozygous for the alternative allele
    :param hets: List of sample names that are heterozygous
    :param samp_ref_dict: Tsv file of three rows -- 'sampleID   0   PG(1-n)' -- Cross between Bison and Cow ;
    0=cow
    1=bison
    2=water buffalo
    @return: sample tally dictionary or None
    """

    total_cows_called = 0
    total_bison_called = 0
    total_wb = 0

    cow_hom_ref_tally = 0
    bison_hom_ref_tally = 0
    wb_hom_ref_tally = 0

    cow_hom_alt_tally = 0
    bison_hom_alt_tally = 0
    wb_hom_alt_tally = 0

    cow_het_tally = 0
    bison_het_tally = 0
    wb_het_tally = 0
    bison_het_pg1 = []

    # Handle homozygous reference samples
    for hom_ref in hom_refs:
        if len(hom_refs) == 0:
            return None
        else:
            if samp_ref_dict[hom_ref][0] == 0:
                total_cows_called += 1
                cow_hom_ref_tally += 1
                continue
            elif samp_ref_dict[hom_ref][0] == 1:
                total_bison_called += 1
                bison_hom_ref_tally += 1
                continue
            elif samp_ref_dict[hom_ref][0] == 2:
                total_wb += 1
                wb_hom_ref_tally += 1
                continue

    # Handle homozygous alternative samples
    for hom_alt in hom_alts:
        if len(hom_alts) == 0:
            return None
        else:
            if samp_ref_dict[hom_alt][0] == 0:
                total_cows_called += 1
                cow_hom_alt_tally += 1
                continue
            elif samp_ref_dict[hom_alt][0] == 1:
                total_bison_called += 1
                bison_hom_alt_tally += 1
                continue
            elif samp_ref_dict[hom_alt][0] == 2:
                total_wb += 1
                wb_hom_alt_tally += 1
                continue

    # Handle heterozygous samples
    for het in hets:
        if len(hets) == 0:
            return None
        else:
            if samp_ref_dict[het][0] == 0:
                total_cows_called += 1
                cow_het_tally += 1
                continue
            elif samp_ref_dict[het][0] == 1:
                total_bison_called += 1
                bison_het_tally += 1
                bison_het_pg1.append(samp_ref_dict[het][1])
                continue
            elif samp_ref_dict[het][0] == 2:
                total_wb += 1
                wb_het_tally += 1
                continue

    tally_dictionary = {
        'cow_hom_ref': cow_hom_ref_tally,
        'bison_hom_ref': bison_hom_ref_tally,
        'wb_hom_ref': wb_hom_ref_tally,
        'cow_hom_alt': cow_hom_alt_tally,
        'bison_hom_alt': bison_hom_alt_tally,
        'wb_hom_alt': wb_hom_alt_tally,
        'cow_het': cow_het_tally,
        'bison_het': bison_het_tally,
        'wb_het': wb_het_tally,
        'total_cows_called': total_cows_called,
        'total_bison_called': total_bison_called,
        'total_wb': total_wb,
        'pg1': bison_het_pg1,
    }

    # Run filter_1 and return
    filter_1 = filter1(tally_dictionary, num_bison, num_cattle)
    if filter_1 is False:
        return None
    elif filter_1 is True:
        return tally_dictionary

######################### Main Function Call ##############################
def allele_frequency_filter(record, samp_ref_df, num_bison, num_cattle):
    try:
        # Initial rejection if call rate less than 50%
        if record.call_rate < 0.5:
            return None

        else:
            # Collect sample names per genotype
            hom_alts = record.get_hom_alts()
            hom_refs = record.get_hom_refs()
            hets = record.get_hets()

            # Put sample names per genotype into lists
            hom_alts_list = [str(data_alt.sample) for data_alt in hom_alts]
            hom_ref_list = [str(data_ref.sample) for data_ref in hom_refs]
            hets_list = [str(data_het.sample) for data_het in hets]

            # tally up the amount of het bison-cows, homozygous reference bison-cows,...etc.
            # Returns False is record not introgressed, True if it is.
            st = sample_tally(
                hom_ref_list,
                hom_alts_list,
                hets_list,
                samp_ref_df,
                num_bison,
                num_cattle,
            )

            if not st:
                return None
            else:
                output_list = [
                    str(record.CHROM),
                    int(record.POS),
                    str(record.REF),
                    str(record.ALT).strip("[]"),  # We are assuming only bi-allelic positions
                    float(round((1 - record.aaf[0]), 4)),
                    float(round(record.aaf[0], 4)),
                    float(round(record.call_rate, 4)),
                    st['cow_hom_ref'],
                    st['bison_hom_ref'],
                    st['wb_hom_ref'],
                    st['cow_hom_alt'],
                    st['bison_hom_alt'],
                    st['wb_hom_alt'],
                    st['cow_het'],
                    st['bison_het'],
                    st['wb_het'],
                    st['total_cows_called'],
                    st['total_bison_called'],
                    hets_list[0],
                    st['pg1'],
                ]
                return output_list
    except ZeroDivisionError:
        logging.debug(f'ZeroDivisionError - Skip record')
        return None
