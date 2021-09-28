"""
Author: Andrew Harris
Python 3.8.2
This script cross-checks the introgressed positions from filter 1 and checks
to see if any cattle show the alternative allele. If so, the introgressed postion
is removed.
"""
import pandas as pd

def vcf_check(
        introgressed_position_files,
        cattle_vcf_file,
        introgressed_output_dir,
        dropped_snp_output_dir,
        vcf_chrom,
        logger,
):
    for f in introgressed_position_files:
        output_filename = introgressed_output_dir / f"{f.stem}.vcfCheck.tsv"
        df = pd.read_csv(f, sep="\t")
        # print(df)
        output_df = df.copy()
        before_drop_len = len(output_df)

        drop_info_output_filename = dropped_snp_output_dir / f"{f.stem}.droppedSNPInfo.tsv"
        drop_info_df = pd.DataFrame(columns=['ReasonDropped', 'Chromosome', 'Position', 'BisonRefBase', 'BisonAltBase', 'CattleRefBase', 'CattleAltBase', 'NumHetCattle', 'HetCattle'])
        drop_idx = 0
        outputFile = True
        for row in df.itertuples():
            # Introgressed record information
            index = row.Index
            chromosome = row.Chromosome
            position = row.Position
            ref_base = row.Reference_Base
            alt_base = row.Alternative_Base.strip("[]")

            try:
                # Fetch corresponding record in cattle VCF file
                records = [record for record in cattle_vcf_file.fetch(str(chromosome), (int(position) - 1), (int(position)))]
            except ValueError:
                outputFile = False
                break 

            # Find position in cattle vcf
            for cattle_record in records:
                cattle_position = cattle_record.POS
                cattle_ref = cattle_record.alleles[0]
                cattle_alt = cattle_record.alleles[1:]
                hetCattle = [s.sample for s in cattle_record.get_hets()]

                # Ensure same positions
                if position != cattle_position:
                    continue
                # Ensure ref bases are same in both files
                elif ref_base != cattle_ref:
                    continue
                # Ensure only using biallelic records
                elif len(cattle_alt) != 1:
                    continue
                # Ensure introgressed alt == cattle_alt
                elif alt_base != cattle_alt[0]:
                    continue
                # If call rate is below 50% skip
                elif cattle_record.call_rate < 0.5:
                    drop_info_df.at[drop_idx, 'ReasonDropped'] = 'LowCallRate'
                    drop_info_df.at[drop_idx, 'Chromosome'] = chromosome
                    drop_info_df.at[drop_idx, 'Position'] = position
                    drop_info_df.at[drop_idx, 'BisonRefBase'] = ref_base
                    drop_info_df.at[drop_idx, 'BisonAltBase'] = alt_base
                    drop_info_df.at[drop_idx, 'CattleRefBase'] = cattle_ref
                    drop_info_df.at[drop_idx, 'CattleAltBase'] = cattle_alt
                    drop_info_df.at[drop_idx, 'NumHetCattle'] = len(hetCattle)
                    if len(hetCattle) > 0:
                        drop_info_df.at[drop_idx, 'HetCattle'] = ':'.join(hetCattle)
                    else:
                        drop_info_df.at[drop_idx, 'HetCattle'] = "NoHetCattle"
                    drop_idx += 1
                    output_df.drop(index=index, inplace=True)
                    break
                # All cattle should be homozygous for reference
                elif (1 - cattle_record.aaf[0]) != 1:
                    drop_info_df.at[drop_idx, 'ReasonDropped'] = 'AltCattlePresent'
                    drop_info_df.at[drop_idx, 'Chromosome'] = chromosome
                    drop_info_df.at[drop_idx, 'Position'] = chromosome
                    drop_info_df.at[drop_idx, 'BisonRefBase'] = ref_base
                    drop_info_df.at[drop_idx, 'BisonAltBase'] = alt_base
                    drop_info_df.at[drop_idx, 'CattleRefBase'] = cattle_ref
                    drop_info_df.at[drop_idx, 'CattleAltBase'] = cattle_alt
                    drop_info_df.at[drop_idx, 'NumHetCattle'] = len(hetCattle)
                    if len(hetCattle) > 0:
                        drop_info_df.at[drop_idx, 'HetCattle'] = ':'.join(hetCattle)
                    else:
                        drop_info_df.at[drop_idx, 'HetCattle'] = "NoHetCattle"
                    drop_idx += 1
                    output_df.drop(index=index, inplace=True)
                    break
                # If not rejected by any criteria, the introgressed position is kept.
                else:
                    break
            continue
        if outputFile:
            after_len = len(output_df)
            logger.info(f"File: {f.name}\nPre-filter total: {before_drop_len:,}\nTotal positions dropped: {(before_drop_len - after_len):,}\nPositions remaining: {(before_drop_len - (before_drop_len - after_len)):,}")
            # Write output_df
            output_df.to_csv(output_filename, sep="\t", index=False)
            drop_info_df.to_csv(drop_info_output_filename, sep='\t', index=False)
        else:
            continue
    return


