"""
SRA metadata extraction from NCBI databases.

This module extracts sample and sequencing metadata for SRA runs using a hybrid
approach that combines pre-built NCBI metadata tables with batch API queries.

Workflow:
    1. Extract SRA IDs from viridian metadata
    2. Map SRA IDs to BioSample IDs from NCBI SRA_Accessions.tab (streaming)
    3. Batch query BioSample database for sample metadata (collection date, location, etc.)
    4. Batch query SRA database for experiment design descriptions
    5. Merge all metadata into final comprehensive TSV

"""

rule get_metadata_from_taxonium:
    """Get metadata stored in viridian tree and add reference sequence information"""
    input:
        taxonium = rules.params.input.viridian_taxonium
    output:
        metadata = "results/viridian/metadata.tsv"
    params:
        name_reference = "MN908947",
        date_reference = "2019-12-26", # hospital admission date
        country_reference = "China",
        pango_reference = "B"
    conda:
        "../envs/bte.yaml"
    shell:
        """
        python3 scripts/mat-manipulation/get_metadata_from_taxonium.py \
            --input_taxonium {input.taxonium} \
            --output_metadata {output.metadata} \
            --name_reference {params.name_reference} \
            --date_reference {params.date_reference} \
            --country_reference {params.country_reference} \
            --pango_reference {params.pango_reference}
        """


rule extract_sra_ids:
    """
    Extract SRA IDs from viridian metadata
    """
    input:
        viridian_metadata = rules.get_metadata_from_taxonium.output.metadata
    output:
        sra_ids = "results/sra_metadata/sra_ids.txt"
    conda:
        "../envs/sra-metadata.yaml"
    shell:
        """
        # Extract SRA IDs from viridian metadata 'strain' column using tsv-select
        tsv-select -H -f strain {input.viridian_metadata} | \
            tail -n +2 > {output.sra_ids}
        """

rule filter_viridian_sra_metadata:
    """
    Filter Victor's SRA metadata to only include SRA IDs present in viridian tree
    """
    input:
        ncbi_metadata = "results/sra_metadata/sra_metadata_from_victor.csv",
        sra_ids = rules.extract_sra_ids.output.sra_ids,
    output:
        viridian_sra_metadata = "results/sra_metadata/sra_metadata_from_viridian.tsv",
    conda:
        "../envs/sra-metadata.yaml"
    shell:
        """
        csv2tsv {input.ncbi_metadata} |
            tsv-join -H --filter-file {input.sra_ids} -k 1 \
             > {output.viridian_sra_metadata}
        """
         
rule extract_biosample_ids:
    """
    Extract BioSample IDs from NCBI SRA_Accessions.tab.

    Streams the large NCBI SRA_Accessions.tab file (29 GB) and extracts only
    the Accession and BioSample columns for requested SRA IDs.

    Educational notes:
        - Set-based filtering: O(1) lookup for millions of IDs
    """
    input:
        sra_ids = rules.extract_sra_ids.output.sra_ids
    output:
        sra_biosample_map = "results/sra_metadata/sra_biosample_map.txt"
    params:
        ncbi_url = rules.params.params.ncbi_sra_url
    conda:
        "../envs/sra-metadata.yaml"
    shell:
        """
        # Extract BioSample IDs from NCBI SRA_Accessions.tab
        python3 scripts/sra-metadata/extract_biosample_ids.py \
            --sra-ids {input.sra_ids} \
            --output {output.sra_biosample_map} \
            --ncbi-url {params.ncbi_url} \
            --columns Accession BioSample
        """

rule get_sra_ids_not_found_metadata:
    """
    Identify SRA IDs from viridian metadata that have no corresponding metadata in NCBI tables
    """
    input:
        sra_ids = rules.extract_sra_ids.output.sra_ids,
        viridian_sra_metadata = rules.filter_viridian_sra_metadata.output.viridian_sra_metadata
    output:
        sra_ids_no_metadata = "results/sra_metadata/sra_ids_no_metadata.txt"
    conda:
        "../envs/sra-metadata.yaml"
    shell:
        """
        tsv-join --filter-file <(tsv-select -H -f 1 {input.viridian_sra_metadata} | tail -n +2) \
            --exclude {input.sra_ids} \
            > {output.sra_ids_no_metadata}
        """

rule get_biosample_sra_missing_metadata:
    """
    Map SRA IDs with no metadata back to BioSample IDs (for second download)
    """
    input:
        sra_ids_no_metadata = rules.get_sra_ids_not_found_metadata.output.sra_ids_no_metadata,
        sra_biosample_map = rules.extract_biosample_ids.output.sra_biosample_map
    output:
        biosamples_no_metadata = "results/sra_metadata/biosamples_no_metadata.txt"
    conda:
        "../envs/sra-metadata.yaml"
    shell:
        """
        tsv-join --filter-file {input.sra_ids_no_metadata} \
            -k 1 {input.sra_biosample_map} | \
            tsv-select -f 2 | tail -n +2 > {output.biosamples_no_metadata}
        """

rule query_missing_biosamples_entrez:
    """
    Query BioSample database in batches for missing biosample metadata.

    Uses -stop parameter to fetch all records (up to 10000).
    """
    input:
        biosamples_no_metadata = rules.get_biosample_sra_missing_metadata.output.biosamples_no_metadata
    output:
        missing_biosamples_xml = "results/sra_metadata/biosamples_missing_metadata.xml",
    conda:
        "../envs/sra-metadata.yaml"
    shell:
        """
        # Check if there are any biosamples to query
        if [ -s {input.biosamples_no_metadata} ] && grep -q -E '^SAM[NDEA]' {input.biosamples_no_metadata}; then
            # Count how many biosamples we have
            count=$(grep -c -E '^SAM[NDEA]' {input.biosamples_no_metadata})
            echo "Fetching $count biosamples from NCBI..."

            cat {input.biosamples_no_metadata} | grep -E '^SAM[NDEA]' \
                |  epost -db biosample \
                | efetch -format xml -start 1 -stop $count \
                > {output.missing_biosamples_xml}
        else
            # Create empty XML output if no biosamples to query
            echo '<?xml version="1.0" encoding="UTF-8"?>' > {output.missing_biosamples_xml}
            echo '<BioSampleSet></BioSampleSet>' >> {output.missing_biosamples_xml}
        fi
        """

rule parse_biosample_xml:
    """
    Parse BioSample XML to extract metadata fields
    """
    input:
        missing_biosamples_xml = rules.query_missing_biosamples_entrez.output.missing_biosamples_xml
    output:
        biosample_metadata = "results/sra_metadata/biosample_metadata.tsv"
    conda:
        "../envs/sra-metadata.yaml"
    shell:
        """
        echo -e "biosample_id\tcollection_date\tgeo_loc_name\tgeo_loc_region\tisolation_source\tpurpose_of_sequencing" \
            > {output.biosample_metadata}

        python3 scripts/sra-metadata/parse_concatenated_biosample.py {input.missing_biosamples_xml} \
            >> {output.biosample_metadata}
        """

rule merge_sra_metadata:
    """
    Merge metadata from Victor's CSV with newly queried BioSample metadata.

    Strategy:
    1. Start with Victor's comprehensive metadata (4.4M+ entries)
    2. For SRA IDs missing from Victor's data, add BioSample metadata via sra_biosample_map
    3. Create unified output with standardized columns

    This produces complete metadata coverage for all sequences in the viridian tree.
    """
    input:
        viridian_sra_metadata = rules.filter_viridian_sra_metadata.output.viridian_sra_metadata,
        biosample_metadata = rules.parse_biosample_xml.output.biosample_metadata,
        sra_biosample_map = "results/sra_metadata/sra_biosample_map.txt",
        sra_ids_no_metadata = rules.get_sra_ids_not_found_metadata.output.sra_ids_no_metadata
    output:
        merged_metadata = "results/sra_metadata/sra_metadata_merged.tsv"
    conda:
        "../envs/sra-metadata.yaml"
    shell:
        """
        # Step 1: Start with Victor's metadata and add isolation_source, purpose_of_sequencing, host_age, host_sex columns
        # Victor's format: acc, biosample, releasedate, collection_date, country, geo_loc_name, jattr (7 cols)
        # New format: add isolation_source, purpose_of_sequencing, host_age, host_sex before jattr (11 cols)
        awk -F'\t' 'BEGIN {{OFS="\t"}}
        NR==1 {{
            # Update header to include new columns
            print $1, $2, $3, $4, $5, $6, "isolation_source", "purpose_of_sequencing", "host_age", "host_sex", $7;
            next;
        }}
        {{
            # Add empty columns for isolation_source, purpose_of_sequencing, host_age, host_sex
            print $1, $2, $3, $4, $5, $6, "", "", "", "", $7;
        }}' {input.viridian_sra_metadata} > {output.merged_metadata}

        # Step 2: For missing SRA IDs, join with BioSample metadata
        #   - Join sra_ids_no_metadata with sra_biosample_map to get BioSample IDs
        #   - Join result with biosample_metadata to get metadata fields
        #   - Output 11 columns to match new format

        if [ -s {input.sra_ids_no_metadata} ]; then
            # Join missing SRA IDs -> BioSample IDs -> metadata
            # Step 2a: Get SRA-BioSample pairs for missing SRA IDs
            # Step 2b: Join with biosample metadata on BioSample ID
            #   Note: Use awk instead of tsv-join to avoid duplicate key issues
            tail -n +2 {input.sra_biosample_map} | \
                tsv-join -k 1 --filter-file {input.sra_ids_no_metadata} | \
                awk -F'\t' 'NR==FNR {{
                    # Load biosample metadata into array (key=biosample_id)
                    biosample_map[$1] = $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6;
                    next;
                }}
                {{
                    # Join SRA-BioSample pairs with metadata
                    sra=$1;
                    biosample=$2;
                    if (biosample in biosample_map) {{
                        split(biosample_map[biosample], fields, "\t");
                        collection_date = fields[1];
                        geo_loc = fields[2];
                        geo_region = fields[3];
                        isolation = fields[4];
                        purpose = fields[5];

                        # Extract country from geo_loc_name
                        split(geo_loc, parts, ":");
                        country = (parts[1] != "N/A") ? parts[1] : "";

                        # Create minimal jattr JSON
                        jattr = sprintf("{{\\"collection_date_sam\\": [\\"%s\\"], \\"geo_loc_name_sam\\": [\\"%s\\"], \\"isolation_source_sam\\": [\\"%s\\"], \\"purpose_of_sequencing_sam\\": [\\"%s\\"]}}",
                            collection_date, geo_loc, isolation, purpose);

                        # Output 11 columns: acc, biosample, releasedate, collection_date, country, geo_loc, isolation_source, purpose_of_sequencing, host_age, host_sex, jattr
                        print sra "\t" biosample "\t" "N/A" "\t" collection_date "\t" country "\t" geo_loc "\t" isolation "\t" purpose "\t" "" "\t" "" "\t" jattr;
                    }}
                }}' <(tail -n +2 {input.biosample_metadata}) - \
                >> {output.merged_metadata}
        fi

        echo "Merged metadata created with $(tail -n +2 {output.merged_metadata} | wc -l) total entries"
        """

rule extract_dates_from_json:
    """
    Extract collection dates from jattr JSON when the collection_date_sam column is empty.

    Many entries in Victor's metadata have collection_date_sam empty in column 4,
    but the date is available in the jattr JSON field (column 7).

    This rule fills in missing dates by parsing the JSON.
    """
    input:
        merged_metadata = rules.merge_sra_metadata.output.merged_metadata
    output:
        metadata_with_dates = "results/sra_metadata/sra_metadata_complete.tsv"
    conda:
        "../envs/sra-metadata.yaml"
    shell:
        """
        # Use awk to extract dates from JSON when column 4 is empty
        # Now handles 11 columns: acc, biosample, releasedate, collection_date, country, geo_loc, isolation_source, purpose_of_sequencing, host_age, host_sex, jattr
        awk -F'\t' 'BEGIN {{OFS="\t"}}
        NR==1 {{
            # Print header as-is
            print $0;
            next;
        }}
        {{
            acc = $1;
            biosample = $2;
            releasedate = $3;
            collection_date = $4;
            country = $5;
            geo_loc = $6;
            isolation_source = $7;
            purpose_of_sequencing = $8;
            host_age = $9;
            host_sex = $10;
            jattr = $11;

            # If collection_date is empty or whitespace-only, try to extract from JSON
            if (collection_date == "" || collection_date ~ /^[[:space:]]*$/) {{
                # Extract collection_date_sam from JSON: "collection_date_sam": ["2021"]
                # Use simpler regex that works with all awk versions
                if (match(jattr, /"collection_date_sam"[[:space:]]*:[[:space:]]*\["[^"]+"]/)) {{
                    date_str = substr(jattr, RSTART, RLENGTH);
                    if (match(date_str, /\["[^"]+"\]/)) {{
                        date_val = substr(date_str, RSTART+2, RLENGTH-4);
                        collection_date = date_val;
                    }}
                }}
                # Also try without array brackets: "collection_date_sam": "2021"
                else if (match(jattr, /"collection_date_sam"[[:space:]]*:[[:space:]]*"[^"]+"/)) {{
                    date_str = substr(jattr, RSTART, RLENGTH);
                    if (match(date_str, /"[^"]+"/)) {{
                        date_val = substr(date_str, RSTART+1, RLENGTH-2);
                        if (date_val != "collection_date_sam") collection_date = date_val;
                    }}
                }}
            }}

            # If isolation_source is empty, try to extract from JSON
            if (isolation_source == "" || isolation_source ~ /^[[:space:]]*$/) {{
                if (match(jattr, /"isolation_source_sam"[[:space:]]*:[[:space:]]*\["[^"]+"]/)) {{
                    iso_str = substr(jattr, RSTART, RLENGTH);
                    if (match(iso_str, /\["[^"]+"\]/)) {{
                        iso_val = substr(iso_str, RSTART+2, RLENGTH-4);
                        isolation_source = iso_val;
                    }}
                }}
            }}

            # If purpose_of_sequencing is empty, try to extract from JSON
            if (purpose_of_sequencing == "" || purpose_of_sequencing ~ /^[[:space:]]*$/) {{
                # Try array format first: "purpose_of_sequencing_sam": ["value"]
                if (match(jattr, /"purpose_of_sequencing_sam"[[:space:]]*:[[:space:]]*\["[^"]+"]/)) {{
                    pur_str = substr(jattr, RSTART, RLENGTH);
                    if (match(pur_str, /\["[^"]+"\]/)) {{
                        pur_val = substr(pur_str, RSTART+2, RLENGTH-4);
                        purpose_of_sequencing = pur_val;
                    }}
                }}
                # Try plain string format: "purpose_of_sequencing_sam": "value"
                else if (match(jattr, /"purpose_of_sequencing_sam"[[:space:]]*:[[:space:]]*"[^"]+"/)) {{
                    pur_str = substr(jattr, RSTART, RLENGTH);
                    # Match from colon onwards to skip the field name
                    if (match(pur_str, /:[[:space:]]*"[^"]+"/)) {{
                        val_str = substr(pur_str, RSTART, RLENGTH);
                        if (match(val_str, /"[^"]+"/)) {{
                            pur_val = substr(val_str, RSTART+1, RLENGTH-2);
                            purpose_of_sequencing = pur_val;
                        }}
                    }}
                }}
            }}

            # If host_age is empty, try to extract from JSON
            # host_age_sam uses array format: ["56"], ["64"], etc.
            if (host_age == "" || host_age ~ /^[[:space:]]*$/) {{
                if (match(jattr, /"host_age_sam"[[:space:]]*:[[:space:]]*\["[^"]+"]/)) {{
                    age_str = substr(jattr, RSTART, RLENGTH);
                    if (match(age_str, /\["[^"]+"\]/)) {{
                        age_val = substr(age_str, RSTART+2, RLENGTH-4);
                        host_age = age_val;
                    }}
                }}
            }}

            # If host_sex is empty, try to extract from JSON
            # host_sex_sam uses array format: ["male"], ["female"], etc.
            if (host_sex == "" || host_sex ~ /^[[:space:]]*$/) {{
                if (match(jattr, /"host_sex_sam"[[:space:]]*:[[:space:]]*\["[^"]+"]/)) {{
                    sex_str = substr(jattr, RSTART, RLENGTH);
                    if (match(sex_str, /\["[^"]+"\]/)) {{
                        sex_val = substr(sex_str, RSTART+2, RLENGTH-4);
                        host_sex = sex_val;
                    }}
                }}
            }}

            print acc, biosample, releasedate, collection_date, country, geo_loc, isolation_source, purpose_of_sequencing, host_age, host_sex, jattr;
        }}' {input.merged_metadata} > {output.metadata_with_dates}

        # Report how many fields were filled in from JSON
        echo "=== Metadata extraction from JSON ==="

        # Collection dates
        empty_dates_before=$(awk -F'\t' 'NR>1 && ($4=="" || $4~/^[[:space:]]*$/)' {input.merged_metadata} | wc -l)
        empty_dates_after=$(awk -F'\t' 'NR>1 && ($4=="" || $4~/^[[:space:]]*$/)' {output.metadata_with_dates} | wc -l)
        filled_dates=$((empty_dates_before - empty_dates_after))
        echo "Collection dates: Filled $filled_dates (was $empty_dates_before empty, now $empty_dates_after empty)"

        # Isolation source
        empty_iso_before=$(awk -F'\t' 'NR>1 && ($7=="" || $7~/^[[:space:]]*$/)' {input.merged_metadata} | wc -l)
        empty_iso_after=$(awk -F'\t' 'NR>1 && ($7=="" || $7~/^[[:space:]]*$/)' {output.metadata_with_dates} | wc -l)
        filled_iso=$((empty_iso_before - empty_iso_after))
        echo "Isolation source: Filled $filled_iso (was $empty_iso_before empty, now $empty_iso_after empty)"

        # Purpose of sequencing
        empty_pur_before=$(awk -F'\t' 'NR>1 && ($8=="" || $8~/^[[:space:]]*$/)' {input.merged_metadata} | wc -l)
        empty_pur_after=$(awk -F'\t' 'NR>1 && ($8=="" || $8~/^[[:space:]]*$/)' {output.metadata_with_dates} | wc -l)
        filled_pur=$((empty_pur_before - empty_pur_after))
        echo "Purpose of sequencing: Filled $filled_pur (was $empty_pur_before empty, now $empty_pur_after empty)"

        # Host age
        empty_age_before=$(awk -F'\t' 'NR>1 && ($9=="" || $9~/^[[:space:]]*$/)' {input.merged_metadata} | wc -l)
        empty_age_after=$(awk -F'\t' 'NR>1 && ($9=="" || $9~/^[[:space:]]*$/)' {output.metadata_with_dates} | wc -l)
        filled_age=$((empty_age_before - empty_age_after))
        echo "Host age: Filled $filled_age (was $empty_age_before empty, now $empty_age_after empty)"

        # Host sex
        empty_sex_before=$(awk -F'\t' 'NR>1 && ($10=="" || $10~/^[[:space:]]*$/)' {input.merged_metadata} | wc -l)
        empty_sex_after=$(awk -F'\t' 'NR>1 && ($10=="" || $10~/^[[:space:]]*$/)' {output.metadata_with_dates} | wc -l)
        filled_sex=$((empty_sex_before - empty_sex_after))
        echo "Host sex: Filled $filled_sex (was $empty_sex_before empty, now $empty_sex_after empty)"
        """

rule identify_sra_requiring_manual_input:
    """
    Identify SRA IDs still missing after Victor's metadata + NCBI query.

    These are SRA IDs that:
    - Are in the Viridian tree
    - NOT in sra_metadata_complete (i.e., not found in Victor's data or NCBI)

    These will require manual metadata curation.
    """
    input:
        viridian_metadata = rules.get_metadata_from_taxonium.output.metadata,
        metadata_complete = rules.extract_dates_from_json.output.metadata_with_dates
    output:
        sra_ids_requiring_manual_input = "results/sra_metadata/sra_ids_requiring_manual_input.txt"
    conda:
        "../envs/sra-metadata.yaml"
    shell:
        """
        # Find SRA IDs in Viridian tree but not in metadata_complete
        tsv-join -H --filter-file {input.metadata_complete} --exclude -k 1 {input.viridian_metadata} | \
            tail -n +2 | tsv-select -f 1 > {output.sra_ids_requiring_manual_input}

        count=$(wc -l < {output.sra_ids_requiring_manual_input})
        echo "=== Identified $count SRA IDs requiring manual metadata input ==="
        """

rule tidy_final_metadata:
    """
    Tidy and filter metadata to prepare final dataset.
    """
    input:
        metadata = rules.extract_dates_from_json.output.metadata_with_dates,
        manual_input_metadata = "data/missing_metadata.txt",
        strains_viridian = rules.extract_sra_ids.output.sra_ids,
    output:
        strains_to_keep = "results/sra_metadata/sra_strains_to_keep.txt",
        subset_metadata = "results/sra_metadata/sra_metadata_final.tsv"
    conda:
        "../envs/sra-metadata.yaml"
    shell:
        """
        Rscript scripts/sra-metadata/tidy_final_metadata.R \
            --metadata {input.metadata} \
            --missing_metadata {input.manual_input_metadata} \
            --strains_viridian {input.strains_viridian} \
            --output_strains_to_keep {output.strains_to_keep} \
            --output_subset_metadata {output.subset_metadata}
        """


rule clean_dates:
    """
    Remove entries with dates that are too old
    """
    input:
        metadata = rules.tidy_final_metadata.output.subset_metadata,
    output:
        metadata = "results/sra_metadata/sra_metadata_just_clean_dates.tsv"
    conda:
        "../envs/sra-metadata.yaml"
    shell:
        """
        Rscript scripts/sra-metadata/remove_incorrect_dates.R \
            --metadata {input.metadata} \
            --output_metadata {output.metadata}
        """


