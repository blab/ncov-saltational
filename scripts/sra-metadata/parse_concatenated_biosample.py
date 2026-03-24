#!/usr/bin/env python3
"""
Parse concatenated BioSample XML documents.

Handles cases where efetch returns multiple XML documents concatenated together.
Splits by </BioSampleSet> boundaries and parses each section separately.
"""

import sys
import re
import xml.etree.ElementTree as ET

def split_concatenated_xml(xml_file):
    """
    Split concatenated XML documents at BioSampleSet boundaries.

    Returns list of XML strings, each containing one complete BioSampleSet document.
    """
    with open(xml_file, 'r') as f:
        content = f.read()

    # Split at </BioSampleSet> and reconstruct complete documents
    parts = content.split('</BioSampleSet>')

    documents = []
    for i, part in enumerate(parts[:-1]):  # Skip last empty part
        # Find where this document starts (at <BioSampleSet>)
        biosampleset_start = part.rfind('<BioSampleSet>')

        if biosampleset_start == -1:
            continue

        # Extract everything from <BioSampleSet> onwards
        doc = part[biosampleset_start:] + '</BioSampleSet>'

        # Add XML declaration if not present
        if not doc.strip().startswith('<?xml'):
            doc = '<?xml version="1.0" encoding="UTF-8" ?>\n<!DOCTYPE BioSampleSet>\n' + doc

        documents.append(doc)

    return documents

def parse_biosample_xml(xml_string):
    """Parse a single BioSampleSet XML document and extract metadata."""
    results = []

    try:
        root = ET.fromstring(xml_string)

        # Iterate through all BioSample elements
        for biosample in root.findall('.//BioSample'):
            # Get BioSample accession
            biosample_id = biosample.get('accession')

            if not biosample_id:
                id_element = biosample.find('.//Id[@db="BioSample"]')
                if id_element is not None and id_element.text:
                    biosample_id = id_element.text

            # Initialize metadata dictionary
            metadata = {
                'biosample_id': biosample_id or 'N/A',
                'collection_date': 'N/A',
                'geo_loc_name': 'N/A',
                'geo_loc_region': 'N/A',
                'isolation_source': 'N/A',
                'purpose_of_sequencing': 'N/A'
            }

            # Extract attributes with support for multiple naming conventions
            for attribute in biosample.findall('.//Attribute'):
                attr_name = attribute.get('attribute_name')
                attr_value = attribute.text

                if not attr_value:
                    continue

                # Collection date: support both underscore and space variants
                if attr_name in ('collection_date', 'collection date'):
                    if metadata['collection_date'] == 'N/A':
                        metadata['collection_date'] = attr_value

                # Geographic location (country): multiple naming conventions
                elif attr_name in ('geo_loc_name', 'geographic location (country and/or sea)'):
                    if metadata['geo_loc_name'] == 'N/A':
                        metadata['geo_loc_name'] = attr_value

                # Geographic location (region/locality): more detailed location info
                elif attr_name == 'geographic location (region and locality)':
                    if metadata['geo_loc_region'] == 'N/A':
                        metadata['geo_loc_region'] = attr_value

                # Isolation source: multiple naming conventions
                elif attr_name in ('isolation_source', 'isolation source',
                                  'isolation source host associated',
                                  'isolation source host-associated',
                                  'isolation source-host associated'):
                    if metadata['isolation_source'] == 'N/A':
                        metadata['isolation_source'] = attr_value

                # Purpose of sequencing
                elif attr_name == 'purpose_of_sequencing':
                    if metadata['purpose_of_sequencing'] == 'N/A':
                        metadata['purpose_of_sequencing'] = attr_value

            results.append(metadata)

    except Exception as e:
        print(f"Error parsing XML document: {e}", file=sys.stderr)

    return results

def write_tsv_rows(records):
    """Write BioSample records as TSV rows to stdout."""
    fields = ['biosample_id', 'collection_date', 'geo_loc_name', 'geo_loc_region',
              'isolation_source', 'purpose_of_sequencing']

    for record in records:
        row = '\t'.join(str(record[field]) for field in fields)
        print(row)

def main():
    if len(sys.argv) != 2:
        print("Usage: python parse_concatenated_biosample.py <biosample.xml>", file=sys.stderr)
        sys.exit(1)

    xml_file = sys.argv[1]

    # Split concatenated XML documents
    print(f"Splitting concatenated XML documents...", file=sys.stderr)
    documents = split_concatenated_xml(xml_file)
    print(f"Found {len(documents)} XML documents", file=sys.stderr)

    # Parse each document
    all_records = []
    for i, doc in enumerate(documents, 1):
        print(f"Parsing document {i}/{len(documents)}...", file=sys.stderr)
        records = parse_biosample_xml(doc)
        all_records.extend(records)
        print(f"  Extracted {len(records)} BioSample records", file=sys.stderr)

    # Write TSV rows to stdout
    write_tsv_rows(all_records)

    print(f"\nTotal: Processed {len(all_records)} BioSample records", file=sys.stderr)

if __name__ == "__main__":
    main()
