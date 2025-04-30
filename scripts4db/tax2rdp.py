#!/usr/bin/env python3

import sys
import os
import re
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Convert RDP output style tax file to RDP training input and UTAX output styles, "
                    "and process a corresponding FASTA file.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("taxonomy_file", help="Input taxonomy file (RDP output style: gid<TAB>lineage)")
    parser.add_argument("fasta_file", help="Input sequences file (FASTA format)")
    parser.add_argument(
        "prefix",
        nargs='?',
        default="taxonomy",
        help="Prefix for output files (default: 'taxonomy')"
    )
    args = parser.parse_args()

    prefix = args.prefix
    ranks = ["rootrank", "kingdom", "phylum", "class", "order", "family", "genus", "species",
             "subspecies", "strain", "sub1", "sub2", "sub3", "sub4", "sub5", "sub6"]

    full_lineage = {}
    entries = {}
    id_counter = 0
    gid_to_taxid_map = {}

    print(f"Processing taxonomy file: {args.taxonomy_file}")
    try:
        with open(args.taxonomy_file, 'r', encoding='utf-8') as infile:
            for line_num, line in enumerate(infile):
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t', 1)
                if len(parts) != 2:
                    print(f"Warning: Skipping malformed line {line_num + 1} in taxonomy file: {line}", file=sys.stderr)
                    continue
                gid_str, lineage_str = parts

                gid = str(gid_str)
                full_lineage[gid] = lineage_str

                levels = lineage_str.split(';')
                parent_id = -1
                current_lineage_key = ""

                for i, level_name in enumerate(levels):
                    if i == 0:
                         current_lineage_key = level_name
                    else:
                         current_lineage_key += ";" + level_name

                    if current_lineage_key not in entries:
                        entries[current_lineage_key] = {
                            'id': id_counter,
                            'name': level_name,
                            'parent': parent_id,
                            'depth': i
                        }
                        id_counter += 1

                    parent_id = entries[current_lineage_key]['id']

                if levels:
                    last_level_key = lineage_str
                    if last_level_key in entries:
                        gid_to_taxid_map[gid] = entries[last_level_key]['id']
                    else:
                        print(f"Warning: Could not find entry for last level of {gid}: {last_level_key}", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: Taxonomy file not found: {args.taxonomy_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing taxonomy file {args.taxonomy_file}: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Processed {len(full_lineage)} taxonomy entries, found {len(entries)} unique rank nodes.")

    rdp_tax_filename = f"{prefix}.rdp.tax"
    utax_tax_filename = f"{prefix}.utax.tax"

    print(f"Writing RDP taxonomy to: {rdp_tax_filename}")
    print(f"Writing UTAX taxonomy to: {utax_tax_filename}")

    try:
        with open(rdp_tax_filename, 'w', encoding='utf-8') as rdp_out, \
             open(utax_tax_filename, 'w', encoding='utf-8') as utax_out:

            sorted_keys = sorted(entries.keys(), key=lambda k: entries[k]['id'])

            for key in sorted_keys:
                entry = entries[key]
                depth = entry['depth']
                rank_name = ranks[depth] if depth < len(ranks) else f"rank_{depth}"

                rdp_line = f"{entry['id']}*{entry['name']}*{entry['parent']}*{entry['depth']}*{rank_name}\n"
                rdp_out.write(rdp_line)

                utax_line = f"{entry['id']}\t{entry['parent']}\t{entry['name']}\t{rank_name}\n"
                utax_out.write(utax_line)

    except IOError as e:
        print(f"Error writing taxonomy output files: {e}", file=sys.stderr)
        sys.exit(1)

    map_filename = f"{prefix}.gid_taxid.map"

    print(f"Writing GID-to-TaxID map to: {map_filename}")
    try:
        with open(map_filename, 'w', encoding='utf-8') as map_out:
            for gid in sorted(gid_to_taxid_map.keys()):
                taxid = gid_to_taxid_map[gid]
                map_out.write(f"{gid}\t{taxid}\n")
    except IOError as e:
        print(f"Error writing mapping file: {e}", file=sys.stderr)
        sys.exit(1)

    rdp_fa_filename = f"{prefix}.rdp.fa"
    utax_fa_filename = f"{prefix}.utax.fa"

    print(f"Processing FASTA file: {args.fasta_file}")
    print(f"Writing RDP FASTA to: {rdp_fa_filename}")
    print(f"Writing UTAX FASTA to: {utax_fa_filename}")

    try:
        with open(args.fasta_file, 'r', encoding='utf-8') as fa_in, \
             open(rdp_fa_filename, 'w', encoding='utf-8') as rdp_fa_out, \
             open(utax_fa_filename, 'w', encoding='utf-8') as utax_fa_out:

            print_ok = False
            current_gid = None

            for line in fa_in:
                if line.startswith('>'):
                    print_ok = False
                    current_gid = None
                    header_content = line[1:].strip()

                    match = re.match(r'^(\S+)', header_content)
                    if match:
                        gid = match.group(1)
                        if gid in full_lineage and gid in gid_to_taxid_map:
                            print_ok = True
                            current_gid = gid
                            lineage = full_lineage[gid]
                            taxid = gid_to_taxid_map[gid]

                            rdp_fa_out.write(f">{gid} {lineage}\n")
                            utax_fa_out.write(f">{gid};tax={taxid}; {lineage}\n")
                        else:
                            print(f"No tax information found for ID: {gid}, skipping sequence...", file=sys.stderr)
                    else:
                         print(f"Warning: Could not parse ID from FASTA header: {line.strip()}", file=sys.stderr)

                elif print_ok and current_gid is not None:
                    sequence_line = line.strip()
                    if sequence_line:
                        rdp_fa_out.write(sequence_line + "\n")
                        utax_fa_out.write(sequence_line + "\n")

    except FileNotFoundError:
        print(f"Error: FASTA file not found: {args.fasta_file}", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"Error reading/writing FASTA files: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
         print(f"An unexpected error occurred during FASTA processing: {e}", file=sys.stderr)
         sys.exit(1)

    print("Processing complete.")

if __name__ == "__main__":
    main()
