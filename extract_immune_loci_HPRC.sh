#!/usr/bin/env bash
set -euo pipefail

# ============
# Input args
# ============
input="$1"      # e.g. merged.fasta.gz i.e. concat mat/pat.fa.gz or hap1/hap2.fa.gz 
sample="$2"     # e.g. sample name
threads=36

# ============
# Step 1. Decompress input if needed
# ============
if [[ "$input" == *.gz ]]; then
    echo "Decompressing $input ..."
    gunzip -f "$input"
    input="${input%.gz}"
fi

# ============
# Step 2. Build BLAST database
# ============
echo "Building BLAST database from $input ..."
makeblastdb -in "$input" -dbtype nucl -out subject_db

# ============
# Step 3. Run BLAST
# ============
echo "Running blastn ..."
blastn \
    -query flank_sequences.fasta \ #immune loci specific flanking sequences
    -db subject_db \
    -outfmt "6 length pident nident mismatch gapopen gaps qseqid qstart qend qlen sseqid sstart send slen sstrand" \
    -num_threads "$threads" \
    > all_ref_blast_result.txt

# ============
# Step 4. Filter BLAST results
# ============
echo "Filtering BLAST hits ..."
cat all_ref_blast_result.txt \
    | sort -n -r -k 1 \
    | awk '{ if ($1 > 1000 && $2 >= 97) print $0 }' \
    > all_ref_blast_result_filtered.txt

# ============
# Step 5. Python inline script to extract coordinates
# ============
echo "Extracting coordinates ..."
python3 - <<'PYCODE' all_ref_blast_result_filtered.txt > coords.txt
import sys

input_file = sys.argv[1]
temp_dict = {}

with open(input_file) as f:
    for line in f:
        line = line.strip()
        parts = line.split('\t')
        loc = parts[6]
        sseqid = parts[10]
        sstart, send, slen, strand = parts[11], parts[12], parts[13], parts[14]

        if 'igk' not in loc:
            key = f"{sseqid}|{loc}"
        else:
            key = sseqid

        if key not in temp_dict:
            temp_dict[key] = [[sstart], [send], slen]
        else:
            if strand == 'minus':
                temp_dict[key][0].append(send)
                temp_dict[key][1].append(sstart)
            else:
                temp_dict[key][0].append(sstart)
                temp_dict[key][1].append(send)

coord_dict = {
    k: [min(map(int, v[0])), max(map(int, v[1]))]
    for k, v in temp_dict.items()
}

for key, (start, end) in coord_dict.items():
    print(f"{key.split('|')[0]}:{start}-{end}")
PYCODE

# ============
# Step 6. Extract sequences with samtools faidx
# ============
echo "Extracting sequences ..."
samtools faidx "$input"
while read coords; do
    samtools faidx "$input" "$coords"
done < coords.txt > "${sample}_all_immune_loci.fa"

# ============
# Step 7. Cleanup intermediate files
# ============
echo "Cleaning up intermediate files ..."
rm -f subject_db.* all_ref_blast_result.txt all_ref_blast_result_filtered.txt coords.txt
echo "Done! Output FASTA: ${sample}_all_immune_loci.fa"
