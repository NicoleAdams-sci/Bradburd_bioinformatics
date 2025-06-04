#!/bin/bash

check_invariant_sites() {
    local vcf_file="$1"
    
    # Check if the file exists
    if [ ! -f "$vcf_file" ]
    then
        echo "Error: $vcf_file does not exist."
        return 1
    fi
    
    # If file is gzipped, use zcat, otherwise use cat
    if [[ "$vcf_file" == *.gz ]]
    then
        READ_COMMAND="zcat"
    else
        READ_COMMAND="cat"
    fi
    
    # Count invariant sites - look for sites where REF exists but ALT is '.'
    invariant_count=$($READ_COMMAND "$vcf_file" | \
        grep -v '^#' | \
        awk '{if ($5 == ".") print}' | \
        wc -l)
    
    # Count total variant sites (excluding header)
    total_sites=$($READ_COMMAND "$vcf_file" | grep -v '^#' | wc -l)
    
    if [ "$invariant_count" -gt 0 ]
    then
        echo "Found $invariant_count invariant sites out of $total_sites total sites"
        return 0
    else
        echo "No invariant sites found in $total_sites total sites"
        return 1
    fi
}

# Usage
if [ "$#" -ne 1 ]
then
    echo "Usage: $0 <vcf_file>"
    exit 1
fi

check_invariant_sites "$1"
exit $?
