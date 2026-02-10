#!/bin/bash

INROOT="data/domain_large"
OUTROOT="data/domain_large_regridded"

# Grid description file
GRIDFILE="data/grid"

# Create output root folder if it doesn't exist
if [ ! -d "$OUTROOT" ]; then
    echo "Output root folder '$OUTROOT' not found. Creating it..."
    mkdir -p "$OUTROOT"
fi

# Loop through all experiment subfolders
for EXPDIR in "$INROOT"/*/; do
    # Extract experiment name (strip trailing slash)
    EXPNAME=$(basename "$EXPDIR")

    echo "-------- Processing experiment: $EXPNAME --------"

    # Create output subfolder
    OUTDIR="$OUTROOT/$EXPNAME"
    mkdir -p "$OUTDIR"

    # Loop through all NetCDF files in this experiment folder
    for FILE in "$EXPDIR"/*.nc; do
        # Skip if no .nc files exist
        [ -e "$FILE" ] || continue

        FILENAME=$(basename "$FILE")
        OUTFILE="$OUTDIR/$FILENAME"

        echo "  Regridding $FILENAME..."

        # Run CDO bilinear remapping
        cdo -P 10 remapbil,"$GRIDFILE" "$FILE" "$OUTFILE"
    done
done

echo "-------- Regridding Complete --------"