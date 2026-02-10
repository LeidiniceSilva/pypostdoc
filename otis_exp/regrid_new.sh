#!/bin/bash


INROOT="data/domain_small"
OUTROOT="data/domain_small_regridded"


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


    # Loop through all NetCDF files starting with uas_ or vas_
    for FILE in "$EXPDIR"/uas_*.nc "$EXPDIR"/vas_*.nc; do
        # Skip if no matching .nc files exist
        [ -e "$FILE" ] || continue


        FILENAME=$(basename "$FILE")
        OUTFILE="$OUTDIR/$FILENAME"


        echo "  Regridding $FILENAME..."


        # Run CDO bilinear remapping
        cdo -P 10 remapbil,"$GRIDFILE" "$FILE" "$OUTFILE"
    done
done


echo "-------- Regridding Complete --------"
