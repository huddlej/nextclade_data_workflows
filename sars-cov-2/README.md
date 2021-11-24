# Using Nextclade as a pango classifier

## Plan

1. Move trait-reconstruction script and integrate into snakemake workflow
1. Download newest alias key
1. Download newest designations
1. Run pangolin for all tips
1. Reconstruct pango states for internal nodes
1. Turn pango labels into node json
1. Use pango node json instead of clade node json in export
