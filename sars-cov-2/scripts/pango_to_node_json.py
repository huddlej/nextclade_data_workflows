import click
import json
import pandas as pd

@click.command()
@click.argument('input_file', type=click.File("r"))
@click.argument('output_file', type=click.File("w"))
def convert(input_file, output_file):
    pango_df = pd.read_csv(input_file)[["taxon", "lineage"]]

    nodes = {"nodes": {}}
    nodes_nodes = nodes["nodes"]
    for row in pango_df.itertuples():
        nodes_nodes[row[1]] = {"clade_membership": row[2]}
    json.dump(nodes, output_file)

if __name__ == "__main__":
    convert()