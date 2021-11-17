import pandas as pd
import click

@click.command()
@click.option('--number', default=20, help='Number of frameshifts')
@click.option('--input-file', default="pre-processed/frameshifts.tsv", help='Input file')
def format(number, input_file):
    df = pd.read_csv(input_file, sep="\t")

    #%%
    for row in df[1:number].sort_values(by=['frame_shifts']).itertuples():
        fs = row.frame_shifts.split(",")[0].split(":")
        gene = fs[0]
        loc = fs[1].split("-")
        # beware: subtract 1 because in qc.json it's 0 indexed as opposed to nextclade output
        start = int(loc[0])-1
        end = int(loc[1])
        click.echo(f'{{ "geneName": "{gene}", "codonRange": {{"begin": {start}, "end": {end} }} }},')

if __name__ == '__main__':
    format()