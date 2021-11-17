import click

@click.command()
@click.option('--input-file', default="pre-processed/aa_substitutions.tsv", type=click.File('r'), help='Input file')
@click.option('--output-file', default="pre-processed/stops.tsv", type=click.File('w'), help='Output file')
def filter_stops(input_file, output_file):
    output = "stops\n"
    for line in input_file:
        subs = line[:-1].split(",")
        stops = []
        for sub in subs:
            if sub[-1] == "*":
                output += sub + "\n"

    output_file.write(output)

if __name__ == '__main__':
    filter_stops()