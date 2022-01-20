import csv
import json
from collections import OrderedDict

import click
import funcy as fy
import pandas as pd
from Bio import Phylo


def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names


def get_terminal_names(clade):
    terminal_list = []
    for tip in clade.get_terminals():
        terminal_list.append(tip.name)
    return terminal_list


def get_lineage_counts(terminal_names, meta):
    lineages = meta.loc[meta.index.isin(terminal_names), "lineage"].value_counts().to_dict()
    return lineages


def dealias_lineage_name(name, alias_dict):
    name_split = name.split(".")
    letter = name_split[0]
    dealiased_letter = alias_dict[letter]
    if type(dealiased_letter) == list:
        return letter
    if len(name_split) > 2:
        return dealiased_letter + "." + ".".join(name_split[1:])
    if len(name_split) == 2:
        return dealiased_letter + "." + name_split[1]
    else:
        return dealiased_letter


def get_candidate_lineage_counts(value_counts):
    candidates = OrderedDict()
    for lineage in value_counts:
        lineage_split = lineage.split(".")
        for i in range(len(lineage_split) + 1):
            candidates[".".join(lineage_split[:i])] = (
                candidates.get(".".join(lineage_split[:i]), 0) + value_counts[lineage]
            )
    return candidates


def realias(dealiased, alias_dict):
    if dealiased == "":
        return "A/B"
    name_split = dealiased.split(".")
    if len(name_split) <= 4:
        return dealiased
    to_lookup = ".".join(name_split[:4])
    alias = [k for k, v in alias_dict.items() if v == to_lookup][0]
    if len(name_split) > 5:
        return alias + "." + ".".join(name_split[4:])
    else:
        return alias + "." + name_split[4]


def get_consensus_lineage(candidates):
    # TODO: #2 figure out how to deal with recombinants like XA
    # Do not dealias Xs
    # Remove Xs if it's not all descendent from an X
    count = candidates[""]
    # recombinant_count = recombinant_count(candidates)
    consensus_lineages = fy.select_values(lambda x: x == count, candidates)
    consensus_lineage = next(reversed(consensus_lineages))
    return consensus_lineage


@click.command()
@click.option("-t", "--tree", type=click.File("r"), default="builds/pango/tree.nwk")
@click.option("-d", "--designations", type=click.File("r"), default="pre-processed/pango_designations_nextstrain_names.csv")
@click.option(
    "-o", "--outfile", type=click.File("w"), default="builds/nextclade/pango_reconstructed.json"
)
@click.option("-a", "--aliases", type=click.File("r"), default="pre-processed/alias.json")
def main(tree, designations, outfile, aliases):
    alias_dict = json.load(aliases)
    alias_dict["A"] = "A"
    alias_dict["B"] = "B"
    tree = Phylo.read(tree, "newick")
    meta = pd.read_csv(designations)[["strain", "lineage"]]
    meta.set_index("strain", inplace=True)
    output = {"nodes": {}}
    for tip in tree.get_terminals():
        try:
            output["nodes"][tip.name] = {"clade_membership": meta.loc[tip.name]["lineage"]}
        except:
            print(f"{tip.name} not found in designations")

    meta["lineage"] = meta["lineage"].apply(lambda x: dealias_lineage_name(x, alias_dict))
    
    lookup_dict = lookup_by_names(tree)
    internals = []
    for internal in tree.get_nonterminals():
        internals.append(internal.name)
    internal_lineages = []

    for internal in internals:
        try:
            terminal_names = get_terminal_names(lookup_dict[internal])
            lineage_counts = get_lineage_counts(terminal_names, meta)
            candidate_lineage_counts = get_candidate_lineage_counts(lineage_counts)
            consensus_lineage = get_consensus_lineage(candidate_lineage_counts)
            lineage = realias(consensus_lineage, alias_dict)
            internal_lineages.append({"node": internal, "lineage": lineage})
            output["nodes"][internal] = {"clade_membership": lineage}
        except:
            print(terminal_names, lineage_counts, candidate_lineage_counts, consensus_lineage)


    json.dump(output, outfile, indent=4)


if __name__ == "__main__":
    main()
