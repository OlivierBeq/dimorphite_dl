# -*- coding: utf-8 -*-

"""This entry point allows the use of dimorphite-dl as a CLI."""

import argparse
from typing import Optional

import click

from .dimorphite_dl import *


examples = """Examples:

    python dimorphite_dl.py --smiles_file sample_molecules.smi
    python dimorphite_dl.py --smiles "CCC(=O)O" --min_ph -3.0 --max_ph -2.0
    python dimorphite_dl.py --smiles "CCCN" --min_ph -3.0 --max_ph -2.0 --output_file output.smi
    python dimorphite_dl.py --smiles_file sample_molecules.smi --pka_precision 2.0 --label_states
    python dimorphite_dl.py --test"""


class Mutex(click.Option):
    def __init__(self, *args, **kwargs):
        """Custom class allowing click.Options to be
        required if other click.Options are not set.

        Derived from: https://stackoverflow.com/a/61684480
        """
        self.not_required_if: list = kwargs.pop("not_required_if")

        assert self.not_required_if, "'not_required_if' parameter required"
        assert isinstance(self.not_required_if, list), "'not_required_if' mut be a list"
        kwargs["help"] = (kwargs.get("help", "") + ' NOTE: This argument is mutually exclusive with ' + ", ".join(
            self.not_required_if) + ".").strip()
        super(Mutex, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        """Override base method."""
        current_opt: bool = self.consume_value(ctx, opts)[0]  # Obtain set value
        for other_param in ctx.command.get_params(ctx):
            if other_param is self:  # Ignore the current parameter
                continue
            # Other argument's name or declaration in self.not_required_if
            if other_param.human_readable_name in self.not_required_if or any(
                    opt.lstrip('-') in self.not_required_if for opt in other_param.opts) or any(
                    opt.lstrip('-') in self.not_required_if for opt in other_param.secondary_opts):
                # Get value assigned to the other option
                other_opt: bool = other_param.consume_value(ctx, opts)[0]
                if other_opt:
                    if current_opt:
                        raise click.UsageError(
                            "Illegal usage: '" + str(self.name)
                            + "' is mutually exclusive with "
                            + str(other_param.human_readable_name) + "."
                        )
                    else:
                        self.required = None  # Override requirement
        return super(Mutex, self).handle_parse_result(ctx, opts, args)


class SpecialEpilog(click.Group):
    """Custom class with a custom display of the epilog."""
    def format_epilog(self, ctx, formatter):
        if self.epilog:
            formatter.width = float('inf')
            formatter.write_paragraph()
            for line in self.epilog.split('\n'):
                formatter.write_text(line)

    def get_help(self, ctx):
        """ standard get help, but without rstrip """
        formatter = ctx.make_formatter()
        self.format_help(ctx, formatter)
        return formatter.getvalue()

    def format_usage(self, ctx, formatter):
        pieces = ["[-h]", "[--min_ph MIN]", "[--max_ph MAX]", "[--pka_precision PRE]", "[--smiles SMI]",
                  "[--smiles_file FILE]", "[--output_file FILE]", "[--max_variants MXV]", "[--label_states]",
                  "[--silent]", "[--test]"]
        formatter.write_usage(ctx.command_path, " ".join(pieces))
        formatter.write_text(" ")
        formatter.write_text(HEADER)
        formatter.write_text(" \n")


@click.command(help="Dimorphite 1.2.4: Creates models of appropriately protonated small molecules.\n"
                    "Apache 2.0 License. Copyright 2020 Jacob D. Durrant.",
               epilog=examples,
               cls=SpecialEpilog, invoke_without_command=True,
               context_settings={'help_option_names': ['-h', '--help']})
@click.option("--min_ph", metavar="MIN", type=float, default=6.4, required=False,
              help="minimum pH to consider (default: 6.4)")
@click.option("--max_ph", metavar="MAX", type=float, default=8.4, required=False,
              help="maximum pH to consider (default: 8.4)")
@click.option("--pka_precision", metavar="PRE", type=float, default=1.0, required=False,
              help="pKa precision factor (number of standard devations, default: 1.0)")
@click.option("--smiles", metavar="SMI", type=str, required=True, cls=Mutex, not_required_if=['smiles_file'],
              help="SMILES string to protonate")
@click.option("--smiles_file", metavar="FILE", type=str, required=True, cls=Mutex, not_required_if=['smiles'],
              help="file that contains SMILES strings to protonate")
@click.option("--output_file", metavar="FILE", type=str, required=False,
              help="output file to write protonated SMILES (optional)")
@click.option("--max_variants", metavar="MXV", type=int, default=128, required=False,
              help="limit number of variants per input compound (default: 128)")
@click.option("--label_states", is_flag=True, required=False, default=False,
              help="label protonated SMILES with target state (i.e., \"DEPROTONATED\", \"PROTONATED\", or \"BOTH\").")
@click.option("--silent", is_flag=True, required=False, default=False,
              help="do not print any messages to the screen")
@click.option("--test", is_flag=True, required=False, default=False,
              help="run unit tests (for debugging)")
def main(min_ph: float = 6.4,
         max_ph: float = 8.4,
         pka_precision: float = 1.0,
         smiles: Optional[str] = None,
         smiles_file: Optional[str] = None,
         output_file: Optional[str] = None,
         max_variants: int = 128,
         label_states: bool = False,
         silent: bool = False,
         test: bool = False):
    if test:
        TestFuncs.test()
        return
    # Main part
    if smiles is not None:
        if isinstance(smiles, str):
            smiles_file = StringIO(smiles)
    # Read in data in SMI format
    smiles_and_data = LoadSMIFile(smiles_file, silent=silent)
    # Separate SMILES from data
    all_smis, all_data = zip(*[[smi_and_data['smiles'], smi_and_data['data']] for smi_and_data in smiles_and_data])
    all_smis, all_data = list(all_smis), list(all_data)
    # Run protonation
    protonated_smiles = run(all_smis, min_ph=min_ph, max_ph=max_ph, pka_precision=pka_precision,
                            max_variants=max_variants, label_states=label_states,
                            silent=silent)
    if label_states:
        for i in range(len(protonated_smiles)):
            if protonated_smiles[i] is not None:
                label = protonated_smiles[i][1]
                data = all_data[i]
                for smi in protonated_smiles[i][0]:
                    line = f"{smi}"
                    if len(data) > 0:
                        line += '\t' + '\t'.join(data)
                    if label is not None:
                        if isinstance(label, list):
                            line += '\t' + '\t'.join(label)
                        else:
                            line += f"\t{label}"
                    line += '\n'
                    if output_file is not None:
                        # An output file was specified, so write to that.
                        with open(output_file, "a") as oh:
                            print(line, file=oh)
                    else:
                        print(line, end='')
    else:
        for i in range(len(protonated_smiles)):
            if protonated_smiles[i] is not None:
                data = all_data[i]
                line = [protonated_smiles[i]] if isinstance(protonated_smiles[i], str)\
                    else protonated_smiles[i]
                if len(data) > 0:
                    line = '\n'.join([f"{smi}\t" + '\t'.join(data) for smi in line])
                else:
                    line = '\n'.join(line)
                line += '\n'
                if output_file is not None:
                    # An output file was specified, so write to that.
                    with open(output_file, "a") as oh:
                        print(line, file=oh)
                else:
                    print(line, end='')
