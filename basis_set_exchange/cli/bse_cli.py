'''
Command line interface for the basis set exchange
'''

import argparse
import argcomplete
from basis_set_exchange import version
from basis_set_exchange.cli.bse_handlers import bse_cli_handle_subcmd
from basis_set_exchange.cli.check import cli_check_normalize_args
from basis_set_exchange.cli.complete import (cli_case_insensitive_validator, cli_family_completer, cli_role_completer,
                                             cli_bsname_completer, cli_write_fmt_completer, cli_read_fmt_completer,
                                             cli_reffmt_completer)


def run_bse_cli():
    ################################################################################################
    # NOTE: I am deliberately not using the 'choices' argument in add_argument. I could use it
    # for formats, etc, however I wouldn't want to use it for basis set names. Therefore, I handle
    # all of that manually so that error output is consistent and clean
    ################################################################################################

    ########################################
    # Main global options
    ########################################
    parser = argparse.ArgumentParser(description='Description of your program')
    parser.add_argument('-V', action='version', version='basis_set_exchange ' + version())
    parser.add_argument('-d', '--data-dir', metavar='PATH', help='Override which data directory to use')
    parser.add_argument('-o', '--output', metavar='PATH', help='Output to given file rather than stdout')

    subparsers = parser.add_subparsers(metavar='subcommand', dest='subcmd')
    subparsers.required = True  # https://bugs.python.org/issue9253#msg186387

    ########################################
    # Listing of data-independent info
    ########################################
    # list-formats subcommand
    subp = subparsers.add_parser('list-formats',
                                 help='Output a list of basis set formats that can be used with obtaining a basis set')
    subp.add_argument('-n', '--no-description', action='store_true', help='Print only the format names')

    # list-writer-formats subcommand
    subp = subparsers.add_parser('list-writer-formats',
                                 help='Output a list available basis set formats that can be written')
    subp.add_argument('-n', '--no-description', action='store_true', help='Print only the format names')

    # list-reader-formats
    subp = subparsers.add_parser('list-reader-formats', help='Output a list of basis set formats that can be read')
    subp.add_argument('-n', '--no-description', action='store_true', help='Print only the format names')

    # list-ref-formats subcommand
    subp = subparsers.add_parser('list-ref-formats',
                                 help='Output a list all available reference formats and descriptions')
    subp.add_argument('-n', '--no-description', action='store_true', help='Print only the reference format names')

    # list-roles subcommand
    subp = subparsers.add_parser('list-roles', help='Output a list all available roles and descriptions')
    subp.add_argument('-n', '--no-description', action='store_true', help='Print only the role names')

    ########################################
    # Listing of general info and metadata
    ########################################
    # get-data-dir
    subparsers.add_parser('get-data-dir', help='Output the default data directory of this package')

    # list-basis-sets subcommand
    subp = subparsers.add_parser('list-basis-sets', help='Output a list all available basis sets and descriptions')
    subp.add_argument('-n', '--no-description', action='store_true', help='Print only the basis set names')
    subp.add_argument('-f', '--family',
                      help='Limit the basis set list to only the specified family').completer = cli_family_completer
    subp.add_argument('-r', '--role',
                      help='Limit the basis set list to only the specified role').completer = cli_role_completer
    subp.add_argument('-s',
                      '--substr',
                      help='Limit the basis set list to only basis sets whose name contains the specified substring')
    subp.add_argument('-e',
                      '--elements',
                      help='Limit the basis set list to only basis sets that contain all the given elements')

    # list-families subcommand
    subparsers.add_parser('list-families', help='Output a list all available basis set families')

    # lookup-by-role
    subp = subparsers.add_parser('lookup-by-role', help='Lookup a companion/auxiliary basis by primary basis and role')
    subp.add_argument(
        'basis', help='Name of the primary basis we want the auxiliary basis for').completer = cli_bsname_completer
    subp.add_argument('role', help='Role of the auxiliary basis to look for').completer = cli_role_completer

    #################################
    # Output of info
    #################################
    # get-basis subcommand
    subp = subparsers.add_parser('get-basis', help='Output a formatted basis set')
    subp.add_argument('basis', help='Name of the basis set to output').completer = cli_bsname_completer
    subp.add_argument('fmt', help='Which format to output the basis set as').completer = cli_write_fmt_completer
    subp.add_argument('--elements',
                      help='Which elements of the basis set to output. Default is all defined in the given basis')
    subp.add_argument('--version', help='Which version of the basis set to output. Default is the latest version')
    subp.add_argument('--noheader', action='store_true', help='Do not output the header at the top')
    subp.add_argument('--unc-gen', action='store_true', help='Remove general contractions')
    subp.add_argument('--unc-spdf', action='store_true', help='Remove combined sp, spd, ... contractions')
    subp.add_argument('--unc-seg', action='store_true', help='Remove general contractions')
    subp.add_argument('--rm-free', action='store_true', help='Remove free primitives')
    subp.add_argument('--opt-gen', action='store_true', help='Optimize general contractions')
    subp.add_argument('--make-gen', action='store_true', help='Make the basis set as generally-contracted as possible')
    subp.add_argument('--aug-steep', type=int, default=0, help='Augment with n steep functions')
    subp.add_argument('--aug-diffuse', type=int, default=0, help='Augment with n diffuse functions')
    subp.add_argument('--get-aux',
                      type=int,
                      default=0,
                      help='Instead of the orbital basis, get an automatically formed auxiliary basis')

    # get-refs subcommand
    subp = subparsers.add_parser('get-refs', help='Output references for a basis set')
    subp.add_argument('basis',
                      help='Name of the basis set to output the references for').completer = cli_bsname_completer
    subp.add_argument('reffmt', help='Which format to output the references as').completer = cli_reffmt_completer
    subp.add_argument('--elements',
                      help='Which elements to output the references for. Default is all defined in the given basis.')
    subp.add_argument('--version', help='Which version of the basis set to get the references for')

    # get-info subcommand
    subp = subparsers.add_parser('get-info', help='Output general info and metadata for a basis set')
    subp.add_argument('basis', help='Name of the basis set to output the info for').completer = cli_bsname_completer

    # get-notes subcommand
    subp = subparsers.add_parser('get-notes', help='Output the notes for a basis set')
    subp.add_argument('basis', help='Name of the basis set to output the notes for').completer = cli_bsname_completer

    # get-family subcommand
    subp = subparsers.add_parser('get-family', help='Output the family of a basis set')
    subp.add_argument('basis', help='Name of the basis set to output the family for').completer = cli_bsname_completer

    # get-versions subcommand
    subp = subparsers.add_parser('get-versions', help='Output a list all available versions of a basis set')
    subp.add_argument('basis', help='Name of the basis set to list the versions of').completer = cli_bsname_completer
    subp.add_argument('-n', '--no-description', action='store_true', help='Print only the version numbers')

    # get-family-notes subcommand
    subp = subparsers.add_parser('get-family-notes', help='Get the notes of a family of basis sets')
    subp.add_argument('family', type=str.lower,
                      help='The basis set family to the get the notes of').completer = cli_family_completer

    #################################
    # Converting basis sets
    #################################
    subp = subparsers.add_parser('convert-basis', help='Convert basis set files from one format to another')
    subp.add_argument('input_file', type=str, help='Basis set file to convert')
    subp.add_argument('output_file', type=str, help='Converted basis set file')
    subp.add_argument(
        '--in-fmt', type=str, default=None,
        help='Input format (default: autodetected from input filename').completer = cli_read_fmt_completer
    subp.add_argument(
        '--out-fmt', type=str, default=None,
        help='Output format (default: autodetected from output filename').completer = cli_write_fmt_completer

    #################################
    # Creating bundles
    #################################
    subp = subparsers.add_parser('create-bundle', help='Create a bundle of basis sets')
    subp.add_argument('fmt', help='Which format to output the basis set as').completer = cli_write_fmt_completer
    subp.add_argument('reffmt', help='Which format to output the references as').completer = cli_reffmt_completer
    subp.add_argument('bundle_file', help='Bundle/Archive file to create')
    subp.add_argument('--archive-type', help='Override the type of archive to create (zip or tbz)')

    #############################
    # DONE WITH SUBCOMMANDS
    #############################

    # setup autocomplete
    argcomplete.autocomplete(parser, validator=cli_case_insensitive_validator)

    # Now parse and handle the args
    args = parser.parse_args()

    # Check and make sure basis sets, roles, etc, are valid
    args = cli_check_normalize_args(args)

    # Actually generate the output
    output = bse_cli_handle_subcmd(args)

    if args.output:
        with open(args.output, 'w', encoding='utf-8') as outfile:
            outfile.write(output)
    else:
        print(output)

    return 0


if __name__ == "__main__":
    run_bse_cli()
