'''
Handlers for command line subcommands
'''

from .. import api, bundle, readers, writers, refconverters, convert
from ..misc import compact_elements
from .common import format_columns


def _bse_cli_list_basis_sets(args):
    '''Handles the list-basis-sets subcommand'''
    metadata = api.filter_basis_sets(args.substr, args.family, args.role, args.elements, args.data_dir)

    if args.no_description:
        liststr = [x['display_name'] for x in metadata.values()]
    else:
        liststr = format_columns([(v['display_name'], v['description']) for k, v in metadata.items()])

    return '\n'.join(liststr)


def _bse_cli_list_families(args):
    '''Handles the list-families subcommand'''
    families = api.get_families(args.data_dir)
    return '\n'.join(families)


def _bse_cli_list_writer_formats(args):
    '''Handles the list-writer-formats subcommand'''
    all_formats = writers.get_writer_formats()

    if args.no_description:
        liststr = all_formats.keys()
    else:
        liststr = format_columns(all_formats.items())

    return '\n'.join(sorted(liststr))


def _bse_cli_list_reader_formats(args):
    all_formats = readers.get_reader_formats()

    if args.no_description:
        liststr = all_formats.keys()
    else:
        liststr = format_columns(all_formats.items())

    return '\n'.join(liststr)


def _bse_cli_list_formats(args):
    all_formats = api.get_formats()

    if args.no_description:
        liststr = all_formats.keys()
    else:
        liststr = format_columns(all_formats.items())

    return '\n'.join(liststr)


def _bse_cli_list_ref_formats(args):
    '''Handles the list-ref-formats subcommand'''
    all_refformats = refconverters.get_reference_formats()

    if args.no_description:
        liststr = all_refformats.keys()
    else:
        liststr = format_columns(all_refformats.items())

    return '\n'.join(liststr)


def _bse_cli_list_roles(args):
    '''Handles the list-roles subcommand'''
    all_roles = api.get_roles()

    if args.no_description:
        liststr = all_roles.keys()
    else:
        liststr = format_columns(all_roles.items())

    return '\n'.join(liststr)


def _bse_cli_get_data_dir(args):
    '''Handles the get-data-dir subcommand'''

    return api.get_data_dir()


def _bse_cli_lookup_by_role(args):
    '''Handles the lookup-by-role subcommand'''
    return api.lookup_basis_by_role(args.basis, args.role, args.data_dir)


def _bse_cli_get_basis(args):
    '''Handles the get-basis subcommand'''

    return api.get_basis(name=args.basis,
                         elements=args.elements,
                         version=args.version,
                         fmt=args.fmt,
                         uncontract_general=args.unc_gen,
                         uncontract_spdf=args.unc_spdf,
                         uncontract_segmented=args.unc_seg,
                         remove_free_primitives=args.rm_free,
                         make_general=args.make_gen,
                         optimize_general=args.opt_gen,
                         augment_diffuse=args.aug_diffuse,
                         augment_steep=args.aug_steep,
                         get_aux=args.get_aux,
                         data_dir=args.data_dir,
                         header=not args.noheader)


def _bse_cli_get_refs(args):
    '''Handles the get-refs subcommand'''
    return api.get_references(basis_name=args.basis,
                              elements=args.elements,
                              version=args.version,
                              fmt=args.reffmt,
                              data_dir=args.data_dir)


def _bse_cli_get_info(args):
    '''Handles the get-info subcommand'''

    bs_meta = api.get_metadata(args.data_dir)[args.basis]
    ret = []
    ret.append('-' * 80)
    ret.append(args.basis)
    ret.append('-' * 80)
    ret.append('    Display Name: ' + bs_meta['display_name'])
    ret.append('     Description: ' + bs_meta['description'])
    ret.append('            Role: ' + bs_meta['role'])
    ret.append('          Family: ' + bs_meta['family'])
    ret.append('  Function Types: ' + ','.join(bs_meta['function_types']))
    ret.append('  Latest Version: ' + bs_meta['latest_version'])
    ret.append('')

    aux = bs_meta['auxiliaries']
    if not aux:
        ret.append('Auxiliary Basis Sets: None')
    else:
        ret.append('Auxiliary Basis Sets:')
        ret.extend(format_columns(list(aux.items()), '    '))

    ver = bs_meta['versions']
    ret.append('')
    ret.append('Versions:')

    # Print 4 columns - version, date, elements, revision description
    version_lines = format_columns([(k, v['revdate'], compact_elements(v['elements']), v['revdesc'])
                                    for k, v in ver.items()], '    ')
    ret.extend(version_lines)

    return '\n'.join(ret)


def _bse_cli_get_notes(args):
    '''Handles the get-notes subcommand'''
    return api.get_basis_notes(args.basis, args.data_dir)


def _bse_cli_get_family(args):
    '''Handles the get-family subcommand'''
    return api.get_basis_family(args.basis, args.data_dir)


def _bse_cli_get_versions(args):
    '''Handles the get-versions subcommand'''
    name = args.basis.lower()
    metadata = api.get_metadata(args.data_dir)
    if name not in metadata:
        raise KeyError(
            "Basis set {} does not exist. For a complete list of basis sets, use the 'list-basis-sets' command".format(
                name))

    version_data = {k: v['revdesc'] for k, v in metadata[name]['versions'].items()}

    if args.no_description:
        liststr = version_data.keys()
    else:
        liststr = format_columns(version_data.items())

    return '\n'.join(liststr)


def _bse_cli_get_family_notes(args):
    '''Handles the get-family-notes subcommand'''
    return api.get_family_notes(args.family, args.data_dir)


def _bse_cli_convert_basis(args):
    '''Handles the convert-basis subcommand'''

    # We convert file -> file
    convert.convert_formatted_basis_file(args.input_file, args.output_file, args.in_fmt, args.out_fmt)
    return "Converted {} -> {}".format(args.input_file, args.output_file)


def _bse_cli_create_bundle(args):
    '''Handles the create-bundle subcommand'''
    bundle.create_bundle(args.bundle_file, args.fmt, args.reffmt, args.archive_type, args.data_dir)
    return "Created " + args.bundle_file


def bse_cli_handle_subcmd(args):
    handler_map = {
        'list-formats': _bse_cli_list_writer_formats,
        'list-writer-formats': _bse_cli_list_writer_formats,
        'list-reader-formats': _bse_cli_list_reader_formats,
        'list-ref-formats': _bse_cli_list_ref_formats,
        'list-roles': _bse_cli_list_roles,
        'get-data-dir': _bse_cli_get_data_dir,
        'list-basis-sets': _bse_cli_list_basis_sets,
        'list-families': _bse_cli_list_families,
        'lookup-by-role': _bse_cli_lookup_by_role,
        'get-basis': _bse_cli_get_basis,
        'get-refs': _bse_cli_get_refs,
        'get-info': _bse_cli_get_info,
        'get-notes': _bse_cli_get_notes,
        'get-family': _bse_cli_get_family,
        'get-versions': _bse_cli_get_versions,
        'get-family-notes': _bse_cli_get_family_notes,
        'convert-basis': _bse_cli_convert_basis,
        'create-bundle': _bse_cli_create_bundle
    }

    return handler_map[args.subcmd](args)
