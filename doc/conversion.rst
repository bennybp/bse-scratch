.. _conversion:

Converting basis sets
==============================================


About conversion
----------------

This library (and command-line interface) has functionality for converting
basis sets from one format to another.

If you have a file that does not convert properly, open an issue
(https://github.com/MolSSI-BSE/basis_set_exchange/issues)

Not all formats that can be written by the library can be read by the library.
Functionality for reading and writing are developed separately.

When converting files, the format is attempted to be discovered based on the filename.
If it cannot be detected, an exception is raised.

This functionality can also read/write bzip2-compressed files seamlessly.


Conversion via python
---------------------

Basis sets can be converted using :func:`convert.convert_formatted_basis_file` and
:func:`convert.convert_formatted_basis_str`. The former converts data stored in a plain-text
file, while the latter converts a (python) string.

An alternate location of data to export can be specified with `data_dir`::

   >>> # Convert /tmp/test.nw from NWChem to Gaussian format
   >>> bse.convert_formatted_basis_file('/path/test.nw', '/path/test.gbs')

   >>> # Manually specify the formats
   >>> bse.convert_formatted_basis_file('/path/test_nwchem.txt', '/path/test_gaussian.txt', in_fmt='nwchem', out_fmt='gaussian94')



Conversion via the command line
-------------------------------


The same functionality can be accessed via the command-line interface::

   bse convert-basis /path/basis.nw /path/basis.gbs

Formats can be specified explicitly if autodetection fails::

   bse convert-basis --in-fmt nwchem --out-fmt gaussian94 /path/file-nwchem.txt /path/file-gaussian.txt
