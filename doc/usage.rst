Quickstart and Basic Usage
==============================================

.. testsetup:: *

   import basis_set_exchange


Installation
-------------------

The BSE can be installed with `pip`.

To install from PyPI: `pip install basis_set_exchange`

To install from the local directory, you can use `pip install -e .` inside
the cloned git repository.


Importing
-------------------

All end-user functionality is available by importing the `basis_set_exchange` module.


Determining the library version
-------------------------------

The library version can be determined with :func:`basis_set_exchange.version()`


Getting a basis set from the library
------------------------------------

The main function for getting a basis set is :func:`basis_set_exchange.get_basis`.
Output format is controlled by the `fmt` parameter. By default, a python
dictionary is returned. If a format is specified, a string is returned
instread. This string will contain the newlines needed to properly format
the string into a basis set.

The name of the basis set is not case sensitive.

The available formats are listed at the documentation for :func:`basis_set_exchange.get_basis`
and can be obtained via :func:`basis_set_exchange.get_formats`

.. doctest::
   :pyversion: >= 3.6

   >>> # Get a basis set as a python dictionary
   >>> bs_dict = basis_set_exchange.get_basis('6-31G*')
   >>> bs_dict['name']
   '6-31G*'

   >>> # Basis set names are case insensitive
   >>> bs_dict = basis_set_exchange.get_basis('6-31g*')
   >>> bs_dict['name']
   '6-31G*'

   >>> # Same as above, but in gaussian format (as a string)
   >>> # header=False disables printing an information block
   >>> bs_str = basis_set_exchange.get_basis('6-31G*', fmt='gaussian94', header=False)
   >>> print(bs_str)
   H     0
   S    3   1.00
         0.1873113696D+02       0.3349460434D-01
         0.2825394365D+01       0.2347269535D+00
         0.6401216923D+00       0.8137573261D+00
   ...


   >>> # Available formats are available via get_formats
   >>> basis_set_exchange.get_formats()
   {'nwchem': 'NWChem', 'gaussian94': 'Gaussian', ...


By default, all elements for which the basis set is defined are included. You can specify a 
subset of elements instead by using the `elements` parameter

.. doctest::
   :pyversion: >= 3.6

   >>> # Get only carbon and oxygen
   >>> bs_str = basis_set_exchange.get_basis('aug-cc-pvtz', elements=[6,8], fmt='nwchem', header=False)
   >>> print(bs_str)
   BASIS "ao basis" SPHERICAL PRINT
   #BASIS SET: (11s,6p,3d,2f) -> [5s,4p,3d,2f]
   C    S
         8.236000E+03           5.310000E-04           0.000000E+00          -1.130000E-04           0.000000E+00
   ...

   >>> # Can also use strings with the element symbols (and be mixed with integers)
   >>> # and integers as strings
   >>> bs_str = basis_set_exchange.get_basis('aug-cc-pvtz', elements=['C', 8, 'Ne', '16'], fmt='nwchem', header=False)
   >>> print(bs_str)
   BASIS "ao basis" SPHERICAL PRINT
   #BASIS SET: (11s,6p,3d,2f) -> [5s,4p,3d,2f]
   C    S
         8.236000E+03           5.310000E-04           0.000000E+00          -1.130000E-04           0.000000E+00
   ...


Getting references
------------------

Reference/citations can be obtained via :func:`basis_set_exchange.get_references`. The `elements`
parameter is similar to that in :func:`basis_set_exchange.get_basis`.

The `fmt` parameter controls the output format. By default, the output
is a dictionary. If `fmt` is specified, the output is a string.

The available formats are listed at the documentation for :func:`basis_set_exchange.get_references`
and can be obtained via :func:`basis_set_exchange.get_reference_formats`

.. doctest::
   :pyversion: >= 3.6

   >>> # Get references for 6-31G*, all elements, as a list of dictionaries
   >>> refs = basis_set_exchange.get_references('6-31G*')
   >>> print(refs[0])
   {'reference_info': [{'reference_description': ...

   >>> # As bibtex, restricting to H and F
   >>> bib = basis_set_exchange.get_references('6-31G*', fmt='bib', elements=[1,9])
   >>> print(bib)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % If you downloaded data from the basis set
   % exchange or used the basis set exchange python library, please cite:
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ...
   % H
   %     31G Split-valence basis set for H,He
   %         ditchfield1971a
   %
   % F
   %     6-31G Split-valence basis set
   %         hehre1972a
   %
   %     Polarization for 6-31G split-valence basis set
   %         hariharan1973a
   %
   <BLANKLINE>
   <BLANKLINE>
   @article{ditchfield1971a,
       author = {Ditchfield, R. and Hehre, W. J. and Pople, J. A.},
       title = {Self-Consistent Molecular-Orbital Methods. IX. An Extended Gaussian-Type Basis for Molecular-Orbital Studies of Organic Molecules},
       journal = {J. Chem. Phys.},
       volume = {54},
       pages = {724-728},
       year = {1971},
       doi = {10.1063/1.1674902}
   }
   ...


   >>> # Available formats are available via get_reference_formats
   >>> basis_set_exchange.get_reference_formats()
   {'txt': 'Plain Text', 'bib': 'BibTeX', 'ris': 'RIS', 'endnote': 'EndNote', 'json': 'JSON'}


Versioning
-------------------

Basis sets within this package are versioned. This allows for changes to be made to a
basis set, while keeping the old data accessible for historical purposes.
Versions are specified by integers. By default,
v0 will match the original EMSL BSE data.

Versions are meant to be increased only when there is a material change to the data.
If data is simply being added (new elements), the version will not be incremented.

Both `basis_set_exchange.get_basis` and :func:`basis_set_exchange.get_references` accept a `version` parameter,
which is a string. If `version` is not specified, the latest version is used.

.. doctest::
   :pyversion: >= 3.6

   >>> # Get the latest version of 6-31G*
   >>> basis_set_exchange.get_basis('6-31G*', fmt='gaussian94', header=False)
   'H     0\nS    3   1.00\n      0.1873113696D+02       0.3349460434D-01\n      0.2825394365D+01...

   >>> # Get the original BSE data
   >>> basis_set_exchange.get_basis('6-31G*', version=0, fmt='gaussian94', header=False)
   'H     0\nS    3   1.00\n     18.7311370              0.03349460\n      2.8253937...

   >>> # Versions can also be passed as strings
   >>> basis_set_exchange.get_basis('6-31G*', version='0', fmt='gaussian94', header=False)
   'H     0\nS    3   1.00\n     18.7311370              0.03349460\n      2.8253937...



Listing Basis Sets and Getting Basis Set Metadata
-------------------------------------------------

This package contains metadata for all the basis sets that is in its data directory.
This information can be accessed by the :func:`basis_set_exchange.get_metadata` function

.. note:: Note that the key is the name of the basis set that has been transformed
          into some internal name (see :func:`basis_set_exchange.transform_basis_name`)

A simple list containing all the basis set names can be obtained via :func:`basis_set_exchange.get_all_basis_names`.
A simple list of families can be obtained with :func:`basis_set_exchange.get_families`.

.. doctest::
   :pyversion: >= 3.6

   >>> # Get the metadata
   >>> md = basis_set_exchange.get_metadata()

   >>> # What is the latest version of 6-31G
   >>> md['6-31g']['latest_version']
   '1'

   >>> # All versions of 6-31G
   >>> md['6-31g']['versions'].keys()
   dict_keys(['0', '1'])

   >>> # Elements defined in v0
   >>> md['6-31g']['versions']['0']['elements']
   ['1', '2', '3', '4', '5', '6',...

   >>> # Print all the basis sets known to the BSE
   >>> basis_set_exchange.get_all_basis_names()
   ['2ZaPa-NR', '2ZaPa-NR-CV', '3-21G', '3ZaPa-NR', '3ZaPa-NR-CV',...

   >>> # A list of all families
   >>> basis_set_exchange.get_families()
   ['acvxz-j', 'ahlrichs', 'ahlrichs_dhf', 'ahlrichs_fit', ...


Lookup a basis by Role
----------------------

Many basis sets have auxiliary basis sets for different purposes (density fitting,
for example). These auxiliary basis sets can be queried in the BSE
using the :func:`basis_set_exchange.lookup_basis_by_role`. This function takes the
primary basis set and the role you wish to look up. The function
returns the name of the basis set.

Like the other functions, the basis name and role are not case sensitive.

The available roles are listed at the documentation for :func:`basis_set_exchange.lookup_basis_by_role`
and can be obtained via :func:`basis_set_exchange.get_roles`

.. doctest::
   :pyversion: >= 3.6

   >>> # Find the MP2-fit basis set for cc-pvtz
   >>> basis_set_exchange.lookup_basis_by_role('cc-pvtz', 'rifit')
   'cc-pvtz-rifit'

   >>> # Find the J-fit basis set for def2-TZVP
   >>> basis_set_exchange.lookup_basis_by_role('def2-tzvp', 'jfit')
   'def2-universal-jfit'

   >>> # Available roles are available via get_roles
   >>> basis_set_exchange.get_roles()
   {'orbital': 'Orbital basis', 'jfit': 'J-fitting', 'jkfit': 'JK-fitting', 'rifit': 'RI-fitting',...



Filtering basis sets
--------------------

Basis sets can be searched for via simple filtering with :func:`basis_set_exchange.filter_basis_sets`. All
search parameters are case insensitive. Basis sets match if all criteria are true.
   
.. doctest::
   :pyversion: >= 3.6

   >>> # Find all basis sets with '31g' in the name
   >>> md = basis_set_exchange.filter_basis_sets('31g')
   >>> md.keys()
   dict_keys(['4-31g', ...

   >>> # Find all basis sets with 'aug' in the name that can be used for RI fitting
   >>> md = basis_set_exchange.filter_basis_sets('aug', role='rifit')
   >>> md.keys()
   dict_keys(['aug-cc-pv5z-pp-rifit', 'aug-cc-pv5z-rifit', 'aug-cc-pv6z-rifit', ...

   >>> # All basis sets of the dunning family that have '5z' in the name
   >>> md = basis_set_exchange.filter_basis_sets('5z', family='dunning')
   >>> md.keys()
   dict_keys(['aug-cc-pcv5z', 'aug-cc-pv5z', 'aug-cc-pwcv5z', 'cc-pcv5z', 'cc-pv5z', 'cc-pwcv5z'...


Basis set and family notes
--------------------------

Notes about a basis set or a basis set family can be obtained, also.

.. doctest::

   >>> # Notes from a basis (name is case insensitive)
   >>> basis_set_exchange.get_basis_notes('6-31g')
   '--------------------------------------------------------------------------------\n   Original BSE Contributor: Dr. David Feller...

   >>> # Get the family of a basis set from the metadata
   >>> fam = basis_set_exchange.get_basis_family('6-31G**')
   >>> fam
   'pople'

   >>> # Get family notes (not case sensitive)
   >>> basis_set_exchange.get_family_notes('pople')
   'Notes about Pople basis sets...


Memoization
-----------

By default, the library will memoize/cache some internal data. This has a big effect when,
for example, running :func:`basis_set_exchange.get_basis` with the same basis set name (even if choosing
different elements and options).

For most uses, this can be left enabled - memory usage will still be very low, even if reading
many basis sets. If you wish, it can be disabled by setting :attr:`basis_set_exchange.memo.memoize_enabled` to `False`.
Note that this does not clear any existing cache.


   >>> # Default is enabled
   >>> basis_set_exchange.memo.memoize_enabled
   True

   >>> # Manually disable it
   >>> basis_set_exchange.memo.memoize_enabled = False
   >>> basis_set_exchange.memo.memoize_enabled
   False
