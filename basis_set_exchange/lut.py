'''
Functions for looking up element names, numbers, and symbols
and for converting between angular momentum conventions

This module has functions for looking up element information by
symbol, name, or number. It also has functions for converting
angular momentum between integers (0, 1, 2) and letters (s, p, d).
'''

# Contains the symbols, names, and Z numbers for all known elements
# NOTE: Some Z number have multiple entries. This is to allow for querying
#       based on older, systematic names/symbols (like uuh), however a query
#       of the Z number should return the 'new', official name (livermorium or lv)
#       This also allows for different spellings of some names (cesium/caesium, aluminum/aluminium)
# yapf: disable
_data_table = [
    ("h",     1, "hydrogen"),      ("he",    2, "helium"),      ("li",    3, "lithium"),
    ("be",    4, "beryllium"),     ("b",     5, "boron"),       ("c",     6, "carbon"),
    ("n",     7, "nitrogen"),      ("o",     8, "oxygen"),      ("f",     9, "fluorine"),
    ("ne",   10, "neon"),          ("na",   11, "sodium"),      ("mg",   12, "magnesium"),
    ("al",   13, "aluminum"),      ("al",   13, "aluminium"),   ("si",   14, "silicon"),
    ("p",    15, "phosphorus"),    ("s",    16, "sulfur"),      ("cl",   17, "chlorine"),
    ("ar",   18, "argon"),         ("k",    19, "potassium"),   ("ca",   20, "calcium"),
    ("sc",   21, "scandium"),      ("ti",   22, "titanium"),    ("v",    23, "vanadium"),
    ("cr",   24, "chromium"),      ("mn",   25, "manganese"),   ("fe",   26, "iron"),
    ("co",   27, "cobalt"),        ("ni",   28, "nickel"),      ("cu",   29, "copper"),
    ("zn",   30, "zinc"),          ("ga",   31, "gallium"),     ("ge",   32, "germanium"),
    ("as",   33, "arsenic"),       ("se",   34, "selenium"),    ("br",   35, "bromine"),
    ("kr",   36, "krypton"),       ("rb",   37, "rubidium"),    ("sr",   38, "strontium"),
    ("y",    39, "yttrium"),       ("zr",   40, "zirconium"),   ("nb",   41, "niobium"),
    ("mo",   42, "molybdenum"),    ("tc",   43, "technetium"),  ("ru",   44, "ruthenium"),
    ("rh",   45, "rhodium"),       ("pd",   46, "palladium"),   ("ag",   47, "silver"),
    ("cd",   48, "cadmium"),       ("in",   49, "indium"),      ("sn",   50, "tin"),
    ("sb",   51, "antimony"),      ("te",   52, "tellurium"),   ("i",    53, "iodine"),
    ("xe",   54, "xenon"),         ("cs",   55, "cesium"),      ("cs",   55, "caesium"),
    ("ba",   56, "barium"),        ("la",   57, "lanthanum"),   ("ce",   58, "cerium"),
    ("pr",   59, "praseodymium"),  ("nd",   60, "neodymium"),   ("pm",   61, "promethium"),
    ("sm",   62, "samarium"),      ("eu",   63, "europium"),    ("gd",   64, "gadolinium"),
    ("tb",   65, "terbium"),       ("dy",   66, "dysprosium"),  ("ho",   67, "holmium"),
    ("er",   68, "erbium"),        ("tm",   69, "thulium"),     ("yb",   70, "ytterbium"),
    ("lu",   71, "lutetium"),      ("hf",   72, "hafnium"),     ("ta",   73, "tantalum"),
    ("w",    74, "tungsten"),      ("re",   75, "rhenium"),     ("os",   76, "osmium"),
    ("ir",   77, "iridium"),       ("pt",   78, "platinum"),    ("au",   79, "gold"),
    ("hg",   80, "mercury"),       ("tl",   81, "thallium"),    ("pb",   82, "lead"),
    ("bi",   83, "bismuth"),       ("po",   84, "polonium"),    ("at",   85, "astatine"),
    ("rn",   86, "radon"),         ("fr",   87, "francium"),    ("ra",   88, "radium"),
    ("ac",   89, "actinium"),      ("th",   90, "thorium"),     ("pa",   91, "protactinium"),
    ("u",    92, "uranium"),       ("np",   93, "neptunium"),   ("pu",   94, "plutonium"),
    ("am",   95, "americium"),     ("cm",   96, "curium"),      ("bk",   97, "berkelium"),
    ("cf",   98, "californium"),   ("es",   99, "einsteinium"), ("fm",  100, "fermium"),
    ("md",  101, "mendelevium"),   ("no",  102, "nobelium"),    ("lr",  103, "lawrencium"),
    ("rf",  104, "rutherfordium"), ("ha",  105, "hahnium"),     ("db",  105, "dubnium"),
    ("sg",  106, "seaborgium"),    ("bh",  107, "bohrium"),     ("hs",  108, "hassium"),
    ("mt",  109, "meitnerium"),    ("uun", 110, "ununnilium"),  ("ds",  110, "darmstadtium"),
    ("uuu", 111, "unununium"),     ("rg",  111, "roentgenium"), ("uub", 112, "ununbium"),
    ("cn",  112, "copernicium"),   ("uut", 113, "ununtrium"),   ("nh",  113, "nihonium"),
    ("uuq", 114, "ununquadium"),   ("fl",  114, "flerovium"),   ("uup", 115, "ununpentium"),
    ("mc",  115, "moscovium"),     ("uuh", 116, "ununhexium"),  ("lv",  116, "livermorium"),
    ("uus", 117, "ununseptium"),   ("ts",  117, "tennessine"),  ("uuo", 118, "ununoctium"),
    ("og",  118, "oganesson"),     ("uue", 119, "ununennium"),  ("ubn", 120, "unbinilium")
]
# yapf: enable

# Maps Z to element data
# Note that only the last entry for a particular Z number will
# be stored (which contains the new, official data rather than
# the older systematic names for some elements)
_element_Z_map = {x[1]: x for x in _data_table}

# Maps element symbol to element data
_element_sym_map = {x[0]: x for x in _data_table}

# Maps element name to element data
_element_name_map = {x[2]: x for x in _data_table}

# Maps AM characters to integers (the integer is the
# index of the character in this string)
# This is the 'hik' convention, where AM=7 is k
_amchar_map_hik = 'spdfghiklmnoqrtuvwxyzabce'

# This is the 'hij' convention, where AM=7 is j
_amchar_map_hij = 'spdfghijklmnoqrtuvwxyzabce'


def all_element_names():
    '''Obtain a list of the names of all the elements'''
    return list(_element_name_map.keys())


def element_data_from_Z(Z):
    '''Obtain elemental data given a Z number

    An exception is thrown if the Z number is not found
    '''

    # Z may be a str
    if isinstance(Z, str) and Z.isdecimal():
        Z = int(Z)

    if Z not in _element_Z_map:
        raise KeyError('No element data for Z = {}'.format(Z))
    return _element_Z_map[Z]


def element_data_from_sym(sym):
    '''Obtain elemental data given an elemental symbol

    The given symbol is not case sensitive

    An exception is thrown if the symbol is not found
    '''

    sym_lower = sym.lower()
    if sym_lower not in _element_sym_map:
        raise KeyError('No element data for symbol \'{}\''.format(sym))
    return _element_sym_map[sym_lower]


def element_data_from_name(name):
    '''Obtain elemental data given an elemental name

    The given name is not case sensitive

    An exception is thrown if the name is not found
    '''

    name_lower = name.lower()
    if name_lower not in _element_name_map:
        raise KeyError('No element data for name \'{}\''.format(name))
    return _element_name_map[name_lower]


def element_name_from_Z(Z, normalize=False):
    '''Obtain an element's name from its Z number

    An exception is thrown if the Z number is not found

    If normalize is True, the first letter will be capitalized
    '''

    r = element_data_from_Z(Z)[2]
    if normalize:
        return r.capitalize()
    else:
        return r


def element_Z_from_name(name, as_str=False):
    '''Obtain an element's Z number given its name

    If as_str is True, then a string is returned (ie, '1' for Hydrogen)

    An exception is thrown if the name is not found
    '''

    Z = element_data_from_name(name)[1]
    if as_str:
        Z = str(Z)
    return Z


def element_sym_from_Z(Z, normalize=False):
    '''Obtain an element's symbol from its Z number

    An exception is thrown if the Z number is not found

    If normalize is True, the first letter will be capitalized
    '''

    r = element_data_from_Z(Z)[0]
    if normalize:
        return r.capitalize()
    else:
        return r


def element_Z_from_sym(sym, as_str=False):
    '''Obtain an element's Z-number given its symbol

    If as_str is True, then a string is returned (ie, '1' for Hydrogen)

    An exception is thrown if the symbol is not found
    '''

    Z = element_data_from_sym(sym)[1]
    if as_str:
        Z = str(Z)
    return Z


def amint_to_char(am, hij=False, use_L=False):
    '''Convert an angular momentum integer to a character

    The input is a list (to handle sp, spd, ... orbitals). The return
    value is a string

    For example, converts [0] to 's' and [0,1,2] to 'spd'

    If hij is True, the ordering spdfghijkl is used. Otherwise, the
    ordering will be spdfghikl (skipping j)

    If use_L is True, sp shells ([0,1]) will return l instead
    '''

    if use_L and am == [0, 1]:
        return 'l'
    if hij:
        amchar_map = _amchar_map_hij
    else:
        amchar_map = _amchar_map_hik

    amchar = []

    for a in am:
        if a < 0:
            raise IndexError('Angular momentum must be a positive integer (not {})'.format(a))
        if a >= len(amchar_map):
            raise IndexError('Angular momentum {} out of range. Must be less than {}'.format(a, len(amchar_map)))
        amchar.append(amchar_map[a])

    return ''.join(amchar)


def amchar_to_int(amchar, hij=False):
    '''Convert an angular momentum integer to a character

    The return value is a list of integers (to handle sp, spd, ... orbitals)

    For example, converts 'p' to [1] and 'sp' to [0,1]

    If hij is True, the ordering spdfghijkl is used. Otherwise, the
    ordering will be spdfghikl (skipping j)
    '''

    if hij:
        amchar_map = _amchar_map_hij
    else:
        amchar_map = _amchar_map_hik

    amchar_lower = amchar.lower()

    amint = []

    for c in amchar_lower:
        if c not in amchar_map:
            raise KeyError('Angular momentum character {} is not valid'.format(c))

        amint.append(amchar_map.index(c))

    return amint


def electron_shells_start(nelectrons, max_am=20):
    '''Return the starting principle quantum numbers of electron shells

    For example, an ECP covering 10 electrons will covers 1s, 2s, 2p shells. The
    electrons shells will then start at 3s, 3p, 3d, and 4f (returned as [3, 3, 3, 4])

    If an ECP covers 30 electrons, then the shells will start at [5, 4, 4, 4]

    Only fully-covered shells are counted. If a shell is partly covered, an exception
    is raised.

    The returned list is extended up to max_am.

    Note: Since the main use of this is for ECPs, we only cover what can really be found
    on the periodic table. No excited states!

    Parameters
    ----------
    nelectrons : int
        Number of electrons covered by an ECP
    max_am : int
        Fill out the starting principal quantum numbers up to this am

    Returns
    -------
    list
        The starting principal quantum numbers of s, p, d, and f shells.
    '''

    if nelectrons < 0:
        raise RuntimeError("Can't have a negative number of electrons")
    if nelectrons > 118:
        raise NotImplementedError("Too many electrons for electrons_shells_start")

    # yapf: disable
    # The usual filling order of electrons you learned in high school
    # Tuple of (am, nelec)
    aminfo = ( (0, 2),                           # He
               (0, 2), (1, 6),                   # Ne
               (0, 2), (1, 6),                   # Ar
               (0, 2), (2, 10), (1, 6),          # Kr
               (0, 2), (2, 10), (1, 6),          # Xe
               (0, 2), (3, 14), (2, 10), (1, 6), # Rn
               (0, 2), (3, 14), (2, 10), (1, 6)  # Og
             )

    # But the above doesn't always work - ECPs don't always
    # cover things in such a regular order. So handle some
    # special ones
    special_am = { 28: [0, 0, 0, 1, 1, 2],
                   46: [0, 0, 0, 0, 1, 1, 1, 2, 2],
                   60: [0, 0, 0, 0, 1, 1, 1, 2, 2, 3],
                   68: [0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3],
                   78: [0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3],
                   92: [0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3]
                 }
    # yapf: enable

    if nelectrons in special_am:
        contained_am = special_am[nelectrons]
    else:
        contained_am = []

        for am, n in aminfo:
            if nelectrons >= n:
                nelectrons -= n
                contained_am.append(am)
            else:
                break

        if nelectrons != 0:
            raise RuntimeError("Electrons cover a partial shell. {} electrons left".format(nelectrons))

    start = [
        contained_am.count(0) + 1,
        contained_am.count(1) + 2,
        contained_am.count(2) + 3,
        contained_am.count(3) + 4
    ]
    start.extend(range(5, max_am + 2))

    return start


def function_type_from_am(shell_am, base_type, spherical_type):
    if max(shell_am) <= 1:
        return base_type
    else:
        return base_type + '_' + spherical_type
