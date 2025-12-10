from placeholder_name.utils.constants import AMINO_ACID_SUB_SITES, PROTECTING_GROUPS, SPECIAL_CASES, AA_FULL, N_CAPS, C_CAPS, COUNTER_ACIDS, GREEK_LETTERS, PREFIX_MAP
from placeholder_name.resolvers.opsin_resolver import name_to_smiles_opsin
from placeholder_name.utils.logging_config import logger
from typing import List, Dict, Tuple


def generate_side_chain_protections():
    """Generate the complete side chain protections dictionary"""
    protections = {}
    
    # Generate standard combinations
    for aa_code, (aa_name, site) in AA_FULL.items():
        for pg_code, pg_name in PROTECTING_GROUPS.items():
            key = f"{aa_code}({pg_code})"
            
            # Handle special formatting for different types
            if pg_code.startswith('o'):  # Ester forms
                value = f"{site}-{pg_name}"
            elif aa_code in ['asp', 'glu'] and not pg_code.startswith('o'):
                # For asp/glu without 'o' prefix, assume ester
                value = f"{site}-{pg_name} ester"
            elif aa_code == 'met' and pg_code == 'o':
                value = f"{site}-{pg_name}"
            else:
                value = f"{site}-{pg_name}"
            
            protections[key] = (aa_name, value)
    
    # Add special cases
    protections.update(SPECIAL_CASES)
    
    return protections


# Generate the complete dictionary
SIDE_CHAIN_PROTECTIONS = generate_side_chain_protections()


def split_peptide_shorthand(shorthand: str) -> List[str]:
    """
    Split a peptide shorthand string into individual tokens.

    Tokens are separated by hyphens (-) unless they are inside parentheses, in which case they are treated as a single token.

    Args:
        shorthand (str): A peptide shorthand string

    Returns:
        List of tokens
    """
    tokens = []
    current_token = ""
    paren_depth = 0
    
    for char in shorthand.strip():
        if char == '(':
            paren_depth += 1
            current_token += char
        elif char == ')':
            paren_depth -= 1
            current_token += char
        elif char == '-' and paren_depth == 0:
            # Split here - we're not inside parentheses
            if current_token:
                tokens.append(current_token)
                current_token = ""
        else:
            current_token += char
    
    if current_token:
        tokens.append(current_token)
    
    return tokens


def parse_protected_residue(token: str):
    """
    Parse a protected residue token into its base amino acid and protection information.

    Args:
        token (str): The protected residue token to parse.

    Returns:
        tuple: (base_aa, protection) or (None, None) if the token cannot be parsed.
    """
    if '(' in token and token.endswith(')'):
        if token.lower() in SIDE_CHAIN_PROTECTIONS:
            base_aa, protection = SIDE_CHAIN_PROTECTIONS[token.lower()]
            return base_aa, protection
        else:
            # Try to parse manually
            base = token.split('(')[0]
            protection_part = token[token.find('(')+1:-1]
            if base in AA_FULL:
                return AA_FULL[base][0], f"with {protection_part} protection"
    return None, None


def extract_prefix_modifier(token: str):
    """
    Extract prefix modifiers from amino acid tokens.
    
    Args:
        token (str): The amino acid token to process
    
    Returns:
        tuple: (modifier_string, remaining_token) or (None, original_token)
    """
    token_lower = token.lower().strip().replace(' ', '')
    
    # Check for parenthesized prefixes first (longer matches)
    for prefix, full_name in sorted(PREFIX_MAP.items(), key=len, reverse=True):
        if token_lower.startswith(prefix):
            remaining = token[len(prefix):]
            remaining_lower = remaining.lower().strip().replace(' ', '')
            
            # Verify the remaining part is a valid amino acid
            if remaining_lower in AA_FULL:
                return full_name, remaining, remaining_lower
    
    return None, token, token.lower().strip().replace(' ', '')


def process_amino_acid_token(
    next_token: str, 
    is_last: bool, 
    is_cyclic: bool, 
    stereo_prefix: str = 'l-'
) -> str:
    """
    Process an amino acid token from a peptide shorthand and return its IUPAC name.

    Args:
        next_token (str): The amino acid token to process.
        is_last (bool): Whether this is the last token in the peptide.
        is_cyclic (bool): Whether the peptide is cyclic.
        stereo_prefix (str, optional): The stereochemistry prefix to add to the amino acid name. Defaults to 'l-'.

    Returns:
        str: The IUPAC name of the amino acid.
    """
    next_token_lower_stripped = next_token.lower().strip().replace(' ', '')
    
    # Handle Greek letter prefix
    greek_prefix = None
    if next_token and next_token[0] in GREEK_LETTERS:
        greek_prefix = next_token[0]
        next_token = next_token[1:]
        next_token_lower_stripped = next_token.lower().strip().replace(' ', '')

    # Handle prefix modifiers (methyl, acetyl, etc.)
    prefix_modifier, next_token, next_token_lower_stripped = extract_prefix_modifier(next_token)

    is_last_updated = is_last and not is_cyclic
    
    # Try to parse the token
    base_aa, protection = parse_protected_residue(next_token_lower_stripped)
    if not base_aa and next_token_lower_stripped in AA_FULL:
        base_aa = AA_FULL[next_token_lower_stripped][0]
    
    if base_aa:
        # Known amino acid
        if is_last_updated:
            base = base_aa
        else:
            # Convert to -yl form
            base_key = None
            for key, val in AA_FULL.items():
                if val == base_aa:
                    base_key = key
                    break
            if base_key and base_key.lower().strip().replace(' ', '') in AA_FULL:
                base = AA_FULL[base_key.lower().strip().replace(' ', '')][1]
            else:
                base = base_aa.replace('ine', 'yl').replace('ic acid', 'yl')
        
        # Add stereochemistry prefix (not for glycine)
        if not base_aa.startswith('glycine'):
            base = stereo_prefix + base

        if greek_prefix:
            base = greek_prefix + '-' + base

        if prefix_modifier:
            base = prefix_modifier + base
            
        if protection:
            base = f"{protection}-{base}"
        
        return base
    else:
        # Unknown token
        return f"{stereo_prefix}{next_token}"


def peptide_shorthand_to_iupac(shorthand: str) -> str:
    """
    Converts a peptide shorthand into its IUPAC name.

    The function processes the shorthand by first stripping whitespace and dashes, then
    checking if it's a cyclic peptide. If not cyclic, it removes any counter acid
    suffixes (if present) and processes the shorthand as a regular peptide.

    Parameters:
        shorthand (str): The peptide shorthand to process

    Returns:
        str: The IUPAC name of the peptide
    """
    shorthand = shorthand.strip()
    shorthand = shorthand.strip('-')
    shorthand = shorthand.replace('--', '-')
    
    # Check if it's a cyclic peptide
    is_cyclic = False
    cyclic_patterns = [
        ('cyclo(', ')', 6, -1),
        ('cyclo[', ']', 6, -1),
        ('cyclo-(', ')', 7, -1),
        ('cyclo-[', ']', 7, -1),
        ('cyclo-{', '}', 7, -1),
    ]
    
    for start_pattern, end_pattern, start_idx, end_idx in cyclic_patterns:
        if shorthand.lower().startswith(start_pattern) and shorthand.endswith(end_pattern):
            is_cyclic = True
            shorthand = shorthand[start_idx:end_idx].strip()
            break
    
    if not is_cyclic and shorthand.lower().startswith('cyclo'):
        is_cyclic = True
        shorthand = shorthand[5:].strip()

    items_to_remove = []
    counter_acid_suffix = None
    if len(shorthand.split('.')) == 2:
        for item in shorthand.split('.'):
            if item.lower().strip().replace(' ', '') in COUNTER_ACIDS:
                counter_acid_suffix = COUNTER_ACIDS[item.lower().strip().replace(' ', '')]
                items_to_remove.append(item)
        if len(items_to_remove) == 1:
            new_shorthand_split = shorthand.split('.')
            new_shorthand_split.remove(items_to_remove[0])
            shorthand = '.'.join(new_shorthand_split)
    
    tokens = split_peptide_shorthand(shorthand.strip())

    # N-cap
    prefix = ''
    if tokens and tokens[0].lower() in N_CAPS:
        prefix = N_CAPS[tokens.pop(0).lower()]

    # C-cap  
    suffix = ''
    if tokens and tokens[-1].lower() in C_CAPS:
        suffix = C_CAPS[tokens.pop().lower()]

    if not tokens:
        raise ValueError("No residues found")

    parts = []
    i = 0
    while i < len(tokens):
        t = tokens[i]
        t_lower_stripped = t.lower().strip().replace(' ', '')
        is_last = (i == len(tokens) - 1)
        
        # For cyclic peptides, all residues should be in -yl form except we need special handling
        # since there's no true "last" residue in a cycle
        if is_cyclic:
            is_last = False  # Treat all as intermediate residues

        # Check if this token is a stereochemistry indicator
        if t_lower_stripped in ['d', '(d)']:
            i += 1
            if i >= len(tokens):
                parts.append('d')
                break
            is_last = (i == len(tokens) - 1)
            parts.append(process_amino_acid_token(tokens[i], is_last, is_cyclic, 'd-'))
                
        elif t_lower_stripped in ['l', '(l)']:
            i += 1
            if i >= len(tokens):
                parts.append('l')
                break
            is_last = (i == len(tokens) - 1)
            parts.append(process_amino_acid_token(tokens[i], is_last, is_cyclic, 'l-'))

        elif t_lower_stripped in ['dl', '(dl)', '(d/l)', 'd/l', 'd,l', '(d,l)']: 
            i += 1
            if i >= len(tokens):
                parts.append('d/l')
                break
            is_last = (i == len(tokens) - 1)
            parts.append(process_amino_acid_token(tokens[i], is_last, is_cyclic, 'dl-'))
                
        else:
            # Regular token (not preceded by d or l)
            if t_lower_stripped.startswith(('(d)')):
                base_prefix = 'd-'
                t = t[3:]
            elif t_lower_stripped.startswith(('(l)')):
                base_prefix = 'l-'
                t = t[3:]
            elif t_lower_stripped.startswith(('(d/l)')):
                base_prefix = '(d/l)-'
                t = t[5:]
            else:
               base_prefix = 'l-' 
            parts.append(process_amino_acid_token(t, is_last, is_cyclic, base_prefix)) 
        i+=1

    name = ''
    if is_cyclic:
        name += 'cyclo('
    
    if prefix:
        name += prefix + '-'
        
    name += '-'.join(parts)
    
    if suffix:
        name += ' ' + suffix
        
    if is_cyclic:
        name += ')'

    if counter_acid_suffix:
        name += ' ' + counter_acid_suffix
        
    return name


def looks_like_peptide_shorthand(potential_peptide: str) -> bool:
    """
    Check if a given string looks like a peptide shorthand by checking if it contains
    any of the standard amino acid abbreviations surrounded by hyphens.

    Args:
        potential_peptide (str): The string to check.

    Returns:
        bool: True if the string looks like a peptide shorthand, False otherwise.
    """
    for k,_ in AA_FULL.items():
        if f'-{k}-' in potential_peptide.lower():
            return True
    return False


def name_to_smiles_peptide(compound_name_list: List[str]) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Converts a list of peptide shorthand names to their corresponding SMILES strings by first 
    converting the shorthand to an IUPAC-like name, then resolving that with OPSIN.

    Args:
        compound_name_list (List[str]): A list of peptide shorthand names to be converted.

    Returns:
        Tuple[Dict[str, str], Dict[str, str]]: A tuple containing two dictionaries. The first dictionary maps each peptide shorthand name to its SMILES string, and the second dictionary maps each peptide shorthand name that failed conversion to its error message.
    """
    peptide_iupac_names = []
    peptide_iupac_to_shorthand_mapping = {}
    for compound_name in compound_name_list:
        if not looks_like_peptide_shorthand(compound_name):
            continue
        peptide_iupac = peptide_shorthand_to_iupac(compound_name)
        peptide_iupac_names.append(peptide_iupac)
        peptide_iupac_to_shorthand_mapping[peptide_iupac] = compound_name

    chemical_name_dict, failure_message_dict = name_to_smiles_opsin(peptide_iupac_names)
    chemical_name_dict = {peptide_iupac_to_shorthand_mapping[k]:v for k,v in chemical_name_dict.items()}
    failure_message_dict = {peptide_iupac_to_shorthand_mapping[k]:v for k,v in failure_message_dict.items()}

    return chemical_name_dict, failure_message_dict