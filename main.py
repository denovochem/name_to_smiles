from constants import AMINO_ACID_SUB_SITES, PROTECTING_GROUPS, SPECIAL_CASES, AA_FULL, N_CAPS, C_CAPS, COUNTER_ACIDS, GREEK_LETTERS, PREFIX_MAP

def generate_side_chain_protections():
    """Generate the complete side chain protections dictionary"""
    protections = {}
    
    # Generate standard combinations
    for aa_code, (aa_name, site) in AMINO_ACIDS.items():
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

def split_peptide_shorthand(shorthand):
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

def parse_protected_residue(token):
    """Parse a token that might contain side chain protection"""
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

def extract_prefix_modifier(token, aa_dict):
    """
    Extract prefix modifiers from amino acid tokens.
    
    Args:
        token: The amino acid token to process
        aa_dict: Dictionary of amino acids to check against
    
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
            if remaining_lower in aa_dict:
                return full_name, remaining, remaining_lower
    
    return None, token, token.lower().strip().replace(' ', '')

def process_amino_acid_token(next_token, is_last, is_cyclic, stereo_prefix='l-'):
    """
    Process an amino acid token and return the properly formatted name.
    """
    next_token_lower_stripped = next_token.lower().strip().replace(' ', '')
    
    # Handle Greek letter prefix
    greek_prefix = None
    if next_token and next_token[0] in GREEK_LETTERS:
        greek_prefix = next_token[0]
        next_token = next_token[1:]
        next_token_lower_stripped = next_token.lower().strip().replace(' ', '')

    # Handle prefix modifiers (methyl, acetyl, etc.)
    prefix_modifier, next_token, next_token_lower_stripped = extract_prefix_modifier(next_token, AA_FULL)

    is_last_updated = is_last and not is_cyclic
    
    # Try to parse the token
    base_aa, protection = parse_protected_residue(next_token_lower_stripped)
    if not base_aa and next_token_lower_stripped in AA_FULL:
        base_aa = AA_FULL[next_token_lower_stripped][0]
    
    if base_aa:
        # Known amino acid
        if is_last_updated:
            print(base_aa, is_last_updated)
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

def shorthand_to_iupac(shorthand: str) -> str:
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
