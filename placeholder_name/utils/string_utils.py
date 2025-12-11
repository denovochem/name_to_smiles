from placeholder_name.utils.constants import (
    COMMON_CHARS_WHITELIST,
    NON_LATIN1_REPLACEMENTS,
)

def safe_str(x):
    """Tries to convert to string, returns none upon exception"""
    try:
        return str(x)
    except Exception:
        return None

def filter_strings_by_whitelist(strings, whitelist=COMMON_CHARS_WHITELIST):
    """
    Filters a list of strings, returning only those that contain
    only characters from the whitelist.

    :param strings: List of strings to filter.
    :param whitelist: List or set of allowed characters.
    :return: List of strings that contain only whitelisted characters.
    """
    whitelist_set = set(whitelist)  # Convert to set for faster lookup
    return [s for s in strings if set(s).issubset(whitelist_set)]

def clean_strings(string, chars_to_replace_dict=NON_LATIN1_REPLACEMENTS):
    """
    Removes specified characters from a string.

    Args:
        string: The input string.
        chars_to_remove: A list (or string or set) of characters to remove.

    Returns:
        The cleaned string with the specified characters removed.
    """
    for char, replacement in chars_to_replace_dict.items():
        string = string.replace(char, replacement)
    string = string.replace('\n', '')
    return string

def is_latin1_compatible(s):
    """Check if a string is compatible with the latin-1 codec."""
    try:
        s.encode('latin-1')
        return True
    except UnicodeEncodeError:
        return False

def filter_latin1_compatible(strings):
    """Filter a list of strings to only include those compatible with the latin-1 codec."""
    return [s for s in strings if is_latin1_compatible(s)]