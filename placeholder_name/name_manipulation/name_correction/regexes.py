PATTERNS = [
    # =================================================================
    # 1. CHARACTER SUBSTITUTION (l/I -> 1, O -> 0)
    # =================================================================
    # CASE: 'l' or 'I' acting as '1'
    # Look for: l/I NOT preceded by a letter.
    # Followed by: digit, separator, Space, Indicated Hydrogen (H), or Stereochem (R/S)
    # Example: "l, 2-", "(lS, 2R)", "lH-indole"
    (
        r"(?<![a-zA-Z])[lI](?=\s*(?:[,\-\d]|H\b|[RS]\b))",
        "1",
        "Locant: 'l/I' → '1'",
    ),
    # CASE: 'O' acting as '0'
    # Look for: O surrounded by digits or separators
    # Example: "2O,3-"
    (r"(?<=[\d,])O(?=[\d,\-])", "0", "Locant: 'O' → '0'"),
    (
        # Pattern: Matches 'Z' only when surrounded by hyphens,
        # and preceded specifically by lowercase letters or digits.
        # Matches: "pyrid-Z-yl", "onyloxy-Z-(", "3-Z-chloro"
        # Ignores: "(Z)-2-", "N-Z-glycine"
        r"(?<=[a-z0-9]-)Z(?=-)",
        "2",
        "Locant: 'Z' → '2' (e.g., 'pyrid-Z-yl')",
    ),
    # =================================================================
    # 2. PUNCTUATION CORRECTION (. -> ,)
    # =================================================================
    # CASE: Dots acting as Commas in Locants or Stereochem
    # Look for: dot preceded by Digit, 'R', or 'S'. Followed by Digit.
    # Example: "1.2-dichloro", "(1S.2R)"
    # Note: We use negative lookahead (?!.*]) strictly if you deal with Bicyclo[2.2.2]
    # but for general names, this pattern is usually safe.
    (
        r"(?<=[\dRS])\.(?=\s*\d)(?![^\[]*\])",
        ",",
        "Locant: '.' → ',' (ignoring [x.y.z] descriptors)",
    ),
    # =================================================================
    # 3. MISSING SEPARATORS (12-di -> 1,2-di)
    # =================================================================
    # CASE: Missing comma between digits before a Multiplier
    # Logic: If we see two digits followed by a hyphen and a multiplier (di, tri, bis),
    # it is almost certainly a list of positions, not carbon #12.
    # Example: "12-dichlorobenzene" -> "1,2-dichlorobenzene"
    (
        r"(?<!\d)(\d)(\d)\s*-(?=\s*(?:di|tri|tetra|penta|hexa|bis|tris|tetrakis))",
        r"\1,\2-",
        "Locant: Insert missing comma '12-di' → '1,2-di'",
    ),
    # =================================================================
    # 4. FORMATTING / WHITESPACE
    # =================================================================
    # CASE: Cleanup spaces inside locants
    # Example: "1, 2-di" -> "1,2-di"
    (r"(?<=,)\s+(?=[\dRS])", "", "Locant: Remove space after comma"),
]
