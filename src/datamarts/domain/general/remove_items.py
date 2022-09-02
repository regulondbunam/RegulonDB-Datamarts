# Recursive function for remove empty items and arrays, null fields, etc
def remove_empty_items(d):
    if isinstance(d, dict):
        return dict((k, remove_empty_items(v)) for k, v in d.items()
                    if v or v == 0 or v is False and remove_empty_items(v) is not None)
    elif isinstance(d, list):
        return [remove_empty_items(v) for v in d
                if v or v == 0 or v is False and remove_empty_items(v) is not None]
    else:
        if d or d == 0 or d is False:
            return d
