import sys

_last_collection_name = None  # variable estática global

def print_progress(current, total, collection_name, bar_length=40):
    global _last_collection_name

    fraction = current / total if total else 1
    filled = int(bar_length * fraction)
    bar = "█" * filled + "-" * (bar_length - filled)
    percent = int(fraction * 100)

    # Si la colección cambió, imprimimos en una línea nueva
    if _last_collection_name != collection_name:
        sys.stdout.write("\n")

    # Imprimimos en la misma línea solo si es la misma colección
    sys.stdout.write(
        f"\rProcessing {collection_name}: |{bar}| {percent}% ({current}/{total})"
    )
    sys.stdout.flush()

    _last_collection_name = collection_name