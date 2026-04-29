import importlib
import importlib.util

image_slicer = (
    importlib.import_module("image_slicer")
    if importlib.util.find_spec("image_slicer") is not None
    else None
)


def main() -> None:
    if image_slicer is None:
        raise ImportError("image_slicer is required for png_split")
    image_slicer.slice("Marburger_Schloss_024.jpg", 3)


if __name__ == "__main__":
    main()
