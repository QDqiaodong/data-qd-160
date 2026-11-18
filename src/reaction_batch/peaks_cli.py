"""Console entry point for ``reaction-peaks``."""
from .cli import peaks_main


def main():
    return peaks_main()


if __name__ == '__main__':
    raise SystemExit(main())
