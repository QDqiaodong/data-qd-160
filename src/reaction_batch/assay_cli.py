"""Console entry point for ``reaction-assay``."""
from .cli import assay_main


def main():
    return assay_main()


if __name__ == '__main__':
    raise SystemExit(main())
