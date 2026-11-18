"""Console entry point for ``reaction-lots``."""
from .cli import lots_main


def main():
    return lots_main()


if __name__ == '__main__':
    raise SystemExit(main())
