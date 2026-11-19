"""Console entry point for ``reaction-certify``."""
from .cli import certify_main


def main():
    return certify_main()


if __name__ == '__main__':
    raise SystemExit(main())
