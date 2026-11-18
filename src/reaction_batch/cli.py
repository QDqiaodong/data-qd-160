DESCRIPTION='Exact formula balancing and limiting-reagent batch accounting'
LOTS_DESCRIPTION='Ordered lot material accounting with lineage isolation'
import argparse, json, sys
from pathlib import Path
from .common import read, save
from .lot_model import analyze as analyze_lots
from .model import analyze


def run(source):
    data,digest=read(source)
    return {'schema_version':1,'source_sha256':digest,'result':analyze(data)}


def run_lots(source):
    data,digest=read(source)
    return {'schema_version':1,'source_sha256':digest,'result':analyze_lots(data)}


def _main(description, compute):
    parser=argparse.ArgumentParser(description=description)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    try:
        if args.source.resolve()==args.output.resolve():raise ValueError('input and output must differ')
        report=compute(args.source);save(args.output,report)
    except (ValueError,TypeError,KeyError,OSError,OverflowError,RecursionError) as error:
        print(str(error),file=sys.stderr);return 2
    print(json.dumps(report,sort_keys=True,ensure_ascii=False,allow_nan=False));return 0


def main():
    return _main(DESCRIPTION, run)


def lots_main():
    return _main(LOTS_DESCRIPTION, run_lots)


if __name__=='__main__':
    # ``python -m reaction_batch.cli`` remains the original reaction-batch
    # entry point; the reaction-lots console script uses lots_main.
    raise SystemExit(main())
