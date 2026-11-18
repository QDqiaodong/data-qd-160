DESCRIPTION='Exact formula balancing and limiting-reagent batch accounting'
LOTS_DESCRIPTION='Ordered lot material accounting with lineage isolation'
ASSAY_DESCRIPTION='Exact composition bounds under element molar measurements'
PEAKS_DESCRIPTION='Shared-drift chromatographic peak identity reconciliation'
import argparse, json, sys
from pathlib import Path
from .common import read, save
from .assay import analyze as analyze_assay
from .lot_model import analyze as analyze_lots
from .peaks import analyze as analyze_peaks
from .model import analyze


def run(source):
    data,digest=read(source)
    return {'schema_version':1,'source_sha256':digest,'result':analyze(data)}


def run_lots(source):
    data,digest=read(source)
    return {'schema_version':1,'source_sha256':digest,'result':analyze_lots(data)}


def run_assay(source):
    data,digest=read(source)
    return {'schema_version':1,'source_sha256':digest,'result':analyze_assay(data)}


def run_peaks(source):
    data,digest=read(source)
    return {'schema_version':1,'source_sha256':digest,'result':analyze_peaks(data)}


def _main(description, compute):
    parser=argparse.ArgumentParser(description=description)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    try:
        if args.source.resolve()==args.output.resolve():raise ValueError('input and output must differ')
        report=compute(args.source);save(args.output,report)
    except (ValueError,TypeError,KeyError,OSError,OverflowError,RecursionError,RuntimeError) as error:
        print(str(error),file=sys.stderr);return 2
    print(json.dumps(report,sort_keys=True,ensure_ascii=False,allow_nan=False));return 0


def main():
    return _main(DESCRIPTION, run)


def lots_main():
    return _main(LOTS_DESCRIPTION, run_lots)


def assay_main():
    return _main(ASSAY_DESCRIPTION, run_assay)


def peaks_main():
    return _main(PEAKS_DESCRIPTION, run_peaks)


if __name__=='__main__':
    # ``python -m reaction_batch.cli`` remains the original reaction-batch
    # entry point; the reaction-lots console script uses lots_main.
    raise SystemExit(main())
