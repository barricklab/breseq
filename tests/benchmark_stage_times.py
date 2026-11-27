#!/usr/bin/env python3
"""Time breseq's pipeline stages on fixed datasets, for before/after comparisons.

This is a measuring tool, not a test: nothing discovers or runs it. Run it from the worktree root,
inside the conda env, ON AN OTHERWISE IDLE MACHINE -- per-stage wall times from a shared machine
are not comparable, and that is the whole point of the numbers.

    conda run --no-capture-output -p ../../../env python3 tests/benchmark_stage_times.py \\
        --label baseline --repeats 3 -j 8 --out /path/to/benchmark_dir

Each dataset is run `--repeats` times into `<out>/<label>/<dataset>/run<N>` (the breseq output is
deleted afterwards unless --keep), its stage times are read from output/summary.html, and a table of
the median per stage is printed and written to `<out>/<label>.tsv`. Compare two labels with

    python3 tests/benchmark_stage_times.py --compare <out>/baseline.tsv <out>/after.tsv

The default datasets are three of the long tests' inputs (fetch them with `make test-data`): a 36-bp
single-end library, a 2x101 library and a 2x150 library, all against REL606. Add others with
--dataset NAME='<breseq arguments after -j and -o>'.
"""
import argparse
import html
import os
import re
import shutil
import statistics
import subprocess
import sys
import time

DOWNLOADS = 'tests/data/downloads'
REF = '-r %s/ltee_REL606/REL606.gbk' % DOWNLOADS
DEFAULT_DATASETS = [
    ('ltee_se36', '%s %s/ena_SRR098038/SRR098038.fastq.gz' % (REF, DOWNLOADS)),
    ('ltee_pe101', '%s %s/ena_SRR2584534/SRR2584534_1.fastq.gz %s/ena_SRR2584534/SRR2584534_2.fastq.gz' % (REF, DOWNLOADS, DOWNLOADS)),
    ('ltee_pe150', '%s %s/ena_SRR2588848/SRR2588848_1.fastq.gz %s/ena_SRR2588848/SRR2588848_2.fastq.gz' % (REF, DOWNLOADS, DOWNLOADS)),
]


def seconds(text):
    """'1 minute 58 seconds' / '41 seconds' / '2 hours 3 minutes 1 second' -> seconds."""
    total = 0
    for value, unit in re.findall(r'(\d+)\s+(hour|minute|second)', text):
        total += int(value) * {'hour': 3600, 'minute': 60, 'second': 1}[unit]
    return total


def stage_times(summary_html):
    """[(stage, seconds), ...] in pipeline order, from the Execution Times table."""
    text = open(summary_html).read()
    start = text.find('Execution Times')
    if start < 0:
        raise RuntimeError('no Execution Times table in %s' % summary_html)
    end = text.find('</table>', start)
    rows = re.findall(r'<tr>(.*?)</tr>', text[start:end], re.S)
    out = []
    for row in rows:
        cells = [html.unescape(re.sub(r'<[^>]+>', '', c)).strip()
                 for c in re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', row, re.S)]
        if len(cells) >= 4 and cells[0].lower() != 'step':
            out.append((cells[0], seconds(cells[-1])))
    return out


def run_label(args):
    datasets = list(DEFAULT_DATASETS) if not args.no_default_datasets else []
    for spec in args.dataset or []:
        name, rest = spec.split('=', 1)
        datasets.append((name, rest))
    results = {}   # dataset -> stage -> [seconds per repeat]
    order = {}
    wall = {}
    for name, dataset_args in datasets:
        for repeat in range(1, args.repeats + 1):
            out_dir = os.path.join(args.out, args.label, name, 'run%d' % repeat)
            if os.path.exists(out_dir):
                shutil.rmtree(out_dir)
            os.makedirs(os.path.dirname(out_dir), exist_ok=True)
            cmd = '%s -j %d -o %s %s' % (args.breseq, args.j, out_dir, dataset_args)
            sys.stderr.write('[%s %s run %d] %s\n' % (args.label, name, repeat, cmd))
            t0 = time.time()
            with open(out_dir + '.log', 'w') as log:
                status = subprocess.call(cmd, shell=True, stdout=log, stderr=subprocess.STDOUT)
            elapsed = time.time() - t0
            if status != 0:
                raise RuntimeError('breseq failed (%d): see %s.log' % (status, out_dir))
            times = stage_times(os.path.join(out_dir, 'output', 'summary.html'))
            order[name] = [stage for stage, _ in times]
            for stage, secs in times:
                results.setdefault(name, {}).setdefault(stage, []).append(secs)
            wall.setdefault(name, []).append(elapsed)
            if args.keep:
                # the GD is what a before/after comparison of RESULTS needs
                pass
            else:
                gd = os.path.join(out_dir, 'output', 'output.gd')
                if os.path.exists(gd):
                    shutil.copy(gd, out_dir + '.gd')
                shutil.rmtree(out_dir)

    table = os.path.join(args.out, args.label + '.tsv')
    with open(table, 'w') as f:
        f.write('dataset\tstage\tmedian_seconds\tmin_seconds\tmax_seconds\trepeats\n')
        for name, _ in datasets:
            for stage in order[name]:
                v = results[name][stage]
                f.write('%s\t%s\t%g\t%d\t%d\t%d\n' % (name, stage, statistics.median(v), min(v), max(v), len(v)))
            v = wall[name]
            f.write('%s\t%s\t%.0f\t%.0f\t%.0f\t%d\n' % (name, 'TOTAL (wall clock)', statistics.median(v), min(v), max(v), len(v)))
    sys.stderr.write('wrote %s\n' % table)
    print(open(table).read())


def read_table(path):
    rows = {}
    order = []
    for line in open(path).read().splitlines()[1:]:
        dataset, stage, median = line.split('\t')[:3]
        rows[(dataset, stage)] = float(median)
        order.append((dataset, stage))
    return rows, order


def compare(before_path, after_path):
    before, order = read_table(before_path)
    after, _ = read_table(after_path)
    print('%-12s %-60s %8s %8s %8s' % ('dataset', 'stage', 'before', 'after', 'change'))
    for key in order:
        if key not in after:
            continue
        b, a = before[key], after[key]
        change = ('%+.0f%%' % (100.0 * (a - b) / b)) if b > 0 else ''
        print('%-12s %-60s %8g %8g %8s' % (key[0], key[1][:60], b, a, change))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--label', help='name for this set of runs (e.g. baseline)')
    p.add_argument('--out', help='directory for runs and tables (not inside the repository)')
    p.add_argument('--breseq', default='src/breseq/breseq')
    p.add_argument('-j', type=int, default=8)
    p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--dataset', action='append', help="NAME='breseq arguments'; may be repeated")
    p.add_argument('--no-default-datasets', action='store_true')
    p.add_argument('--keep', action='store_true', help='keep each breseq output directory')
    p.add_argument('--compare', nargs=2, metavar=('BEFORE.tsv', 'AFTER.tsv'))
    args = p.parse_args()
    if args.compare:
        compare(*args.compare)
        return
    if not args.label or not args.out:
        p.error('--label and --out are required to run')
    run_label(args)


if __name__ == '__main__':
    main()
