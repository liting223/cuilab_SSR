import argparse
import os.path
import pyfastx
from itertools import islice


def seq_generator_fun(fa_file):
    fa = pyfastx.Fastx(fa_file, uppercase=True)
    for name, seq in fa:
        yield name, seq


def ssr_arrs_generator_fun(ssr_file):
    with open(ssr_file, 'r') as inf:
        for line in inf:
            data_arr = line.strip().split(',')
            yield data_arr[0] + '+' + data_arr[1] + '+' + data_arr[2]


def outf_name_fun(ssr_file, file_path):
    base_name = os.path.basename(ssr_file)
    if base_name.endswith('.csv'):
        base_name = base_name[:-4]
    outf_name = base_name + '_seq.csv'
    out_file = os.path.join(file_path, outf_name)
    return out_file


def ssr_seq_fun(seq_generator, ssr_generator, out_file):
    with open(out_file, 'w') as outf:
        flask_region = 500
        for ssr_arr in ssr_generator:
            data_arr = ssr_arr.split('+')
            chr = data_arr[0]
            ssr_start = int(data_arr[1])
            ssr_end = int(data_arr[2])
            for name, seq in seq_generator:
                if name == chr:
                    seq_len = len(seq)
                    start = max(0, ssr_start - flask_region - 1)
                    end = min(seq_len, ssr_end + flask_region)
                    seq_flanking = seq[start:end]
                    line = chr + "_" + str(ssr_start) + "_" + str(ssr_end) + "," + seq_flanking + "\n"
                    outf.write(line)


def main():
    parser = argparse.ArgumentParser(description="Extract flanking sequences based on SSR sites")
    parser.add_argument('-f', '--fa_file', required=True, help='input Fasta file')
    parser.add_argument('-s', '--ssr_file', required=True, help='input ssr file')
    args = parser.parse_args()
    file_path = os.path.dirname(args.ssr_file)

    seq_generator = seq_generator_fun(args.fa_file)
    ssr_generator = ssr_arrs_generator_fun(args.ssr_file)
    out_file = outf_name_fun(args.ssr_file, file_path)
    ssr_seq_fun(seq_generator, ssr_generator, out_file)


if __name__ == "__main__":
    main()
