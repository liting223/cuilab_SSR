#! /usr/bin
import argparse
import os.path
import pyfastx
from itertools import islice


def seq_dic_fun(fa_file):
    seq_dic = {}
    fa = pyfastx.Fastx(fa_file, uppercase=True)
    for name, seq in fa:
        seq_dic[name] = seq
    return seq_dic


def ssr_arrs_fun(ssr_file):
    ssr_arrs = []
    inf = open(ssr_file, 'r')
    for line in islice(inf, 0, None):
        chr = line.strip().split(',')[0]  # chr name  1H
        ssr_arrs.append(chr + '+' + line.strip().split(',')[1] + '+' + line.strip().split(',')[2])
    return ssr_arrs


def outf_name_fun(fa_file,file_path):
    base_name = os.path.basename(fa_file)
    if base_name.endswith('.fasta'):
        base_name = base_name[:-6]
    outf_name = base_name + '_seq.csv'
    out_file = os.path.join(file_path, outf_name)
    return out_file


def ssr_seq_fun(seq_dic, ssr_arrs, out_file):
    outf = open(out_file, 'w')
    flask_region = 500
    for ssr_arr in ssr_arrs:
        chr = ssr_arr.split('+')[0]
        ssr_start = int(ssr_arr.split('+')[1])
        ssr_end = int(ssr_arr.split('+')[2])
        if chr not in seq_dic:
            continue
        seq_len = len(seq_dic[chr])
        if ssr_start - flask_region > 0:
            start = ssr_start - flask_region - 1  # index -1
        else:
            start = 0
        if ssr_end + flask_region < seq_len:
            end = ssr_end + flask_region
        else:
            end = seq_len
        if chr in seq_dic:
            seq = seq_dic[chr][start:end]
            line = chr + "_" + str(ssr_start) + "_" + str(ssr_end) + "," + seq + "\n"
            outf.write(line)
    outf.close()


def main():
    parser = argparse.ArgumentParser(description="Extract flanking sequences based on SSR sites")
    parser.add_argument('-f', '--fa_file', required=True, help='input Fasta file')
    parser.add_argument('-s', '--ssr_file', required=True, help='input ssr file')
    args = parser.parse_args()
    file_path = os.path.dirname(args.fa_file)

    seq_dic = seq_dic_fun(args.fa_file)
    ssr_arrs = ssr_arrs_fun(args.ssr_file)
    out_file = outf_name_fun(args.fa_file,file_path)
    ssr_seq_fun(seq_dic, ssr_arrs, out_file)


if __name__ == "__main__":
    main()
