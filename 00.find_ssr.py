import argparse
import os.path
import pytrf
import pyfastx

def outf_name_fun(fa_file,file_path):
    base_name = os.path.basename(fa_file)
    if base_name.endswith('.fa'):
        base_name = base_name[:-3]
    if base_name.endswith('.fasta'):
        base_name = base_name[:-6]
    if base_name.endswith('.fna'):
        base_name = base_name[:-4]
    outf_name = base_name + '.retry.csv'
    out_file = os.path.join(file_path, outf_name)
    return out_file

def find_and_print_ssr(fa_file, out_name):
    with open(out_name, 'w') as outf:
        fa = pyfastx.Fastx(fa_file, uppercase=True)
        for name, seq in fa:
            for ssr in pytrf.STRFinder(name, seq, 15, 10, 8, 5, 5, 5):
                outf.write(ssr.as_string(',') + '\n')

def main():
    parser = argparse.ArgumentParser(description='Find and print SSRs in a FASTA file')
    parser.add_argument('-f', '--fa_file', required=True, help='input FASTA file')
    args = parser.parse_args()
    file_path = os.path.dirname(args.fa_file)

    out_file = outf_name_fun(args.fa_file,file_path)
    find_and_print_ssr(args.fa_file, out_file)

if __name__ == "__main__":
    main()

