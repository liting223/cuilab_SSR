from itertools import islice
import os
import argparse

def outf_name_fun(ssr_file,file_path):
    base_name = os.path.basename(ssr_file)
    if base_name.endswith('.csv'):
        base_name = base_name[:-4]
    outf_name = base_name + '_ssr.gff3'  # 1.horduem.csv  --> 1.hordeum_genic.csv
    out_file = os.path.join(file_path, outf_name)
    return out_file


def read_ssr(ssr_file,out_file):
    i = 0
    inf = open(ssr_file,'r')
    outf = open(out_file,'w')
    source = "wheatssr"
    type = "ssr"
    for line in islice(inf,0,None):
        i += 1
        field = line.strip().split(',')
        chr = field[0]
        start = field[1]
        end = field[2]
        ssr = field[3]
        ssr_type = field[4]
        repeate = field[5]
        length = field[6]

        new_line = chr + "\t" + source + "\t" + type + "\t" + str(start) + "\t" + str(end) + "\t" + "." + "\t" + "+" + "\t" +"." +"\t" + "ID=ssr" + str(i) + ";" + "Motif=" + ssr + ";" + "Type=" + str(ssr_type) + ";" + "Repeate=" + str(repeate) + ";" + "Length=" + str(length) + "\n"
        outf.write(new_line)

def main():
    parser = argparse.ArgumentParser()

    # 添加命令行参数 -f，并指定其类型为字符串
    parser.add_argument("-s", "--ssr_file", type=str, help="input ssr file")
    # 解析命令行参数
    args = parser.parse_args()
    file_path = os.path.dirname(args.ssr_file)

    out_file = outf_name_fun(args.ssr_file,file_path)
    read_ssr(args.ssr_file,out_file)

if __name__ == '__main__':
    main()


