
import argparse
from itertools import islice
import os

# usage python demo.py -g demo.gff3 -s ssr.csv


def get_region_loc(line):
    fields = line.strip().split('\t')
    chr = fields[0]
    region = fields[2]
    start = fields[3]
    end = fields[4]
    region_loc = str(start) + ',' + str(end)
    return chr, region, region_loc

def process_gff(gff_file):
    inf = open(gff_file,'r')
    gene_list = [[] for _ in range(11)]
    cds_list = [[] for _ in range(11)]
    utr_list = [[] for _ in range(11)]
    chr_dic = {}

    for line in islice(inf, 0, None):
        if not line.startswith('#'):
            chr, region, region_loc = get_region_loc(line)
            start = int(region_loc.strip().split(',')[0])
            start_nr = start // 100000000

            if region == 'gene':
                gene_list[start_nr].append(region_loc)
            elif 'UTR' in region:
                utr_list[start_nr].append(region_loc)
            elif region == 'CDS':
                cds_list[start_nr].append(region_loc)

            chr_dic[chr] = [gene_list, cds_list, utr_list]

    return chr_dic

def process_ssr_genic(ssr_file, chr_dic):
    inf = open(ssr_file,'r')
    genic_list = []
    for line in islice(inf, 0, None):
        data = line.strip().split(',')
        ssr_chr = data[0]
        ssr = data[3]
        ssr_type = data[4]
        ssr_repeat = data[5]
        ssr_start = int(data[1])
        ssr_end = int(data[2])
        ssr_len = data[6]

        gene_list = chr_dic[ssr_chr][0]
        ssr_start_nr = ssr_start // 100000000
        for gene_loc in gene_list[ssr_start_nr]:
            start = int(gene_loc.strip().split(',')[0])
            end = int(gene_loc.strip().split(',')[1])
            # start, end = map(int, gene_loc.strip().split(','))
            if ssr_start >= start and ssr_end <= end:
                # genic_line = f"{ssr_chr},{ssr},{ssr_type},{ssr_repeat},{ssr_start},{ssr_end},{ssr_len},genic,{start},{end}\n"
                genic_line = f"{ssr_chr},{ssr},{ssr_type},{ssr_repeat},{ssr_start},{ssr_end},{ssr_len},intron\n"
                genic_list.append(genic_line)
                break
    return genic_list

def process_cds(genic_list,chr_dic):
    for i in islice(range(len(genic_list)), 0, None):
        data = genic_list[i].strip().split(',')
        ssr_chr = data[0]
        ssr = data[1]
        ssr_type = data[2]
        ssr_repeat = data[3]
        ssr_start = int(data[4])
        ssr_end = int(data[5])
        ssr_len = data[6]

        cds_list = chr_dic[ssr_chr][1]
        ssr_start_nr = ssr_start // 100000000
        for cds_loc in cds_list[ssr_start_nr]:
            start = int(cds_loc.strip().split(',')[0])
            end = int(cds_loc.strip().split(',')[1])
            # start, end = map(int, cds_loc.strip().split(','))
            if ssr_start >= start and ssr_end <= end:
                genic_line = f"{ssr_chr},{ssr},{ssr_type},{ssr_repeat},{ssr_start},{ssr_end},{ssr_len},CDS\n"
                genic_list[i] = genic_line
                break

    return  genic_list

def outf_name_fun(ssr_file,file_path):
    base_name = os.path.basename(ssr_file)
    if base_name.endswith('.csv'):
        base_name = base_name[:-4]
    outf_name = base_name + '_genic.csv'  # 1.horduem.csv  --> 1.hordeum_genic.csv
    out_file = os.path.join(file_path, outf_name)
    return out_file

def process_utr(genic_list,chr_dic,out_file):
    outf = open(out_file,'w')
    for i in islice(range(len(genic_list)), 0, None):
        data = genic_list[i].strip().split(',')
        ssr_chr = data[0]
        ssr = data[1]
        ssr_type = data[2]
        ssr_repeat = data[3]
        ssr_start = int(data[4])
        ssr_end = int(data[5])
        ssr_len = data[6]
        type = data[7]

        utr_list = chr_dic[ssr_chr][2]
        ssr_start_nr = ssr_start // 100000000

        for utr_loc in utr_list[ssr_start_nr]:
            start = int(utr_loc.strip().split(',')[0])
            end = int(utr_loc.strip().split(',')[1])
            # start, end = map(int, utr_loc.strip().split(','))
            if ssr_start >= start and ssr_end <= end:
                type = 'UTR'
                break
        genic_line = f"{ssr_chr},{ssr_start},{ssr_end},{ssr},{ssr_type},{ssr_repeat},{ssr_len},{type}\n"
        outf.write(genic_line)

def main():
    parser = argparse.ArgumentParser()

    # 添加命令行参数 -f，并指定其类型为字符串
    parser.add_argument("-g", "--gff_file", type=str, help="input gff3 file")
    parser.add_argument("-s", "--ssr_file", type=str, help="input ssr file")
    # 解析命令行参数
    args = parser.parse_args()
    file_path = os.path.dirname(args.ssr_file)

    chr_dic = process_gff(args.gff_file)
    genic_list = process_ssr_genic(args.ssr_file, chr_dic)
    genic_list = process_cds(genic_list, chr_dic)
    out_file = outf_name_fun(args.ssr_file,file_path)
    process_utr(genic_list, chr_dic, out_file)


if __name__ == '__main__':
    main()



