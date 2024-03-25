import argparse
from dataclasses import dataclass
import glob
import os.path
from typing import List


@dataclass
class Stats:
    sample_name: str
    total: int
    primary_cnt: int
    primary_per: str


def get_file_path(path: str) -> List[str]:
    assert os.path.isdir(path), f"{path} is not found."
    path_file: str = os.path.join(path, "*flagstat*")
    result_list: List[str] = [i for i in glob.glob(path_file) if os.path.isfile(i)]
    return result_list


def aggregate(file_path_list: List[str]) -> List[Stats]:
    result_list: List[Stats] = list()

    for fstat_file in file_path_list:
        sample_name = os.path.basename(fstat_file).replace(".flagstat", "")
        total = ""
        primary_per = ""
        primary_cnt = ""
        with open(fstat_file) as fin:
            for line in fin:
                line: str = line.strip()
                tmp_list: List[str] = line.split("\t")
                if "total (QC-passed reads + QC-failed reads)" == tmp_list[2]:
                    total = int(tmp_list[0])
                if 'primary mapped' == tmp_list[2]:
                    primary_cnt = int(tmp_list[0])
                if 'primary mapped %' == tmp_list[2]:
                    primary_per = tmp_list[0]
            tmp_stat: Stats = Stats(sample_name=sample_name, total=total, primary_cnt=primary_cnt, primary_per=primary_per)
            result_list.append(tmp_stat)
    return result_list


def write_output(stat_list: List[Stats], output_filename: str):
    with open(output_filename, 'w') as fout:
        outline: str = f"sample_name,total_n_reads,primary_n_reads,primary_percent\n"
        fout.write(outline)
        for item in stat_list:
            sample_name: str = item.sample_name
            total_n_reads: int = item.total
            primary_n_reads: int = item.primary_cnt
            primary_percent: str = item.primary_per
            line: str = f"{sample_name},{total_n_reads},{primary_n_reads},{primary_percent}\n"
            fout.write(line)


@dataclass
class Args:
    fstat_path: str
    output_file_name: str


def get_args() -> Args:
    parser = argparse.ArgumentParser(
        description="doh - summarize flagstat files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "-i",
        "--fstat_path",
        help="fstat_path",
        metavar="fstat_path",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-o",
        "--output_file_name",
        help="output_file_name",
        metavar="output_file_name",
        type=str,
        required=True,
    )

    args = parser.parse_args()
    return Args(args.fstat_path, args.output_file_name)


def main():
    args = get_args()
    fstat_path: str = args.fstat_path
    output_file_name: str = args.output_file_name

    fstat_path_list:List[str] = get_file_path(fstat_path)
    fstat_list:List[Stats] = aggregate(fstat_path_list)

    write_output(fstat_list, output_file_name)


if __name__ == '__main__':
    main()
