# example_usage.py

from FastaReader import FastaReader
from FastqReader import FastqReader
from SamReader import SamReader
from VcfReader import VcfReader

import pandas as pd
import matplotlib.pyplot as plt


DATA_DIR = "data"


def demo_fasta():
    print("=== FASTA ===")
    fasta_path = f"{DATA_DIR}/test.fasta"
    reader = FastaReader(fasta_path)

    lengths = []
    for seq_id, seq in reader.read():
        print(f"ID: {seq_id}, length: {len(seq)}")
        lengths.append(len(seq))

    if lengths:
        print(f"Количество последовательностей: {len(lengths)}")
        print(f"Средняя длина: {sum(lengths) / len(lengths):.2f}")
    print()


def demo_fastq():
    print("=== FASTQ ===")
    fastq_path = f"{DATA_DIR}/test.fastq"
    reader = FastqReader(fastq_path)

    lengths = []
    records = list(reader.read())
    for rec in records:
        seq_id = rec["id"]
        seq = rec["sequence"]
        qual = rec["quality"]
        print(f"ID: {seq_id}, length: {len(seq)}, mean qual: {sum(qual)/len(qual):.2f}")
        lengths.append(len(seq))

    if lengths:
        print(f"Количество последовательностей: {len(lengths)}")
        print(f"Средняя длина: {sum(lengths) / len(lengths):.2f}")

    # --- Примеры простых графиков качества ---

    # Per-base sequence quality: среднее качество по каждой позиции
    if records:
        max_len = max(len(r["quality"]) for r in records)
        per_base_q = []
        positions = list(range(1, max_len + 1))
        for pos in range(max_len):
            vals = []
            for r in records:
                qs = r["quality"]
                if pos < len(qs):
                    vals.append(qs[pos])
            if vals:
                per_base_q.append(sum(vals) / len(vals))
            else:
                per_base_q.append(0)

        plt.figure()
        plt.plot(positions, per_base_q)
        plt.xlabel("Позиция в прочтении")
        plt.ylabel("Среднее качество (Phred)")
        plt.title("Per base sequence quality (FASTQ)")
        plt.tight_layout()
        plt.show()

    # Распределение длины прочтений
    if lengths:
        plt.figure()
        plt.hist(lengths, bins=range(1, max(lengths) + 2))
        plt.xlabel("Длина прочтения")
        plt.ylabel("Количество")
        plt.title("Sequence length distribution (FASTQ)")
        plt.tight_layout()
        plt.show()

    print()


def demo_sam():
    print("=== SAM ===")
    sam_path = f"{DATA_DIR}/test.sam"
    reader = SamReader(sam_path)

    # Заголовок
    header = reader.get_header()
    print("Группы заголовков:")
    for tag, lines in header.items():
        print(f"{tag}: {len(lines)} строк")
    print()

    # Выравнивания (в виде DataFrame)
    alignments = list(reader.read())
    df = pd.DataFrame(alignments)
    print("Первые строки таблицы выравниваний:")
    print(df.head())
    print()

    # Количество выравниваний
    print(f"Количество выравниваний: {len(df)}")

    # Статистика: количество выравниваний по хромосоме
    counts_by_chrom = df["RNAME"].value_counts()
    print("\nСтатистика 'количество выравниваний - хромосома':")
    print(counts_by_chrom)
    print()

    # Пример "интерсекта": выравнивания на chr1 в диапазоне [10, 30]
    chrom = "chr1"
    start = 10
    end = 30
    interval_alignments = df[(df["RNAME"] == chrom) & (df["POS"] >= start) & (df["POS"] <= end)]
    print(f"Выравнивания на {chrom}:{start}-{end}:")
    print(interval_alignments)
    print()


def demo_vcf():
    print("=== VCF ===")
    vcf_path = f"{DATA_DIR}/test.vcf"
    reader = VcfReader(vcf_path)

    header = reader.get_header()
    print("Количество meta-строк:", len(header["meta"]))
    print("Колонки:", header["columns"])
    print()

    variants = list(reader.read())
    df = pd.DataFrame(variants)
    print("Первые строки таблицы вариантов:")
    print(df[["CHROM", "POS", "ID", "REF", "ALT", "QUAL"]].head())
    print()

    print(f"Количество вариантов: {len(df)}")

    # Статистика: количество вариантов по региону (для простоты по хромосоме)
    counts_by_chrom = df["CHROM"].value_counts()
    print("\nСтатистика 'количество вариантов - хромосома':")
    print(counts_by_chrom)
    print()

    # Пример выборки вариантов из геномного интервала chr1:10-30
    chrom = "chr1"
    start = 10
    end = 30
    interval_variants = df[(df["CHROM"] == chrom) & (df["POS"] >= start) & (df["POS"] <= end)]
    print(f"Варианты на {chrom}:{start}-{end}:")
    print(interval_variants[["CHROM", "POS", "ID", "REF", "ALT", "QUAL"]])
    print()

    # Пример: получить генотип SAMPLE1 для первого варианта
    if variants:
        first = variants[0]
        gt = reader.get_genotype("SAMPLE1", first)
        print("Генотип SAMPLE1 для первого варианта:", gt)
    print()


def main():
    demo_fasta()
    demo_fastq()
    demo_sam()
    demo_vcf()


if __name__ == "__main__":
    # Просто запускаем демонстрацию без какого-либо CLI.
    main()
