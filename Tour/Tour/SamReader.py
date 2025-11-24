from __future__ import annotations

from typing import Any, Dict, Iterator, List, Optional

from GenomicDataReader import GenomicDataReader


class SamReader(GenomicDataReader):
    """
    Ридер для формата SAM.

    Выравнивание представляется словарём с основными полями SAM.
    """

    # ---------- Реализация абстрактного интерфейса Reader ----------

    def _parse_line(self, line: str) -> Optional[Dict[str, Any]]:
        """
        Преобразовать строку SAM (не заголовок) в словарь-выравнивание.
        Заголовочные строки ('@...') пропускаются.
        """
        if not line or line.startswith("@"):
            return None

        fields = line.split("\t")
        if len(fields) < 11:
            return None

        alignment: Dict[str, Any] = {
            "QNAME": fields[0],
            "FLAG": int(fields[1]),
            "RNAME": fields[2],
            "POS": int(fields[3]),
            "MAPQ": int(fields[4]),
            "CIGAR": fields[5],
            "RNEXT": fields[6],
            "PNEXT": int(fields[7]) if fields[7].isdigit() else 0,
            "TLEN": int(fields[8]) if fields[8].lstrip("-").isdigit() else 0,
            "SEQ": fields[9],
            "QUAL": fields[10],
            "TAGS": fields[11:],
        }
        return alignment

    # ---------- Методы, указанные в UML для SamReader ----------

    def read_alignments(self) -> List[Dict[str, Any]]:
        """
        Прочитать все выравнивания из SAM-файла.
        """
        return list(self.read())

    def get_header(self) -> Dict[str, List[str]]:
        """
        Получить заголовок SAM в виде словаря:
        {
            "HD": [...],
            "SQ": [...],
            "RG": [...],
            "PG": [...],
            "CO": [...],
        }
        """
        header: Dict[str, List[str]] = {}
        with open(self._filename, "r", encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith("@"):
                    break
                line = line.rstrip("\n")
                tag = line[1:3]
                header.setdefault(tag, []).append(line)
        return header

    def filter_alignments(self, flag: int) -> List[Dict[str, Any]]:
        """
        Отфильтровать выравнивания по значению поля FLAG.
        """
        return [aln for aln in self.read() if aln["FLAG"] == flag]

    def calculate_coverage(self, chrom: str) -> Dict[int, int]:
        """
        Простейшая оценка покрытия: число начал выравниваний
        в каждой позиции (игнорируя CIGAR).
        """
        coverage: Dict[int, int] = {}
        for aln in self.read():
            if aln["RNAME"] != chrom:
                continue
            pos = aln["POS"]
            coverage[pos] = coverage.get(pos, 0) + 1
        return coverage

    # ---------- Реализация абстрактного интерфейса GenomicDataReader ----------

    def get_chromosomes(self) -> List[str]:
        chroms: set[str] = set()
        for aln in self.read():
            rname = aln["RNAME"]
            if rname != "*":
                chroms.add(rname)
        return sorted(chroms)

    def get_reference_genome(self) -> str:
        """
        Пытается извлечь название референса из строк @SQ или @RG (поле AS:).
        Если не найдено, возвращает пустую строку.
        """
        header = self.get_header()
        lines = header.get("SQ", []) + header.get("RG", [])
        for line in lines:
            parts = line.split("\t")[1:]
            for field in parts:
                if field.startswith("AS:"):
                    return field[3:]
        return ""

    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """
        Проверяет, что хромосома присутствует в файле и позиция положительна.
        """
        return chrom in self.get_chromosomes() and pos > 0
