# VcfReader.py
from __future__ import annotations

from typing import Any, Dict, Iterator, List, Optional

from GenomicDataReader import GenomicDataReader


class VcfReader(GenomicDataReader):
    """
    Ридер для формата VCF.

    Вариант представляется словарём с основными полями VCF.
    """

    # ---------- Реализация абстрактного интерфейса Reader ----------

    def _parse_line(self, line: str) -> Optional[Dict[str, Any]]:
        """
        Преобразовать строку VCF с вариантом в словарь.
        Заголовочные строки ('##' и '#CHROM') пропускаются.
        """
        if not line or line.startswith("#"):
            return None

        fields = line.split("\t")
        if len(fields) < 8:
            return None

        chrom = fields[0]
        pos = int(fields[1])
        vid = fields[2]
        ref = fields[3]
        alt = fields[4].split(",") if fields[4] != "." else []
        qual_str = fields[5]
        qual: Optional[float]
        try:
            qual = float(qual_str)
        except ValueError:
            qual = None
        flt = fields[6]
        info_str = fields[7]

        info: Dict[str, Any] = {}
        if info_str and info_str != ".":
            for item in info_str.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    info[key] = value
                else:
                    info[item] = True

        variant: Dict[str, Any] = {
            "CHROM": chrom,
            "POS": pos,
            "ID": vid,
            "REF": ref,
            "ALT": alt,
            "QUAL": qual,
            "FILTER": flt,
            "INFO": info,
            "FORMAT": [],
            "SAMPLES": {},
        }

        if len(fields) > 8:
            format_fields = fields[8].split(":")
            variant["FORMAT"] = format_fields
            samples: Dict[str, Dict[str, Any]] = {}
            for idx, sample_field in enumerate(fields[9:], start=1):
                sample_name = f"SAMPLE{idx}"
                sample_values = sample_field.split(":")
                sample_dict: Dict[str, Any] = {}
                for key, value in zip(format_fields, sample_values):
                    sample_dict[key] = value
                samples[sample_name] = sample_dict
            variant["SAMPLES"] = samples

        return variant

    # ---------- Методы, указанные в UML для VcfReader ----------

    def read_variants(self) -> List[Dict[str, Any]]:
        """
        Прочитать все варианты из VCF-файла.
        """
        return list(self.read())

    def get_header(self) -> Dict[str, Any]:
        """
        Получить заголовок VCF в виде словаря:
        {
            "meta": [...],        # строки '##...'
            "columns": [ ... ],   # список названий колонок из строки '#CHROM'
        }
        """
        meta: List[str] = []
        columns: List[str] = []
        with open(self._filename, "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("##"):
                    meta.append(line.rstrip("\n"))
                elif line.startswith("#CHROM"):
                    columns = line.lstrip("#").rstrip("\n").split("\t")
                    break
        return {"meta": meta, "columns": columns}

    def filter_by_quality(self, min_qual: float) -> List[Dict[str, Any]]:
        """
        Отфильтровать варианты по минимальному значению QUAL.
        """
        result: List[Dict[str, Any]] = []
        for var in self.read():
            qual = var["QUAL"]
            if qual is not None and qual >= min_qual:
                result.append(var)
        return result

    def get_genotype(self, sample: str, variant: Dict[str, Any]) -> Optional[str]:
        """
        Получить генотип указанного образца для данного варианта.

        sample – имя образца (например, 'SAMPLE1').
        Возвращает строку GT или None, если информация отсутствует.
        """
        samples: Dict[str, Dict[str, Any]] = variant.get("SAMPLES", {})
        sample_info = samples.get(sample)
        if not sample_info:
            return None

        gt = sample_info.get("GT")
        return gt if isinstance(gt, str) else None

    # ---------- Реализация абстрактного интерфейса GenomicDataReader ----------

    def get_chromosomes(self) -> List[str]:
        chroms: set[str] = set()
        for var in self.read():
            chroms.add(var["CHROM"])
        return sorted(chroms)

    def get_reference_genome(self) -> str:
        """
        Пытается извлечь ссылку на референс из мета-строк '##reference='.
        Если не найдено, возвращает пустую строку.
        """
        header = self.get_header()
        for line in header["meta"]:
            if line.startswith("##reference="):
                return line.split("=", 1)[1]
        return ""

    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """
        Проверяет, что хромосома присутствует в файле и позиция положительна.
        """
        return chrom in self.get_chromosomes() and pos > 0
