# FastqReader.py
from __future__ import annotations

from typing import Dict, Iterator, List, Optional

from SequenceReader import SequenceReader


class FastqReader(SequenceReader):
    """
    Ридер для формата FASTQ.

    Запись представляется словарём:
    {
        "id": str,
        "sequence": str,
        "quality": List[int]  # оценки качества по Фреду
    }
    """

    def _parse_line(self, line: str):
        """
        В этом классе базовый метод read() не используется,
        поэтому _parse_line реализован только для соответствия абстрактному
        интерфейсу. Вся логика чтения вынесена в переопределённый read().
        """
        return line

    def read(self) -> Iterator[Dict[str, object]]:
        """
        Генератор по записям FASTQ.

        На каждой итерации возвращает словарь с ключами "id",
        "sequence" и "quality".
        """
        with open(self._filename, "r", encoding="utf-8") as handle:
            while True:
                header = handle.readline()
                if not header:
                    break  # EOF
                seq = handle.readline()
                plus = handle.readline()
                qual = handle.readline()
                if not (seq and plus and qual):
                    break  # неполная запись в конце файла

                seq_id = header.strip()[1:] if header.startswith("@") else header.strip()
                sequence = seq.strip()
                quality_str = qual.strip()
                quality_scores = [ord(ch) - 33 for ch in quality_str]

                yield {
                    "id": seq_id,
                    "sequence": sequence,
                    "quality": quality_scores,
                }

    def get_quality_scores(self, seq_id: str) -> List[int]:
        """
        Получить список оценок качества для указанной последовательности.
        """
        for record in self.read():
            if record["id"] == seq_id:
                return record["quality"]  # type: ignore[return-value]
        return []

    def get_average_quality(self, seq_id: str) -> float:
        """
        Получить среднее качество по всем позициям для указанной последовательности.
        """
        scores = self.get_quality_scores(seq_id)
        if not scores:
            return 0.0
        return sum(scores) / len(scores)

    # Реализация абстрактных методов SequenceReader

    def get_sequence(self, seq_id: str) -> Optional[str]:
        """
        Получить последовательность по идентификатору из FASTQ-файла.
        """
        for record in self.read():
            if record["id"] == seq_id:
                return record["sequence"]  # type: ignore[return-value]
        return None

    def validate_sequence(self, sequence: str) -> bool:
        """
        Простая проверка: допустимы символы A, C, G, T, N (верхний регистр).
        """
        allowed = set("ACGTN")
        return all(base in allowed for base in sequence.upper())
