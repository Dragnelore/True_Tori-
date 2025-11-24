# FastaReader.py
from __future__ import annotations

from typing import Dict, Iterator, Tuple, Optional

from SequenceReader import SequenceReader


class FastaReader(SequenceReader):
    """
    Ридер для формата FASTA.

    Запись представляется кортежем (seq_id, sequence),
    где seq_id – идентификатор последовательности без '>'.
    """

    def _parse_line(self, line: str):
        """
        В этом классе базовый метод read() не используется,
        поэтому _parse_line реализован только для соответствия абстрактному
        интерфейсу. Вся логика чтения вынесена в переопределённый read().
        """
        return line

    def read(self) -> Iterator[Tuple[str, str]]:
        """
        Генератор по записям FASTA.

        Возвращает кортеж (seq_id, sequence) для каждой последовательности.
        """
        seq_id: Optional[str] = None
        chunks: list[str] = []

        with open(self._filename, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if seq_id is not None:
                        sequence = "".join(chunks)
                        yield seq_id, sequence
                    seq_id = line[1:].strip()
                    chunks = []
                else:
                    chunks.append(line.strip())
            # последняя последовательность
            if seq_id is not None:
                sequence = "".join(chunks)
                yield seq_id, sequence

    def read_sequences(self) -> Dict[str, str]:
        """
        Прочитать все последовательности в словарь {seq_id: sequence}.
        """
        sequences: Dict[str, str] = {}
        for seq_id, sequence in self.read():
            sequences[seq_id] = sequence
        return sequences

    def get_sequence_length(self, seq_id: str) -> int:
        """
        Получить длину последовательности по идентификатору.
        """
        sequence = self.get_sequence(seq_id)
        return len(sequence) if sequence is not None else 0

    # Реализация абстрактных методов SequenceReader

    def get_sequence(self, seq_id: str) -> Optional[str]:
        """
        Получить последовательность по идентификатору из FASTA-файла.
        """
        for cur_id, sequence in self.read():
            if cur_id == seq_id:
                return sequence
        return None

    def validate_sequence(self, sequence: str) -> bool:
        """
        Простая проверка: допустимы символы A, C, G, T, N (верхний регистр).
        """
        allowed = set("ACGTN")
        return all(base in allowed for base in sequence.upper())
