from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

from reader import Reader


class SequenceReader(Reader, ABC):
    """
    Абстрактный ридер последовательностей (fasta/fastq).
    Наследуется от Reader.
    """

    @abstractmethod
    def get_sequence(self, seq_id: str) -> Any:
        """
        Получить последовательность по идентификатору.
        Конкретный тип возвращаемого объекта определяется в наследниках.
        """
        raise NotImplementedError

    @abstractmethod
    def validate_sequence(self, sequence: str) -> bool:
        """
        Проверка корректности символов последовательности.
        Конкретные правила задаются в наследниках.
        """
        raise NotImplementedError