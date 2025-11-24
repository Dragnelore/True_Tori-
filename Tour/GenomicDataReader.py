from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List

from reader import Reader


class GenomicDataReader(Reader, ABC):
    """
    Абстрактный ридер геномных данных (sam/vcf).
    Наследуется от Reader.
    """

    @abstractmethod
    def get_chromosomes(self) -> List[str]:
        """
        Получить список хромосом, представленных в файле.
        """
        raise NotImplementedError

    @abstractmethod
    def get_reference_genome(self) -> str:
        """
        Получить название референсного генома из заголовка файла.
        """
        raise NotImplementedError

    @abstractmethod
    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """
        Проверить корректность координаты (хромосома + позиция).
        Конкретные правила задаются в наследниках.
        """
        raise NotImplementedError