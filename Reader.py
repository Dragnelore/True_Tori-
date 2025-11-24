from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Iterator, Any


class Reader(ABC):
    """
    Абстрактный базовый класс для чтения файлов с биологическими данными.

    Атрибуты
    -------- 
    _filename : str
        Имя (путь) к файлу, из которого читаются данные.
    """
    def __init__(self, filename: str) -> None:
        self._filename: str = filename
        self._file = open(self._filename, mode="r", encoding="utf-8")

    def read(self) -> Iterator[Any]:
        """
        Читает файл построчно и возвращает итератор по объектам Record.

        Для каждой строки вызывает защищённый абстрактный метод `_parse_line`,
        который должен быть реализован в подклассах.
        """
        for line in self._file:
            # Убираем символы перевода строки, чтобы парсер работал с "чистой" строкой
            line = line.rstrip("\n\r")
            record = self._parse_line(line)
            # Если парсер вернул None, просто пропускаем строку
            if record is not None:
                yield record

    def close(self) -> None:
        """Закрывает открытый файл."""
        if not self._file.closed:
            self._file.close()

    @abstractmethod
    def _parse_line(self, line: str) -> Any:
        """
        Преобразует одну строку файла в объект Record.
        Должен быть реализован в каждом конкретном наследнике ридера
        (FastaReader, FastqReader, SamReader, VcfReader и т.п.).
        """
        raise NotImplementedError
