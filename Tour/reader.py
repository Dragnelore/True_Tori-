# reader.py
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Iterator, Any


class Reader(ABC):
    """
    Абстрактный базовый класс для всех ридеров.

    Атрибуты
    --------
    _filename: str
        Путь к файлу, с которым работает ридер.
    """

    def __init__(self, filename: str) -> None:
        self._filename = filename

    def read(self) -> Iterator[Any]:
        """
        Базовый генератор по строкам файла.

        Для каждой строки файла вызывает защищённый метод _parse_line,
        который должен быть реализован в наследниках.

        Возвращает итератор по объектам "Record" (тип задаётся наследниками).
        """
        with open(self._filename, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.rstrip("\n")
                record = self._parse_line(line)
                if record is not None:
                    yield record

    def close(self) -> None:
        """
        Метод оставлен для соответствия UML-диаграмме.
        В текущей реализации файловый дескриптор всегда
        закрывается контекстным менеджером в read().
        """
        # Ничего не делаем, т.к. файл открывается в контекстном менеджере.
        return None

    @abstractmethod
    def _parse_line(self, line: str) -> Any:
        """
        Преобразование строки файла в объект записи.
        Должен быть реализован в наследниках.
        """
        raise NotImplementedError