"""Minimal FASTA I/O (no Biopython required)."""

from __future__ import annotations

from pathlib import Path
from typing import Iterator, List, Tuple


def parse_fasta(path: str | Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    pid: str | None = None
    chunks: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if pid is not None:
                    records.append((pid, "".join(chunks)))
                header = line[1:].strip()
                pid = header.split()[0]
                # Prefer UniProt accession when header is sp|ACC|NAME
                if "|" in pid:
                    parts = header.split("|")
                    if len(parts) >= 2 and parts[0] in {"sp", "tr"}:
                        pid = parts[1]
                chunks = []
            else:
                chunks.append(line)
    if pid is not None:
        records.append((pid, "".join(chunks)))
    return records


def iter_fasta(path: str | Path) -> Iterator[Tuple[str, str]]:
    yield from parse_fasta(path)


def write_fasta(path: str | Path, records: List[Tuple[str, str]], width: int = 60) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for pid, seq in records:
            fh.write(f">{pid}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")
