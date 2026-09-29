"""Restore missing or damaged ANTILLES tables using only Python's standard library."""

import argparse
from pathlib import Path
import os
import tempfile
import zipfile
import zlib


def ensure_antilles(directory):
    archive = directory / "logPkRatio_antilles.zip"
    expected = [f"logPkRatio_antilles-{i}.dat" for i in range(1, 401)]
    restored = 0
    with zipfile.ZipFile(archive) as bundle:
        if len(bundle.namelist()) != 400 or set(bundle.namelist()) != set(expected):
            raise ValueError(f"{archive}: expected exactly ANTILLES tables 1..400")
        for name in expected:
            info = bundle.getinfo(name)
            target = directory / name
            if target.is_file() and target.stat().st_size == info.file_size:
                if zlib.crc32(target.read_bytes()) == info.CRC:
                    continue
            # ZIP reading validates the CRC before any existing table is replaced.
            data = bundle.read(name)
            fd, temporary = tempfile.mkstemp(prefix=".antilles-", dir=directory)
            try:
                with os.fdopen(fd, "wb") as output:
                    output.write(data)
                    os.fchmod(output.fileno(), 0o644)
                os.replace(temporary, target)
            finally:
                if os.path.exists(temporary):
                    os.unlink(temporary)
            restored += 1
    if restored:
        print(f"ANTILLES: extracted {restored} tables into {directory}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    try:
        ensure_antilles(args.directory)
    except (OSError, ValueError, zipfile.BadZipFile, RuntimeError) as error:
        parser.exit(1, f"ANTILLES extraction failed: {error}\n")
