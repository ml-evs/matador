from pathlib import Path

with open(Path(__file__).parent.parent.joinpath("data/words").resolve(), "r") as f:
    WORDS = f.readlines()
with open(Path(__file__).parent.parent.joinpath("data/nouns").resolve(), "r") as f:
    NOUNS = f.readlines()

__all__ = ["NOUNS", "WORDS"]
