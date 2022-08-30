""" This submodule includes some useful utility functions
for exporting data.

"""

import os
import random


def file_writer_function(function):
    """Wrapper for file writers to safely overwrite/hash duplicate files.

    Keyword arguments:
        overwrite (bool): whether or not to overwrite colliding files.
        hash_dupe (bool): whether or not to create a unique filename for
            any colliding files, or just skip writing them.

    """

    from functools import wraps

    @wraps(function)
    def wrapped_writer(doc, path, overwrite=False, hash_dupe=False, **kwargs):
        """Wrap and return the writer function."""
        try:
            overwrite = kwargs.pop("overwrite", overwrite)
            hash_dupe = kwargs.pop("hash_dupe", hash_dupe)
            flines, ext = function(
                doc, path, overwrite=overwrite, hash_dupe=hash_dupe, **kwargs
            )
            if ext is not None and not path.endswith("." + ext):
                path += ".{}".format(ext)

            if os.path.isfile(path):
                if overwrite:
                    os.remove(path)
                elif hash_dupe:
                    path = "{}-{}.{}".format(
                        path.replace(ext, ""), generate_hash(), ext
                    )
                else:
                    return False

            with open(path, "w") as f:
                for line in flines:
                    f.write(line + "\n")

        except Exception as exc:
            raise type(exc)("Failed to write {}: {}".format(path, exc))

    return wrapped_writer


def generate_hash(hash_len=6):
    """Quick hash generator, based on implementation in PyAIRSS by J. Wynn.

    Keyword arguments:
        hash_len (int): desired length of hash.

    """
    hash_chars = [str(x) for x in range(10)] + [chr(97 + i) for i in range(25)]
    _hash = ""
    for _ in range(hash_len):
        _hash += random.choice(hash_chars)
    return _hash


def generate_relevant_path(**args):
    """Generates a suitable path name based on query."""
    dirname = ""
    if args.get("subcmd") is not None:
        dirname += args.get("subcmd") + "-"
    else:
        dirname = "query"
    if args.get("composition") is not None:
        for comp in args["composition"]:
            dirname += comp
    elif args.get("formula") is not None:
        dirname += args.get("formula")[0]
    if args.get("db") is not None:
        dirname += "-" + args.get("db")[0]
    if args.get("swap") is not None:
        for swap in args["swap"]:
            dirname += "-" + swap
        if args.get("hull_cutoff") is not None:
            dirname += "-hull-" + str(args.get("hull_cutoff")) + "eV"
    if args.get("id") is not None:
        dirname += "-" + args.get("id")[0] + "_" + args.get("id")[1]
    dirname = dirname.replace("--", "-")
    return dirname
