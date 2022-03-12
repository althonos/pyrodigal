import collections


class Record(collections.namedtuple("Record", ["id", "seq", "description"])):
    pass


def parse(path):

    if isinstance(path, str):
        file = open(path)
    else:
        file = path

    with file:

        id_ = None
        seq = []
        for line in file:
            l = line.strip()
            if line.startswith(">"):
                if id_ is not None:
                    yield Record(id_, "".join(seq), desc)
                id_ = line[1:].split()[0].strip()
                desc = " ".join(line[1:].split(maxsplit=1))
                seq = []
            elif l:
                seq.append(l)
        if id_ is not None:
            yield Record(id_, "".join(seq), desc)
        elif seq:
            raise ValueError("not in FASTA format")
