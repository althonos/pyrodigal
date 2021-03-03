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
        for line in file:
            if line.startswith(">"):
                if id_ is not None:
                    yield Record(id_, "".join(seq), desc)
                id_ = line[1:].split()[0].strip()
                desc = " ".join(line[1:].split(maxsplit=1))
                seq = []
            else:
                seq.append(line.strip())
        if id_ is not None:
            yield Record(id_, "".join(seq), desc)
