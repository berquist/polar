from math import factorial

def unique(iterable):
    seen = set()
    for x in iterable:
        if x in seen:
            continue
        seen.add(x)
        yield x

def fac(order, term):
    return factorial(order) / (factorial(term.count(0)) * factorial(term.count(1)) * factorial(term.count(2)))
