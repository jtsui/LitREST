import itertools
import progressbar
import pprint

pr = pprint.PrettyPrinter(indent=2)


def serialize_rxn(substrates, products):
    '''
    Serializes a list of substrates and products (either chemical ids or names) 
    into a readable and hashable string. 
    input: substrates = [28248, 2003], products = [10153, 10559]
    output: '2003 + 28248 => 10153 + 10559'
    '''
    return '%s => %s' % (' + '.join([str(x) for x in sorted(substrates)]),
                         ' + '.join([str(x) for x in sorted(products)]))


def deserialize_rxn(rxn):
    '''
    Deserializes reaction string into a tuple of substrate list, product list
    input: '2003 + 28248 => 10153 + 10559'
    output: substrates = [28248, 2003], products = [10153, 10559]
    '''
    subs = [y for y in rxn.split(' => ')[0].split(' + ')]
    prods = [y for y in rxn.split(' => ')[1].split(' + ')]
    return subs, prods


def pbar(size):
    """Returns a new progress bar"""
    bar = progressbar.ProgressBar(maxval=size,
                                  widgets=[progressbar.Bar('=', '[', ']'),
                                           ' ', progressbar.Percentage(),
                                           ' ', progressbar.ETA(),
                                           ' ', progressbar.Counter(),
                                           '/%s' % size])
    return bar


def lines_in_file(file_name):
    num_lines = 0
    with open(file_name, 'rb') as f_in:
        for line in f_in:
            num_lines += 1
    return num_lines


def flatten_list(alist):
    out_list = []
    for l in alist:
        if isinstance(l, list):
            out_list.extend(flatten_list(l))
        else:
            out_list.append(l)
    return out_list


def grouper(iterable, i=2):
    '''
    takes s and returns (s0,s1), (s1,s2), (s2, s3), ...
    '''
    # a list of i iterators
    iters = itertools.tee(iterable, i)
    # loop through the iterators
    for j in range(i):
        # advance next j times
        for k in range(j):
            next(iters[j], None)
    # * unpacks the list turning each element to an arg
    return itertools.izip(*iters)


def lget(alist, index, default):
    try:
        return alist[index]
    except IndexError:
        return default


def all_combinations(alist, min_length=1):
    combs = []
    for i in xrange(min_length, len(alist) + 1):
        els = [list(x) for x in itertools.combinations(alist, i)]
        combs.extend(els)
    return combs


def chunks(l, n):
    for i in xrange(0, len(l), n):
        yield l[i:i + n]


def count_lines(file_name):
    line_count = 0
    with open(file_name, 'rb') as f_in:
        for line in f_in:
            line_count += 1
    return line_count


def main():
    """Put test code here"""
    pass

if __name__ == '__main__':
    main()
