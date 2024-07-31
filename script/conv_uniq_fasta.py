#usage: conv_uniq_fa.py <fa> <mark>


import sys
fa = sys.argv[1] if 1 < len(sys.argv) else None
mark = sys.argv[2] if 2 < len(sys.argv) else None


mark = "" if not mark else "mark:%s;" % mark
def do(f):
    n = 0
    while True:
        head = f.readline().strip("> \t\n")
        seq  = f.readline().strip("> \t\n")
        if (not head) or (not seq):
            break
        sys.stdout.write(">%sindex:%d;%s\n%s\n" % (mark, n, head, seq))
        n += 1


if fa:
    with open(fa) as f:
        do(f)
else:
    do(sys.stdin)
