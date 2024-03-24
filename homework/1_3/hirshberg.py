from anytree import Node, RenderTree
import numpy as np

def levenshtein_distance(a, b, delet, ins, match, mismatch):
    n, m = len(a), len(b)

    d_prev, d_cur = [x * ins for x in range(m + 1)], [0] * (m + 1)

    for i in range(1, n + 1):
        d_cur[0] = d_prev[0] + delet
        for j in range(1, m + 1):
            if a[i - 1] != b[j - 1]:
                is_match = mismatch
            else:
                is_match = match
            d_cur[j] = max(d_prev[j] + delet, d_cur[j - 1] + ins,
                              d_prev[j - 1] + is_match)
        d_prev, d_cur = d_cur, d_prev

    return np.array(d_prev)

def build_tree_hirschberg(a, b, delet, ins, match, mismatch, node):
    n, m = len(a), len(b)
    if n <= 1 or m <= 1:
        return node
    mid = n // 2 - 1
    d_1 = levenshtein_distance(a[:mid + 1], b, delet, ins, match, mismatch)
    d_2 = levenshtein_distance(a[:mid:-1], b[::-1], delet, ins, match, mismatch)
    d = d_1 + d_2[::-1]
    i = np.argmax(d)
    node_l = Node(parent=node, name=(a[:mid + 1], b[:i]))
    node_r = Node(parent=node, name=(a[mid + 1:], b[i:]))
    build_tree_hirschberg(a[:mid + 1], b[:i], delet, ins, match, mismatch, node_l)
    build_tree_hirschberg(a[mid + 1:], b[i:], delet, ins, match, mismatch, node_r)
    return node

if __name__ == "__main__":
    a = 'AGTACGCA'
    b = 'TATGC'
    delet = -2
    ins = -2
    match = 2
    mismatch = -1

    node = Node(name=(a, b))
    tree = build_tree_hirschberg(a, b, delet, ins, match, mismatch, node)
    for pre, fill, node in RenderTree(tree):
        print("%s%s" % (pre, node.name))

    #print(substitution_matrix)

