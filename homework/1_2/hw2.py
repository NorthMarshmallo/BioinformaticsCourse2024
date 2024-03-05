from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices, PairwiseAligner


def align_Needleman_Wunsch(a, b, substitution_matrix, gap_score):

    n, m = len(a) + 1, len(b) + 1
    d = [[0] * m for _ in range(n)]
    for i in range(n):
        d[i][0] = i * gap_score
    for j in range(m):
        d[0][j] = j * gap_score
    for i in range(1, n):
        for j in range(1, m):
            d[i][j] = max(d[i - 1][j] + gap_score, d[i][j - 1] + gap_score,
                              d[i - 1][j - 1] + substitution_matrix[a[i - 1]][b[j - 1]])

    i, j = n - 1, m - 1
    seq_a = []
    seq_b = []
    while i > 0 or j > 0:
        score = d[i][j]
        score_diag = d[i - 1][j - 1]
        score_left = d[i - 1][j]
        if score == score_diag + substitution_matrix[a[i - 1]][b[j - 1]]:
            seq_a.append(a[i - 1])
            seq_b.append(b[j - 1])
            i -= 1
            j -= 1
        elif score == score_left + gap_score:
            seq_a.append(a[i - 1])
            seq_b.append('-')
            i -= 1
        else:
            seq_a.append('-')
            seq_b.append(b[j - 1])
            j -= 1
    seq_a.reverse()
    seq_b.reverse()
    return (d[n - 1][m - 1], ''.join(seq_a), ''.join(seq_b))

def align_affinity(a, b, substitution_matrix, open_gap_score, extend_gap_score):

    n, m = len(a) + 1, len(b) + 1
    d = [[0] * m for _ in range(n)]
    i_x = [[0] * m for _ in range(n)]
    i_y = [[0] * m for _ in range(n)]
    for i in range(1, n):
        d[i][0] = open_gap_score + (i - 1) * extend_gap_score
        i_y[i][0] = -10e6

    for j in range(1, m):
        d[0][j] = open_gap_score + (j - 1) * extend_gap_score
        i_x[0][j] = -10e6

    for i in range(1, n):
        for j in range(1, m):
            i_x[i][j] = max(d[i - 1][j] + open_gap_score, i_x[i - 1][j] + extend_gap_score)
            i_y[i][j] = max(d[i][j - 1] + open_gap_score, i_y[i][j - 1] + extend_gap_score)
            d[i][j] = max(i_x[i][j], i_y[i][j], d[i - 1][j - 1]
                          + substitution_matrix[a[i - 1]][b[j - 1]])

    i, j = n - 1, m - 1
    seq_a = []
    seq_b = []

    while i > 0 or j > 0:
        score = d[i][j]
        score_diag = d[i - 1][j - 1]
        if score == score_diag + substitution_matrix[a[i - 1]][b[j - 1]]:
            seq_a.append(a[i - 1])
            seq_b.append(b[j - 1])
            i -= 1
            j -= 1
        elif score == i_x[i][j] or j == 0:
            while i_x[i][j] != d[i - 1][j] + open_gap_score and i > 1:
                seq_a.append(a[i - 1])
                seq_b.append('-')
                i -= 1
            seq_a.append(a[i - 1])
            seq_b.append('-')
            i -= 1
        elif score == i_y[i][j] or i == 0:
            while i_y[i][j] != d[i][j - 1] + open_gap_score and j > 1:
                seq_a.append('-')
                seq_b.append(b[j - 1])
                j -= 1
            seq_a.append('-')
            seq_b.append(b[j - 1])
            j -= 1
    seq_a.reverse()
    seq_b.reverse()
    return (d[n - 1][m - 1], ''.join(seq_a), ''.join(seq_b))

def __test_task_1(fasta_sequences_1, fasta_sequences_2, substitution_matrix, gap_score):
    print('__test_task_1')
    print('True aligment:')
    algn = PairwiseAligner(substitution_matrix=substitution_matrix, gap_score=gap_score)
    for aligment in algn.align(fasta_sequence_1, fasta_sequence_2):
        print('Result:', aligment.score)
        print(aligment)

    print('Needleman-Wunsch aligment:')
    score, seq_1, seq_2 = align_Needleman_Wunsch(fasta_sequences_1,
                                                 fasta_sequences_2, substitution_matrix, gap_score)
    print('Result:', score)
    print('Sequence_1:', seq_1)
    print('Sequence_2:', seq_2)

def __test_task_2(fasta_sequences_1, fasta_sequences_2, substitution_matrix, open_gap_score, extend_gap_score):
    print('__test_task_2')
    print('True aligment:')
    algn = PairwiseAligner(substitution_matrix=substitution_matrix, open_gap_score=open_gap_score,
                           extend_gap_score=extend_gap_score)
    for aligment in algn.align(fasta_sequence_1, fasta_sequence_2):
        print('Result:', aligment.score)
        print(aligment)

    print('Affinity aligment:')
    score, seq_1, seq_2 = align_affinity(fasta_sequences_1, fasta_sequences_2,
                                         substitution_matrix, open_gap_score, extend_gap_score)

    print('Result:', score)
    print('Sequence_1:', seq_1)
    print('Sequence_2:', seq_2)


if __name__ == "__main__":
    fasta_sequence_1 = 'KEEVLA'
    fasta_sequence_2 = 'EVL'
    substitution_matrix = substitution_matrices.load('BLOSUM62')
    gap_score = -3
    open_gap_score = -4
    extend_gap_score = -3

    __test_task_1(fasta_sequence_1, fasta_sequence_2, substitution_matrix, gap_score)
    __test_task_2(fasta_sequence_1, fasta_sequence_2, substitution_matrix, open_gap_score, extend_gap_score)

    #print(substitution_matrix)

"""
    Теоретическая часть.
    1. Получить выравнивание строк n, m можно или делецией относительно первой строки из строк n - 1, m,
    или инсерцией из строк n, m - 1. Мутации для n, m уже лежат в данных блоках, поэтому рекуррентная
    формула:
    N(n, m) = N(n - 1, m) + N(n, m - 1) 
    
    2. Точная формула для выравниваний - количество сочетаний N(n, m) = C(n, n + m)
    Начальные условия:
    N(1, 0) = 1 = C(0, 1), N(0, 1) = 1 = C(0, 1)
    N(1, 1) = N(1, 0) + N(0, 1) = 2 = C(1, 1 + 1) - база индукции
    В предположении, что N(n - 1, m) = C(n - 1, n + m - 1), N(n, m - 1) = C(n, n + m - 1) получим
    N(n, m) = N(n - 1, m) + N(n, m - 1) = C(n - 1, n + m - 1) + C(n, n + m - 1) = C(n, n + m)
    Т.к.
    (n + m - 1)! / ((n - 1)! (m - 1)!) * (1 / n + 1 / m) = (n + m - 1)! ((n + m) / (n - 1)! (m - 1)! nm) = 
    (n + m)! / (n! m!)
    
    3. Приближение Стирлинга для факториала:
    n! ~ sqrt(2 pi n) (n / e)^n
    Отсюда:
    N(n, m) = (n + m)! / (n! m!) = sqrt(2 pi (n + m)) * (n + m / e)^(n + m) / 
    (sqrt(2 pi n) * (n / e)^n * sqrt(2 pi m) * (m / e)^m)
    
    N(n, m) = sqrt((n + m) / (n * m * 2pi)) * (n + m)^(n + m) / (n^n * m^m)
    """