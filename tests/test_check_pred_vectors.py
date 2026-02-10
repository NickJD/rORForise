import importlib.util
spec = importlib.util.spec_from_file_location('cp','/home/nick/Git/rORForise/src/rORForise/check_pred.py')
cp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cp)

# deterministic vectors: each is (cds_open, cds_close, cds_dir, read_open, read_close, read_dir, pred_start, pred_end, pred_dir, read_seq, expected_answer_names)
VECTORS = [
    # category ('+','+','+'): read covers cds start and stop; pred aligns exactly to start..stop
    (100, 300, '+', 90, 310, '+', 11, 211, '+', 'A'*211, ['correct direction','correct start','correct stop','correct frame']),
    # category ('+','-','-'): reverse read mapping; pred on '-' with appropriate coords gives correct start/stop
    (1000, 1200, '+', 990, 1210, '-', 11, 211, '-', 'A'*121, ['correct direction']),
    # category ('-','+','-'): cds on minus, read on plus, pred on -
    (2000, 2200, '-', 1990, 2210, '+', 11, 211, '-', 'A'*221, ['correct direction']),
    # category ('-','-','+')
    (3000, 3200, '-', 2990, 3210, '-', 11, 211, '+', 'A'*221, ['correct direction']),
]


def lookup_names(ans_set):
    return set(cp.inverse_answers[a] for a,_ in ans_set)


def test_vectors():
    for vec in VECTORS:
        cds_open, cds_close, cds_dir, read_open, read_close, read_dir, pred_start, pred_end, pred_dir, read_seq, expect_names = vec
        answers = cp.check_pred(cds_open, cds_close, cds_dir, read_open, read_close, read_dir, pred_start, pred_end, pred_dir, read_seq)
        names = lookup_names(answers)
        for n in expect_names:
            assert n in names

