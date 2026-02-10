import importlib.util
import os
import tempfile
import csv

spec = importlib.util.spec_from_file_location('gen','/home/nick/Git/rORForise/src/rORForise/generate_testing_pred.py')
gen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gen)

spec2 = importlib.util.spec_from_file_location('cp','/home/nick/Git/rORForise/src/rORForise/check_pred.py')
cp = importlib.util.module_from_spec(spec2)
spec2.loader.exec_module(cp)


def test_map_cds_coord_to_pred_forward_and_reverse():
    # create a simple read that spans genome coords 100..249 (length 150)
    rs, re = 100, 249
    # CDS start at 110, end at 200 on + strand
    cds_start, cds_end, cds_strand = 110, 200, '+'
    # generate read_seq
    read_seq = gen.generate_sequence(re - rs + 1)

    # test mapping for + read
    pstart, pend = gen.map_cds_coord_to_pred if hasattr(gen, 'map_cds_coord_to_pred') else (None,None)
    # call via generator's wrapper: use a predict invocation
    pred_start, pred_end, pred_strand, seq, expected = gen.generate_prediction(rs, re, '+', cds_start, cds_end, cds_strand, read_seq, scenario='correct_start')
    # when read is +, pred_start should be close to cds_start - rs + 1
    assert pred_start >= 1 and pred_start <= (re - rs + 1)

    # test mapping for reverse read
    pred_start_r, pred_end_r, pred_strand_r, seqr, expectedr = gen.generate_prediction(rs, re, '-', cds_start, cds_end, cds_strand, read_seq, scenario='correct_start')
    assert pred_start_r >= 1 and pred_start_r <= (re - rs + 1)


def test_check_pred_consistency_small():
    # small sanity: generate 10 reads and check check_pred can be called without errors
    for i in range(10):
        cds_start, cds_end, cds_strand = 1000 + i*1000, 1000 + i*1000 + 300, '+' if i % 2 == 0 else '-'
        read_start, read_end, read_strand = gen.generate_read_mapping(cds_start, cds_end, cds_strand, read_length=150, scenario='correct_start')
        read_seq = gen.generate_sequence(read_end - read_start + 1)
        pred_start, pred_end, pred_strand, seq, expected = gen.generate_prediction(read_start, read_end, read_strand, cds_start, cds_end, cds_strand, read_seq, scenario='correct_start')
        # call cp.check_pred and ensure returns iterable of tuples
        answers = cp.check_pred(cds_start, cds_end, cds_strand, read_start, read_end, read_strand, pred_start, pred_end, pred_strand, read_seq)
        assert hasattr(answers, '__iter__')
        for a in answers:
            assert isinstance(a, tuple)

