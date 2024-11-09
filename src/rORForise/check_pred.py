# Public stuff to be used from this module:
# -  answers
# -  inverse_answers
# -  check_pred(cds_open, cds_close, cds_direction,
#              read_open, read_close, read_direction,
#              pred_start, pred_end, pred_direction)

# See examples in main at end of file.

# Glossary: pred_cds_start is the cds coordinates of the start of the prediction.
# This may come after the pred_cds_end if the read aligns to the genome in reverse direction.

# This code assumes all positions are 1-based.


#------------------------------------------------------------

start_codons = ["ATG", "GTG", "TTG"]
stop_codons  = ["TAA", "TAG", "TGA"]

def is_ok_start_codon(codon):
     # return True
     return (codon in start_codons)

def is_ok_stop_codon(codon):
     # return True
     return (codon in stop_codons)

# -----------------------------------------------------------
# Dictionaries converting strings to ints and vice versa
# (for memory efficiency)
# These describe the possibilities for whether a prediction is (partially)
# correct or not for a CDS, and in what way.
# We look at direction, frame, and if possible, start and stop

answers = {
    "correct start": 0,
    "correct stop": 1,
    "correct frame": 2,
    "correct direction": 3,
    "incorrect stop": 4,
    "incorrect start": 5,
    "incorrect frame": 6,
    "incorrect direction": 7,
    "prediction ends before cds starts": 8,
    "prediction starts after cds ends": 9,
    "alternative start": 10,
    "middle or alternative start": 11,
    "alternative stop": 12,
    "middle or alternative stop": 13,
    "middle": 14,
}

     


inverse_answers = {v: k for k, v in answers.items()}


#------------------------------------------------------------

# Use the dictionaries to add an answer to the answers so far

def _add_answer(string_answer, codon, answers_so_far):
     answers_so_far.add((answers[string_answer], codon))

     
#------------------------------------------------------------

# This really needs a fuller working
def _rev_comp(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))

#------------------------------------------------------------

# This function is used by all categories that have the correct direction,
# to accumulate the possible answers for a cds-read-pred. 
# Returns a set of answer_vals
def _collect_answers(read_start, read_end, pred_start, pred_end, start_diff, stop_diff, read_captures_cds_start, read_captures_cds_end, pred_before_start, pred_after_end, pred_start_codon, pred_end_codon):

    answer_vals = set()
    
    # We know that the caller of this function has already checked that this
    # prediction is in the correct direction for the gene
    _add_answer("correct direction", None, answer_vals)
    
    if read_captures_cds_start:
        if start_diff == 0:
            _add_answer("correct start", pred_start_codon, answer_vals)
            _add_answer("correct frame", None, answer_vals)
        elif start_diff % 3 == 0:
            _add_answer("correct frame", None, answer_vals)
            if is_ok_start_codon(pred_start_codon):
               _add_answer("alternative start", pred_start_codon, answer_vals)
            else:
               _add_answer("incorrect start", pred_start_codon, answer_vals)
        else: # wrong frame
            _add_answer("incorrect start", None, answer_vals)
            _add_answer("incorrect frame", None, answer_vals)
        if pred_before_start:
            _add_answer("prediction ends before cds starts", None, answer_vals)
            
    if read_captures_cds_end:
        if stop_diff == 0:
            _add_answer("correct stop", pred_end_codon, answer_vals)
            _add_answer("correct frame", None, answer_vals)
        elif stop_diff % 3 == 0:
            _add_answer("correct frame", None, answer_vals)
            if is_ok_stop_codon(pred_end_codon):
               _add_answer("alternative stop", pred_end_codon, answer_vals)
            else:
               _add_answer("incorrect stop", pred_end_codon, answer_vals)
        else: # wrong frame
            _add_answer("incorrect stop", None, answer_vals)
            _add_answer("incorrect frame", None, answer_vals)
        if pred_after_end:
             _add_answer("prediction starts after cds ends", None, answer_vals)

    if not (read_captures_cds_end or read_captures_cds_start): # middle
        # Check start frame
        if start_diff % 3 == 0: # start_diff can't be 0 here
            _add_answer("correct frame", None, answer_vals)
            if abs(pred_start - read_start) >= 3:
                # doesn't fill the read
                if is_ok_start_codon(pred_start_codon):
                   _add_answer("alternative start", pred_start_codon, answer_vals)
                else:
                   _add_answer("incorrect start", pred_start_codon, answer_vals)
            else:
                # fills as much of the read as it can
                if is_ok_start_codon(pred_start_codon):
                   _add_answer("middle or alternative start", pred_start_codon, answer_vals)
                else:
                   _add_answer("middle", None, answer_vals)
        # Check stop frame
        if stop_diff % 3 == 0: # stop_diff can't be 0 here
            _add_answer("correct frame", None, answer_vals)
            if abs(pred_end - read_end) >= 3: 
                # doesn't fill the read
                if is_ok_stop_codon(pred_end_codon):
                   _add_answer("alternative stop", pred_end_codon, answer_vals)
                else:
                   _add_answer("incorrect stop", pred_end_codon, answer_vals) # should really remove any existing "middle" answer
            else:
                # fills as much of the read as it can
                if is_ok_stop_codon(pred_end_codon):
                   _add_answer("middle or alternative stop", pred_end_codon, answer_vals)
                else:
                   _add_answer("middle", None, answer_vals)
        else: # wrong frame
            _add_answer("incorrect start", None, answer_vals)
            _add_answer("incorrect stop",  None, answer_vals)
            _add_answer("incorrect frame", None, answer_vals)

    return answer_vals


#------------------------------------------------------------

# This function determines the category and sets the parameters
# accordingly, then calls collect_answers if correct direction.
# Returns a set of answer_vals.
# Assumes all positions are 1-based.

def check_pred(cds_open, cds_close, cds_direction,
               read_open, read_close, read_direction,
               pred_start, pred_end, pred_direction,
               read_seq):

    category = (cds_direction, read_direction, pred_direction)

    if category == ('+','+','+'):
        #print("category 1")
        read_start = 1
        read_end = read_close - (read_open-1)

        pred_cds_start = read_open + (pred_start - 1)
        pred_cds_end   = read_open + (pred_end - 1)

        start_diff = pred_cds_start - cds_open
        stop_diff  = pred_cds_end - cds_close

        read_captures_cds_start = read_open <= cds_open
        read_captures_cds_end = read_close >= cds_close

        pred_before_start_of_cds = pred_cds_end < cds_open
        pred_after_end_of_cds = pred_cds_start > cds_close 

        pred_start_codon = read_seq[pred_start-1:pred_start-1+3]
        pred_end_codon   = read_seq[pred_end-3:pred_end]

        #print(read_start, read_end, pred_start, pred_end, start_diff, stop_diff, read_captures_cds_start, read_captures_cds_end, pred_before_start_of_cds, pred_after_end_of_cds, pred_start_codon, pred_end_codon)
        #print(read_seq)
        
        answer_vals = _collect_answers(read_start, read_end, pred_start, pred_end, start_diff, stop_diff, read_captures_cds_start, read_captures_cds_end, pred_before_start_of_cds, pred_after_end_of_cds, pred_start_codon, pred_end_codon)


    elif category == ('+','-','-'):
        #print("category 3")
        read_start = 1
        read_end   = read_close - (read_open-1)

        pred_cds_start = read_close - (pred_end - 1)
        pred_cds_end   = read_close - (pred_start - 1)
        
        start_diff = pred_cds_start - cds_open
        stop_diff  = pred_cds_end - cds_close
        
        read_captures_cds_start = read_open <= cds_open
        read_captures_cds_end   = read_close >= cds_close

        pred_before_start_of_cds = pred_cds_end < cds_open
        pred_after_end_of_cds = pred_cds_start > cds_close 

        # read is '-' so in samfile and tsv file the sequence will be reverse complemented
        pred_start_codon = _rev_comp(read_seq)[read_end-pred_end:read_end-pred_end+3] 
        pred_end_codon   = _rev_comp(read_seq)[read_end-pred_start-3+1:read_end-pred_start+1]

        #print(read_open, read_close, pred_cds_start, pred_cds_end, start_diff, stop_diff, read_captures_cds_start, read_captures_cds_end, pred_before_start_of_cds, pred_after_end_of_cds, pred_start_codon, pred_end_codon)
        #print(read_seq)
        
        answer_vals = _collect_answers(read_open, read_close, pred_cds_start, pred_cds_end, start_diff, stop_diff, read_captures_cds_start, read_captures_cds_end, pred_before_start_of_cds, pred_after_end_of_cds, pred_start_codon, pred_end_codon)

        
    elif category == ('-','+','-'):
        #print("category 6 -+-")
        read_start = 1
        read_end   = read_close - (read_open - 1)

        pred_cds_start = read_open + (pred_start - 1)
        pred_cds_end   = read_open + (pred_end - 1)
 
        start_diff = pred_cds_end - cds_close
        stop_diff  = pred_cds_start - cds_open

        read_captures_cds_start = read_close >= cds_close
        read_captures_cds_end = read_open <= cds_open

        pred_before_start_of_cds = pred_cds_start > cds_close
        pred_after_end_of_cds = pred_cds_end < cds_open 
        
        pred_start_codon = _rev_comp(read_seq[pred_end-3:pred_end])
        pred_end_codon   = _rev_comp(read_seq[pred_start-1:pred_start-1+3])

        #print(read_start, read_end, pred_start, pred_end, start_diff, stop_diff, read_captures_cds_start, read_captures_cds_end, pred_before_start_of_cds, pred_after_end_of_cds, pred_start_codon, pred_end_codon)
        #print(read_seq)

        answer_vals = _collect_answers(read_start, read_end, pred_start, pred_end, start_diff, stop_diff, read_captures_cds_start, read_captures_cds_end, pred_before_start_of_cds, pred_after_end_of_cds, pred_start_codon, pred_end_codon)


    elif category == ('-','-','+'):
        #print("category 8")
        read_start = 1
        read_end   = read_close - (read_open -1 )

        pred_cds_start = read_close - (pred_start - 1)
        pred_cds_end   = read_close - (pred_end - 1)

        start_diff = pred_cds_start - cds_close
        stop_diff  = pred_cds_end - cds_open
        
        read_captures_cds_start = read_close >= cds_close
        read_captures_cds_end   = read_open <= cds_open

        pred_before_start_of_cds = pred_cds_end > cds_close
        pred_after_end_of_cds = pred_cds_start < cds_open

        # read is '-' so in samfile and tsv file the sequence will be reverse complemented
        pred_start_codon = read_seq[pred_start-1:pred_start-1+3] 
        pred_end_codon   = read_seq[pred_end-3:pred_end] 

        #print(read_start, read_end, pred_start, pred_end, start_diff, stop_diff, read_captures_cds_start, read_captures_cds_end, pred_before_start_of_cds, pred_after_end_of_cds, pred_start_codon, pred_end_codon)
        #print(read_seq)

        answer_vals = _collect_answers(read_start, read_end, pred_start, pred_end, start_diff, stop_diff, read_captures_cds_start, read_captures_cds_end, pred_before_start_of_cds, pred_after_end_of_cds, pred_start_codon, pred_end_codon)


    else:
        # All other categories are incorrect strand/direction,
        # so don't need to inspect further
        #print("Not a good category")
        answer_vals = set()
        _add_answer("incorrect direction", None, answer_vals)
        

    return answer_vals



#--------------------------------------------------------------------        

# TBD: Examples of how to use this check_pred function
# All examples have 1-based positions

def main():

     pass
    

if __name__ == "__main__": 
    main()




