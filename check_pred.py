# Public stuff to be used from this module:
# -  answers
# -  inverse_answers
# -  check_pred(cds_open, cds_close, cds_direction,
#              read_open, read_close, read_direction,
#              pred_start, pred_end, pred_direction)

# See examples in main at end of file.

# Glossary: pred_cds_start is the cds coordinates of the start of the prediction. This may come after the pred_cds_end if the read aligns to the genome in reverse direction.


# -----------------------------------------------------------
# Dictionaries converting strings to ints and vice versa
# (for memory efficiency)
# These describe the possibilities for whether a prediction is (partially)
# correct or not for a CDS, and in what way.
# We look at direction, frame, and if possible, start and stop

answers = {
    "correct start": 0,
    "alternative start": 1,
    "middle-alternative start": 2, # can't tell if alternative or correct
    "incorrect start": 3,
    "correct stop": 4,
    "alternative stop": 5,
    "middle-alternative stop": 6, # can't tell if alternative or correct
    "incorrect stop": 7,
    "correct frame": 8,
    "incorrect frame": 9,
    "correct direction": 10,
    "incorrect direction": 11,
}

inverse_answers = {v: k for k, v in answers.items()}

#------------------------------------------------------------

# Use the dictionaries to add an answer to the answers so far

def _add_answer(string_answer, answers_so_far):
     answers_so_far.add(answers[string_answer])

     
#------------------------------------------------------------

# This function is used by all categories that have the correct direction,
# to accumulate the possible answers for a cds-read-pred. 
# Returns a set of answer_vals

def _collect_answers(read_start, read_end, pred_start, pred_end, start_diff, stop_diff, read_captures_cds_start, read_captures_cds_end):

    answer_vals = set()
    
    # We know that the caller of this function has already checked that this
    # prediction is in the correct direction 
    _add_answer("correct direction", answer_vals)
    
    if read_captures_cds_start:
        if start_diff == 0:
            _add_answer("correct start", answer_vals)
            _add_answer("correct frame", answer_vals)
        elif start_diff % 3 == 0:
            _add_answer("alternative start", answer_vals)
            _add_answer("correct frame",     answer_vals)
        else: # wrong frame
            _add_answer("incorrect start", answer_vals)
            _add_answer("incorrect frame", answer_vals)

    if read_captures_cds_end:
        if stop_diff == 0:
            _add_answer("correct stop",  answer_vals)
            _add_answer("correct frame", answer_vals)
        elif stop_diff % 3 == 0:
            _add_answer("alternative stop", answer_vals)
            _add_answer("correct frame",    answer_vals)
        else: # wrong frame
            _add_answer("incorrect stop",  answer_vals)
            _add_answer("incorrect frame", answer_vals)

    if not (read_captures_cds_end or read_captures_cds_start): # middle
        # Check frame
        if start_diff % 3 == 0: # start_diff can't be 0 here
            _add_answer("correct frame", answer_vals)
            if abs(pred_start - read_start) > 3:
                # doesn't fill the read
                _add_answer("alternative start", answer_vals)
            else:
                # we can't tell if this would be correct start or alternative stop
                _add_answer("middle-alternative start", answer_vals)
        if stop_diff % 3 == 0: # stop_diff can't be 0 here
            _add_answer("correct frame", answer_vals)
            if abs(pred_end - read_end) > 3:
                # doesn't fill the read
                _add_answer("alternative stop", answer_vals)
            else:
                # we can't tell if this would be correct stop or alternative stop
                _add_answer("middle-alternative stop", answer_vals)
        else: # wrong frame
            _add_answer("incorrect start", answer_vals)
            _add_answer("incorrect stop",  answer_vals)
            _add_answer("incorrect frame", answer_vals)

    return answer_vals


#------------------------------------------------------------

# This function determines the category and sets the parameters
# accordingly, then calls collect_answers if correct direction.
# Returns a set of answer_vals.

def check_pred(cds_open, cds_close, cds_direction,
               read_open, read_close, read_direction,
               pred_start, pred_end, pred_direction):

    category = (cds_direction, read_direction, pred_direction)

    if category == ('+','+','+'):
        print("category 1")
        read_start = 1
        read_end   = read_close - (read_open - 1)

        pred_cds_start = read_open + (pred_start - 1)
        pred_cds_end   = read_open + (pred_end - 1) 

        start_diff = pred_cds_start - cds_open
        stop_diff  = pred_cds_end - cds_close

        read_captures_cds_start = read_open <= cds_open
        read_captures_cds_end   = read_close >= cds_close

        answer_vals = _collect_answers(read_start, read_end, pred_start, pred_end, start_diff, stop_diff, read_captures_cds_start, read_captures_cds_end)
        

    elif category == ('+','-','-'):
        print("category 3")
        read_start = 1
        read_end   = read_close - (read_open - 1)

        pred_cds_start = read_close - (pred_end - 1)
        pred_cds_end   = read_close - (pred_start - 1)
        
        start_diff = pred_cds_start - cds_open
        stop_diff  = pred_cds_end - cds_close
        
        read_captures_cds_start = read_open <= cds_open
        read_captures_cds_end   = read_close >= cds_close

        answer_vals = _collect_answers(read_start, read_end, pred_cds_start, pred_cds_end, start_diff, stop_diff, read_captures_cds_start, read_captures_cds_end)

        
    elif category == ('-','+','-'):
        print("category 6")
        read_start = 1
        read_end   = read_close - (read_open - 1)

        pred_cds_start = read_open + (pred_start - 1)
        pred_cds_end   = read_open + (pred_end - 1)

        start_diff = pred_cds_end - cds_close
        stop_diff  = pred_cds_start - cds_open

        read_captures_cds_start = read_close >= cds_close
        read_captures_cds_end = read_open <= cds_open

        answer_vals = _collect_answers(read_start, read_end, pred_start, pred_end, start_diff, stop_diff, read_captures_cds_start, read_captures_cds_end)

        
    elif category == ('-','-','+'):
        print("category 8")
        read_start = 1
        read_end   = read_close - (read_open - 1)

        pred_cds_start = read_close - (pred_start - 1)
        pred_cds_end   = read_close - (pred_end - 1)
        
        start_diff = pred_cds_start - cds_close
        stop_diff  = pred_cds_end - cds_open
        
        read_captures_cds_start = read_close >= cds_close
        read_captures_cds_end   = read_open <= cds_open
        
        answer_vals = _collect_answers(read_start, read_end, pred_cds_start, pred_cds_end, start_diff, stop_diff, read_captures_cds_start, read_captures_cds_end)


    else:
        # All other categories are incorrect strand/direction,
        # so don't need to inspect further
        print("Not a good category")
        answer_vals = set()
        _add_answer("incorrect direction", answer_vals)
        

    return answer_vals



#--------------------------------------------------------------------        

# Examples of how to use this check_pred function

def main():
    answer_vals = check_pred(686,1828,'+', 596,811,'+', 91,213,'+') # cat 1, left
    print(answer_vals) # correct start
    for a in answer_vals:
        print(inverse_answers[a])
    print()
    
    answer_vals = check_pred(686,1828,'+', 614,821,'-', 2,136,'-') # cat 3, left
    print(answer_vals) # correct start
    for a in answer_vals:
        print(inverse_answers[a])
    print()
    
    answer_vals = check_pred(12701,13564,'-', 12649,12858,'+', 53,208,'-') # cat 6, left
    print(answer_vals) # correct stop
    for a in answer_vals:
        print(inverse_answers[a])
    print()
    
    answer_vals = check_pred(12701,13564,'-', 12612,12816,'-', 2,203,'+') # cat 8, left
    print(answer_vals) # alternative stop
    for a in answer_vals:
        print(inverse_answers[a])
    print()

    answer_vals = check_pred(1,2,'-', 1,2,'-', 1,2,'-') # incorrect direction
    print(answer_vals) # incorrect direction
    for a in answer_vals:
        print(inverse_answers[a])
    print()


if __name__ == "__main__": 
    main()
