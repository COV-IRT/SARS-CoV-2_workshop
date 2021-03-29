#!/usr/bin/env python3
#
# Subsample sequences
#
#
# Author: Josh Cherry

import os, sys, getopt, itertools
import datetime, time, hashlib, random, pickle


# Read fasta (from course.py)
def fasta_read(fname_or_fid):
    """
    Read sequences from a file in FASTA format.
    Returns a dictionary whose values are sequences
    and keys are sequence identifiers.
    
    *fname_or_fid* is a string (file name) or an open file
    """
    seq_lines = {}
    names = []
    if type(fname_or_fid) == str:  # were passed a file name
        fid = open(fname_or_fid)
    else:                          # otherwise, treat as file-like object
        fid = fname_or_fid
    for line in fid:
        if line.startswith('>'):
            name = line[1:].split()[0]
            if name in seq_lines:
                raise RuntimeError('Duplicate name "%s" in fasta file %s'
                                   % (name, fname_or_fid))
            seq_lines[name] = []
            names.append(name)
        else:
            seq_lines[name].append(line.strip())

    rv = {}
    for name in names:
        rv[name] = ''.join(seq_lines[name]) # .upper()

    return [rv[name] for name in names], names


# Crude fasta writing
def writefasta(seqs, ids, fname, comments=None):
    line_len = 60
    if len(seqs) != len(ids):
        raise RuntimeError('lenth of seqs not equal to length of ids')
    with open(fname, 'w') as fid:
        for i in range(len(seqs)):
            if comments and i < len(comments):
                fid.write('>%s  %s\n' % (ids[i], comments[i]))
            else:
                fid.write('>%s\n' % ids[i])
            for pos in range(0, len(seqs[i]), line_len):
                fid.write(seqs[i][pos:pos+line_len]); fid.write('\n')

                

# Group things (into lists) according to the result of supplied function
# (which becomes a dictionary key)
def GroupBy(l, func):
    d = {}
    for el in l:
        key = func(el)
        if not key in d:
            d[key] = []
        d[key].append(el)
    return d

    
# iterator (generator) that "flattens" a list of lists (or other iterables),
# iterating over all the elements of all the elements.
def flatten(l):
    for t in l:
        for x in t:
            yield x   
            

            
# return randomly shuffled version of iterable as a list
def shuffled(it):
    import random
    l = list(it)
    random.shuffle(l)
    return l


# Given a list or tuple, make a dictionary of counts of the keys
def Counts(l):
    rv = {}
    for el in l:
        if not el in rv:
            rv[el] = 0
        rv[el] += 1
    return rv


_times = []
def tick():
    import os
    global _times
    _times.append(os.times())

    
tt_messages = []
def tock(msg=''):
    import os
    global _times
    global tt_messages
    new_times = os.times()
    tic_times = _times.pop()
    indent = (2 * len(_times)) * '  '
    st = '%s%s  %s' % (indent,
                       ' '.join(['%.3f' % (new_times[i] - tic_times[i])
                               for i in [0, 1, 4]]),
                     msg)
    print(st)
    tt_messages.append(st)


def GetInputFname(prompt, type_text=''):
    '''
    Get path to file from user, guaranteeing that it exists
    and is readable, and tolerating blank line input.
    '''
    while True:
        fname = input(prompt)
        if not fname.strip():
            # just white space
            pass
        if fname.strip() == '.':
                try:
                    fname = Gui.FileDialog('r', prompt, type_text)
                except Exception as e:
                    print('Problem running graphical file dialog: %s'
                          % e.args[0])
                    print('Returning to text interface')
                    print()
                    continue
                if fname == None:
                    print('File dialog cancelled; enter file name at prompt')
                    print()
                    continue
                else:
                    return fname
        elif not os.access(fname, os.F_OK):
            print()
            print('%s does not exist; try again' % fname)
        elif not os.access(fname, os.R_OK):
            print()
            print('%s exists but is not readable; try again' % fname)
        else:
            # good fname
            return fname


def GetPrefsWithinCat():
    '''
    Get preference crieteria in order of priority.
    '''
    print('Specify the criteria for choosing within a category, ')
    print('in order of importance.  The first criterion takes priority, ')
    print('so the others only matter when candidates are "tied" for the')
    print('first.  Similarly, the second takes priority over all but ')
    print('the first, and so on.  Randomization applies in all cases.')
    print()
    print('Specify one or more letters (either case), separated by whitespace, commas, or nothing:')
    print()
    print('N  No preferences.  Purely random choice.  Must be only letter specified.')
    print('U  Seek uniformity of date distribution.')
    print('   If used, this must come first in priorities,')
    print('   and the date (in effect, year) must be part of')
    print('   the category definition.')
    print('L  Maximize sequence length '
          '(includes internal gaps, but not leading and trailing gaps)')
    print('D  Completeness of date.  YMD > YM > Y')
    print('M  Presence of month; no preference for YMD over YM.')
    print('   It can be meaningful and useful to provide both D and M,')
    print('   provided that M comes first and something comes between, e.g., MLD.')
    print('G  Minimize number of internal gaps with lengths not divisible by 3 and less than 5')
    print('A  Minimize number of ambiguity characters (N, etc.).')
    print()
    while True:
        line = input('Specify preferences for isolate selection in order of priority: ')
        if not line.strip():
            # just white space
            pass
        else:
            # Get the letters as a string with nothing else
            s = ''.join(line.replace(',', ' ').split()).upper()
            if len(set(s)) < len(s):
                print('Duplicate letters not allowed')
                continue
            if 'N' in s:
                if len(s) > 1:
                    print('N must be the only letter if specified')
                    continue
            else:
                if not set(s).issubset('ULDMGA'):
                    print('Response contained non-allowed character(s): %s'
                            % set(s).difference('ULDMGA'))
                    continue
                if 'U' in s[1:]:
                    print('U must come first if specified.')
                    continue
                if 'D' in s and 'M' in s:
                    idx_d = s.find('D')
                    idx_m = s.find('M')
                    if idx_d < idx_m:
                        print('M is redundant with D unless M comes first (has higher priority)')
                        continue
                    if idx_d == idx_m + 1:
                        print('M is redundant with D unless something comes between them')
                        continue
            return s
            
                
            
def GetYN(prompt, default=None):
    '''
    '''
    while True:
        s = input(prompt)
        if not s:
            # just hit return
            if default is not None:
                return default
        else:
            s_mod = s.strip().lower()
            if s_mod and len(s_mod) <= 3:
                if 'yes'.startswith(s_mod):
                    return True
                if 'no'.startswith(s_mod):
                    return False
        print()
        if default is not None:
            print('please type y(es) or n(o) (or RETURN for %s)'
                  % ('yes' if default else 'no'))
        else:
            print('please type y(es) or n(o)')


def GetOutputFname(prompt, type_text=''):
    '''
    Get path to file from user, tolerating blank line input.
    
    '.' causes file dialog to launch, if possible (requires wxPython)
    '''
    while True:
        fname = input(prompt)
        if not fname.strip():
            # just white space
            pass
        else:
            if fname.strip() == '.':
                try:
                    fname = Gui.FileDialog('w', prompt, type_text)
                except Exception as e:
                    print('Problem running graphical file dialog: %s' % e.args[0])
                    print('Returning to text interface')
                    print()
                    continue
                if fname == None:
                    print('File dialog cancelled; enter file name at prompt')
                    print()
                    continue
            # good fname
            return fname


def GetFieldChoices(num_fields):
    '''
    Get one or more field indices.
    '''
    while True:
        s = input('Enter one or more field numbers '
                  'separated by commas or whitespace, '
                  'or S to load spreadsheet: ')
        if not s.strip():
            # just white space
            pass
        else:
            if s.strip().upper() == 'S':
                return 'S'
            field_strs = s.replace(',', ' ').split()
            field_indices = []
            for fs in field_strs:
                try:
                    n = int(fs)
                except ValueError:
                    print('"%s" is not a valid integer' % fs)
                    break
                if n < 0 or n >= num_fields:
                    print('%d out of range' % n)
                    break
                field_indices.append(n)
            else:  # unless break from above
                if len(set(field_indices)) < len(field_indices):
                    print('Error: duplicate field indices')
                else:
                    return field_indices


def GetPosInt(prompt):
    '''
    Get a positive integer from the user.
    '''
    while True:
        s = input(prompt)
        if not s.strip():
            # just white space
            pass
        else:
            try:
                n = int(s)
            except ValueError:
                print('"s" not a valid integer' % s)
                continue
            if n <= 0:
                print('Number must be greater than zero')
                continue
            print()
            return n

            
def GetSingleCharacter(prompt, default):
    '''
    Get single character from user, with default value
    '''
    while True:
        ch = input(prompt)
        if not ch:
            # just hit return
            return default
        elif len(ch) > 1:
            print()
            print('%s invalid; must be single character' % ch)
        else:
            return ch


def StrFTime(fmt, t=time.localtime()):
    '''
    Fix for long (non-abbreviated) time zone string under Windows
    '''
    s = time.strftime(fmt, t)
    if '%Z' in fmt:
        tz = time.strftime('%Z', t)
        if len(tz.split()) == 3:
            abbr = ''.join([x[0] for x in tz.split()])
            s = s.replace(tz, abbr)
    return s
                  
           
def TimeRange(st_early, st_late):
    '''
    Acts on time.struct_time instances.
    Makes a string represnting the range of time.
    '''
    if st_early.tm_year != st_late.tm_year:
        return '%s-%s' % (StrFTime('%X %Z %d %b %Y', st_early), StrFTime('%X %Z %d %b %Y', st_late))
    elif st_early.tm_mon != st_late.tm_mon or st_early.tm_mday != st_late.tm_mday:
        # Not the same day, but same year
        return '%s-%s' % (StrFTime('%X %Z %d %b', st_early), StrFTime('%X %Z %d %b %Y', st_late))
    elif st_early.tm_isdst != st_late.tm_isdst:
        # Unlikely: same day, but DST has changed
        return '%s-%s' % (StrFTime('%X %Z', st_early), StrFTime('%X %Z %d %b %Y', st_late))
    elif st_early.tm_hour != st_late.tm_hour:
        # Same day, and DST has not changed.
        return '%s-%s' % (StrFTime('%X', st_early), StrFTime('%X %Z %d %b %Y', st_late))
    elif st_early.tm_min != st_late.tm_min:
        return '%s-%s' % (StrFTime('%X', st_early), StrFTime(':%M:%S %Z %d %b %Y', st_late))
    else:
        return '%s-%s' % (StrFTime('%X', st_early), StrFTime(':%S %Z %d %b %Y', st_late))
            
            
def Quit(val):
    #Gui.fini()
    if sys.platform.startswith('win32'):
        if 'PROMPT' not in os.environ:  # Not run from command prompt?
            tmp = input('Press return to exit')
    exit(val)

    
def ShowExampleTagVals(names, delimiter):
    '''
    Show example tag field values to aid user
    '''
    def processTag(tag):
        '''
        If tag length above threshold, replace with
        truncation plus '...'
        '''
        if len(tag) > 30:
            return tag[:27] + '...'
        return tag
    print('Example tag values:')
    print()
    tags1 = names[0].split(delimiter)
    tags2 = names[-1].split(delimiter)
    if len(tags1) != len(tags2):
        print('Number of fields in first and last entries not equal')
        Quit(1)
    for i in range(len(tags1)):
        print('%d  %-30s  %-30s'
              % (i, processTag(tags1[i]), processTag(tags2[i])))
    print()    


def FindDateField(names, delim):
    def CouldBeDate(s):
        return (len(s) >= 4 and s[:4].isdecimal()
                and set(s).issubset('0123456789-')
                and (len(s) < 5 or s[4] == '-'))
    fields = [n.split(delim) for n in names]
    date_fields = [i for i in range(len(fields[0]))
                    if all(CouldBeDate(t[i]) for t in fields)]
    if len(date_fields) != 1:
        raise RuntimeError('Could not determine date field')
    return date_fields[0]


def FieldTuple(name, delim, fields, date_index):
    l = name.split(delim)
    if date_index in fields:
        l[date_index] = l[date_index][:4]
    return tuple(l[i] for i in fields)


def DateCompleteness(name, delim, idx):
    '''
    0 if only year
    1 if year + month
    2 if year + month + date
    '''
    s = name.split(delim)[idx]
    if set(s[4:]).issubset('-'):
        return 0
    if s[4] == '-' and s[5:7].isdecimal():
        # at least have month
        if set(s[7:]).issubset('-'):
            return 1
        if s[7] == '-' and s[8:10].isdecimal() and len(s) == 10:
            return 2
    raise ValueError('Invalid date string: %s' % s)  # Should have returned


def NumGapsIndivisibleThreeLT5(s):
    '''
    For covid, count only gaps of length < 5
    '''
    s = s.strip('-')  # lose leading and terminal "gaps"
    # replace everything but '-' with ' '
    for ch in set(s):
        if ch != '-':
            s = s.replace(ch, ' ')
    gaps = s.split()
    return len([g for g in gaps if len(g) % 3 != 0 and len(g) < 5])
    
     

all_distrs_cache = {}
scores_cache = {}
def AllDistrs(n, ncat, cat_maxes):
    if n == 0:
        return set([ncat * (0,)])
    if ncat == 1:
        return set([(n,)])
    if type(cat_maxes) != tuple:
        cat_maxes = tuple(cat_maxes)
    if (n, ncat, cat_maxes) in all_distrs_cache:
        return all_distrs_cache[(n, ncat, cat_maxes)]
    res = []
    half = ncat // 2
    for i in range(max(0, n - sum(cat_maxes[half:])),
                   min(n, sum(cat_maxes[:half])) + 1):
        first = AllDistrs(i, len(cat_maxes[:half]), cat_maxes[:half])
        second = AllDistrs(n - i, len(cat_maxes[half:]), cat_maxes[half:])
        res.extend(f + s for f in first for s in second)
    all_distrs_cache[(n, ncat, cat_maxes)] = res
    return res


all_best_distrs_cache = {}
def AllBestDistrs(n, ncat, cat_maxes):
    '''
    Returns only "best" distributions with regard to
    numbers on same date
    ''' 
    if n == 0:
        return [ncat * (0,)]
    if ncat == 1:
        return [(n,)]
    if type(cat_maxes) != tuple:
        cat_maxes = tuple(cat_maxes)
    if (n, ncat, cat_maxes) in all_best_distrs_cache:
        return all_best_distrs_cache[(n, ncat, cat_maxes)]
    res = []
    best = None  # will contain sorted best distr
    half = ncat // 2
    for i in range(max(0, n - sum(cat_maxes[half:])),
                   min(n, sum(cat_maxes[:half])) + 1):
        first = AllBestDistrs(i, len(cat_maxes[:half]), cat_maxes[:half])
        second = AllBestDistrs(n - i, len(cat_maxes[half:]), cat_maxes[half:])
        if first and second:
            cand_best = sorted(first[0] + second[0], reverse=True)
            if cand_best == best:
                res.extend(f + s for f in first for s in second)
            elif best is None or cand_best < best:
                res = [f + s for f in first for s in second]
                best = cand_best
    all_best_distrs_cache[(n, ncat, cat_maxes)] = res
    return res


# Load some pre-computed results for speed, if they can be found
mec_basename = 'subsamp_data.pic'
prog_dir = os.path.dirname(os.path.realpath(__file__))
mec_fname = os.path.join(prog_dir, mec_basename)
if os.path.exists(mec_fname):
    min_energy_cache = pickle.load(open(mec_fname, 'rb'))
else:
    msg = ''
    msg += '\n'
    msg += '** WARNING: file %s not found in %s\n' % (mec_basename, prog_dir)
    msg += '** If the \'U\' criterion is used, computations may be\n'
    msg += '** slower than they could be.\n'
    msg += '** Results will be correct nonetheless.\n'
    msg += '\n'
    sys.stderr.write(msg)
    min_energy_cache = {}

def MinEnergy(n, ncat, cat_maxes):
    if n == 0:
        return set([ncat * (0,)])
    if ncat == 1:
        return set([(n,)])
    if type(cat_maxes) != tuple:
        cat_maxes = tuple(cat_maxes)
    if (n, cat_maxes) in min_energy_cache:
        return min_energy_cache[(n, cat_maxes)]
    half = ncat // 2
    en_min = Energy((n,) + (ncat-1)*(0,))
    res_min = []
    for i in range(max(0, n - sum(cat_maxes[half:])),
                   min(n, sum(cat_maxes[:half])) + 1):
        first = AllDistrs(i, len(cat_maxes[:half]),
                          cat_maxes[:half])
        second = AllDistrs(n - i, len(cat_maxes[half:]),
                           cat_maxes[half:])
        first = list(first)
        second = list(second)
        energy_first = [Energy(f, ncat, n) for f in first]
        energy_second = [Energy(s, ncat, n) for s in second]
        #for i1, d1 in enumerate(first):
        for i1 in sorted(range(len(energy_first)),
                         key=lambda x: energy_first[x]):
            if energy_first[i1] > en_min:
                continue
            d1 = first[i1]
            for i2 in sorted(range(len(energy_second)),
                             key=lambda x: energy_second[x]):
                if energy_first[i1] + energy_second[i2] > en_min:
                    continue
                d2 = second[i2]
                d = d1 + d2
                en = Energy(d)
                if en <= en_min:
                    if en < en_min:
                        en_min = en
                        res_min = [d]
                    else:
                        res_min.append(d)
    min_energy_cache[n, cat_maxes] = res_min
    return res_min


def MinEnergy2(n, ncat, cat_maxes):
    if n == 0:
        return set([ncat * (0,)])
    if ncat == 1:
        return set([(n,)])
    if type(cat_maxes) != tuple:
        cat_maxes = tuple(cat_maxes)
    if (n, cat_maxes) in min_energy_cache:
        return min_energy_cache[(n, cat_maxes)]
    half = ncat // 2
    en_min = Energy((n,) + (ncat-1)*(0,))
    res_min = []
    for i in range(max(0, n - sum(cat_maxes[half:])),
                   min(n, sum(cat_maxes[:half])) + 1):
        first = AllBestDistrs(i, len(cat_maxes[:half]),
                              cat_maxes[:half])
        second = AllBestDistrs(n - i, len(cat_maxes[half:]),
                               cat_maxes[half:])
        first = list(first)
        second = list(second)
        energy_first = [Energy(f, ncat, n) for f in first]
        energy_second = [Energy(s, ncat, n) for s in second]
        #for i1, d1 in enumerate(first):
        for i1 in sorted(range(len(energy_first)),
                         key=lambda x: energy_first[x]):
            if energy_first[i1] > en_min:
                continue
            d1 = first[i1]
            #for i2, d2 in enumerate(second):
            for i2 in sorted(range(len(energy_second)),
                             key=lambda x: energy_second[x]):
                if energy_first[i1] + energy_second[i2] > en_min:
                    continue
                d2 = second[i2]
                d = d1 + d2
                en = Energy(d)
                if en <= en_min:
                    if en < en_min:
                        en_min = en
                        res_min = [d]
                    else:
                        res_min.append(d)
    min_energy_cache[n, cat_maxes] = res_min
    return res_min


def EnergyKernel(t, len_tot):
    lcm = LCM(len_tot - 1)
    kern = (len_tot - len(t)) * [0]
    for i, n in enumerate(t):
        if not n:
            continue
        for j in range(len(kern)):
            d = len(t) - i + j
            kern[j] += lcm // d * t[i]
    return kern


def combinations(items, k):
    n = len(items)
    if n == k:
        yield items
    elif k==0:
        yield tuple()
    else:
        for i in range(n - k + 1):
            for t in combinations(items[i+1:], k - 1):
                yield (items[i],) + t


def nchoosek(n, k):
    k = max(k, n - k)
    num = 1
    for i in range(k + 1, n + 1):
        num *= i
    den = 1
    for i in range(1, n - k + 1):
        den *= i
    res = num // den
    if num / den != res:
        raise RuntimeError('%d != %f; num=%d, den=%d'
                           % (res, num / den, num, den))
    return res


def inside_out(it):
    l = list(it)
    l1 = l[:len(l)//2]
    l2 = l[len(l)//2:]
    for i in range(len(l2)):
        yield l2[i]
        if i < len(l1):
            yield l1[-(i+1)]

            
def MinEnergy3(n, cat_maxes):
    ncat = len(cat_maxes)
    if n == 0:
        return set([ncat * (0,)])
    if ncat == 1:
        return set([(n,)])
    if type(cat_maxes) != tuple:
        cat_maxes = tuple(cat_maxes)
    if sum(cat_maxes) == n:
        return [cat_maxes]
    if (n, cat_maxes) in min_energy_cache:
        return min_energy_cache[(n, cat_maxes)]

    mx = max(cat_maxes)
    poss = tuple(i for i, x in enumerate(cat_maxes)  # places where one might go
                 if x == mx)
    base = tuple(min(x, mx - 1) for x in cat_maxes)  # definitely placed
    
    k = n - sum(base)  # remaining to be placed
    res_min = [tuple(base[i] + (i in poss[:k]) for i in range(len(base)))]
    en_min = Energy(res_min[0]) + 1

    if verbose:
        print('Choosing %d from %d positions' % (k, len(poss)))

    ncombs = nchoosek(len(poss), k)

    if ncombs < 20000:
        for t in combinations(poss, k):
            cand = list(base)
            for i in t:
                cand[i] += 1
            cand = tuple(cand)
            en = Energy(cand)
            if en <= en_min:
                if en < en_min:
                    en_min = en
                    res_min = []
                res_min.append(cand)
    else:
        # This should be faster than above, at least for hard cases

        # get an upper bound on best energy
        ens = [Energy(MinEnergyHeur(n, cat_maxes))
               for i in range(50)]
        en_ub = min(ens)

        # find candidates whose lower bound is <= upper bound
        half = len(poss) // 2
        poss1 = poss[:half]
        poss2 = poss[half:]
        ncands = 0
        
        en_min = en_ub + 1
        res_min = []
        for k1 in inside_out(range(max(k - len(poss2), 0),
                                   min(k, len(poss1)) + 1)):
            k2 = k - k1
            cands1 = [tuple(base[i] + (i in t)
                            for i in range(poss1[-1]+1))
                      for t in combinations(poss1, k1)]
            cands2 = [tuple(base[i] + (i in t)
                            for i in range(poss1[-1]+1, len(cat_maxes)))
                      for t in combinations(poss2, k2)]
            ens1 = [Energy(t, n_tot=n, len_tot=ncat) for t in cands1]
            ens2 = [Energy(t, n_tot=n, len_tot=ncat) for t in cands2]

            far2 = tuple(base[i] + (i in poss2[-k2:])
                         for i in range(poss1[-1]+1, len(cat_maxes)))
            en_far2 = Energy(far2, n_tot=n, len_tot=ncat)
            min_en_comb1 = [Energy(t + far2) - ens1[i] - en_far2
                            for i, t in enumerate(cands1)]

            far1 = tuple(base[i] + (i in poss1[:k1])
                         for i in range(poss1[-1]+1))
            en_far1 = Energy(far1, n_tot=n, len_tot=ncat)
            min_en_comb2 = [Energy(far1 + t) - ens2[i] - en_far1
                            for i, t in enumerate(cands2)]

            en1_min = min(ens1)
            en2_min = min(ens2)
            cands1 = [(i, c) for i, c in enumerate(cands1)
                      if min_en_comb1[i] + ens1[i] + en2_min <= en_ub]
            cands2 = [(i, c) for i, c in enumerate(cands2)
                      if min_en_comb2[i] + ens2[i] + en1_min <= en_ub]
            
            for i1, c1 in cands1:
                kern = EnergyKernel(c1, ncat)
                en_c1_base = sum(k * base[len(c1)+i]
                                 for i, k in enumerate(kern))
                for i2, c2 in cands2:
                    # conjunction is faster than using max()?
                    if (ens1[i1] + ens2[i2] + min_en_comb1[i1] <= en_min
                        and
                        ens1[i1] + ens2[i2] + min_en_comb2[i2] <= en_min):
                        ncands += 1
                        en = ens1[i1] + ens2[i2] + en_c1_base
                        for i, c in enumerate(c2):
                            if c == mx:
                                en += kern[i]
                                if en > en_min:
                                    en = None
                                    break
                        if en is None:
                            continue
                        if en <= en_min:
                            cand = c1 + c2
                            if en < en_min:
                                en_min = en
                                res_min = []
                            res_min.append(cand)

        if verbose:
            print('Reduced from %d to %d candidates'
                  % (ncombs, ncands))

    min_energy_cache[n, cat_maxes] = res_min
    return res_min


def Minimize(t, cat_maxes):
    t = list(t)
    en_min = Energy(tuple(t))
    while True:
        for i in range(len(t)):
            if t[i] == 0:
                continue
            for j in range(len(t)):
                if j == i or t[j] == cat_maxes[j] or t[j] >= t[i]:
                    continue
                cand = list(t)
                cand[i] -= 1
                cand[j] += 1
                en_cand = Energy(tuple(cand))
                if en_cand < en_min:
                    t = cand
                    en_min = en_cand
                    break
            else:
                continue
            break
        else:
            # occupied = [i for i, n in enumerate(t) if n]
            # for i in occupied[:-1]:
            #     for j in occupied:
            #         if j <= i or j == len(t) - 1:
            #             continue
            break
    return t


def MinEnergyHeur(n, cat_maxes):
    ncat = len(cat_maxes)
    if n == 0:
        return set([ncat * (0,)])
    if ncat == 1:
        return set([(n,)])
    if type(cat_maxes) != tuple:
        cat_maxes = tuple(cat_maxes)

    # Initialize distribution
    dist_ini = len(cat_maxes) * [0]
    for it in range(n):
        mn = min(c for i, c in enumerate(dist_ini)
                 if c < cat_maxes[i])
        tmp = [i for i, c in enumerate(dist_ini)
               if c == mn and c < cat_maxes[i]]
        idx = random.choice(tmp)
        dist_ini[idx] += 1

    res = Minimize(dist_ini, cat_maxes)
    res = tuple(res)
    return res


def Primes(n):
    '''
    Returns all primes <= n
    '''
    primes = []
    for i in range(2, n+1):
        if all([i % p for p in primes]):
            primes.append(i)
    return primes


lcms = {}
def LCM(n):
    '''
    Least common multiple of {1, 2, 3...n}
    '''
    if n in lcms:
        return lcms[n]
    primes = Primes(n)
    powers = []
    for p in primes:
        for i in itertools.count():
            if p**i > n:
                powers.append(i-1)
                break
    res = 1
    for i, p in enumerate(primes):
        res *= p**powers[i]
    lcms[n] = res
    return res


energy_cache = {}
def Energy(t, len_tot=None, n_tot=None):
    '''
    n_tot is the total number that will be chosen, which may
    be greater than sum(t) if t contains just part of the
    distribution.  It affects the energy associated
    with elements of t that are greater than one.
    By default (i.e., if None), it is sum(t).

    Similarly, len_tot is the total length of the distribution,
    and affects the scale of the result (lcm.  By default it
    is len(t)
    '''
    if n_tot is None:
        n_tot = sum(t)
    if len_tot is None:
        len_tot = len(t)
    if (t, len_tot, n_tot) in energy_cache:
        return energy_cache[(t, len_tot, n_tot)]
    lcm = LCM(len_tot - 1)
    indices = set([i for i, x in enumerate(t) if x])
    n = sum(t)
    rv = sum([(lcm // d if d else n_tot*lcm) * sum([t[i]*t[i+d]
                                                    for i in indices
                                                    if i+d in indices
                                                    and i+d < len(t)])
              for d in range(len(t))])
    energy_cache[(t, len_tot, n_tot)] = rv
    return rv


def cumsum(l):
    rv = []
    cs = 0
    for x in l:
        cs += x
        rv.append(cs)
    return rv


def UniformitySum(t):
    n = sum(t)
    l = len(t)

    # We make cumulative and expected (under uniformity) cumulative
    # integers by expressing relative to 1 / (n * l)
    cum = [l * x for x in cumsum(t)]
    cum_uni = cumsum(l*[n])

    return sum([abs(cum[i] - cum_uni[i]) for i in range(len(cum))])



def FindMostUni(n, cat_maxes, cats_left=0, n_left=0, cats_right=0, n_right=0):
    cats_tot = len(cat_maxes) + cats_left + cats_right
    n_tot = n + n_left + n_right
    ex = [(cats_left + 1 + i) * n_tot      # expected under uniformity
          for i in range(len(cat_maxes))] 

    if len(cat_maxes) == 1:
        if n > cat_maxes[0]:
            raise RuntimeError('n (%d) > cat_maxes[0] (%d)' % (n, cat_maxes[0]))
        return [(n,)], abs((n + n_left) * cats_tot - ex[0]) + cats_tot*n_tot*n*(n-1)

    if sum(cat_maxes) == 0 or n == 0:
        return ([len(cat_maxes) * (0,)],
                sum([abs(n_left * cats_tot - x) for x in ex]))
    
    half = len(cat_maxes) // 2
    score_best = (len(cat_maxes) * n_tot * cats_tot   # any real score will be
                  + cats_tot*n_tot*n*(n-1))           # less than this
    dists_best = []
    for n_first in range(0, n + 1):
        n_second = n - n_first
        if n_first > sum(cat_maxes[:half]) or n_second > sum(cat_maxes[half:]):
            # impossible to satisfy
            continue
        dists1, score1 = FindMostUni(n_first, cat_maxes[:half],
                                     cats_left, n_left,
                                     cats_right + len(cat_maxes) - half,
                                     n_right + n_second)
        dists2, score2 = FindMostUni(n_second, cat_maxes[half:],
                                     cats_left + half,
                                     n_left + n_first,
                                     cats_right, n_right)
        score = score1 + score2
        if score <= score_best:
            dists = [d1 + d2 for d1 in dists1 for d2 in dists2]
            if score < score_best:
                score_best = score
                dists_best = []
            dists_best.extend(dists)

    return dists_best, score_best


days_in_mo_common = [None, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
days_in_mo_leap   = [None, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

def SubSample(ids, number, priorities, delim, date_index, aln):
    global verbose
    if verbose:
        print('Subsample %d from %d isolates' % (number, len(ids)))
        print(list(ids)[0])
    order_ecs = [shuffled(ids)]  # Equivalence classes under ordering, in order
    if len(ids) <= number: # could come before shuffling above for speed
        return list(ids)  # new list is defensive
    if priorities != 'N':  # 'N' means "none", i.e., no preferences.
        # Sort by lowest priority first, moving eventually to highest.
        # Relies on the fact that sort is a stabile sort.
        for ch in priorities:
            if ch == 'L':
                #l.sort(key=lambda name: len(aln[name].strip('-')), reverse=True)
                tmp = [GroupBy(l, lambda name: len(aln[name].strip('-')))
                       for l in order_ecs]
                order_ecs = [[d[val] for val in sorted(d, reverse=True)]
                             for d in tmp]
                order_ecs = list(flatten(order_ecs))
            elif ch == 'D':
                tmp = [GroupBy(l, lambda name: DateCompleteness(name,
                                                                delim,
                                                                date_index))
                       for l in order_ecs]
                order_ecs = [[d[val] for val in sorted(d, reverse=True)]
                             for d in tmp]
                order_ecs = list(flatten(order_ecs))
            elif ch == 'M':
                # Calling bool() on completeness
                # produces True iff month is present
                tmp = [GroupBy(l,
                               lambda name: bool(DateCompleteness(name,
                                                                  delim,
                                                                  date_index)))
                       for l in order_ecs]
                order_ecs = [[d[val] for val in sorted(d, reverse=True)]
                             for d in tmp]
                order_ecs = list(flatten(order_ecs))
            elif ch == 'G':
                tmp = [GroupBy(l,
                               lambda name: NumGapsIndivisibleThreeLT5(aln[name]))
                       for l in order_ecs]
                order_ecs = [[d[val] for val in sorted(d)]
                             for d in tmp]
                order_ecs = list(flatten(order_ecs))
            elif ch == 'A':
                tmp = [GroupBy(l, lambda name: sum([aln[name].count(nuc)
                                                    for nuc in 'AGCT-']))
                       for l in order_ecs]
                order_ecs = [[d[val] for val in sorted(d, reverse=True)]
                             for d in tmp]
                order_ecs = list(flatten(order_ecs))
            elif ch == 'N':
                raise ValueError('N cannot occur with additional letters '
                                 'in priority string')
            elif ch == 'U':
                continue
            else:
                raise ValueError('Unrecognized character '
                                 'in priorities string: %s' % ch)
            
    l = list(flatten(order_ecs))
    
    if priorities[0] != 'U':
        return l[:number]

    # Remainder of function deals with date uniformity ('U')
    dc = [DateCompleteness(id, delim, date_index) for id in l]
    with_month = [i for i, x in enumerate(dc) if x > 0]
    #sad = ScoresAllDistrs(min(number, len(with_month)), 12)
    month_counts = 12 * [0]
    month_counts_w_date = 12 * [0]
    for idx in with_month:
        id = l[idx]
        d = id.split(delim)[date_index]
        month = int(d.split('-')[1])
        month_counts[month - 1] += 1
        if dc[idx] == 2:
            month_counts_w_date[month - 1] += 1

    if sum(month_counts) <= number:
        # must use all with month, so cannot do any unifority optimization
        return l[:number]

    for i in itertools.count(1):
        # at most i from each month
        mc_red = [min(c, i) for c in month_counts]
        if sum(mc_red) >= number or mc_red == month_counts:
            break

    if verbose:
        print('Computing distribution across months')

    all_distrs = sorted(AllBestDistrs(min(number, sum(month_counts)),
                                      len(mc_red), mc_red))
    by_score = GroupBy(all_distrs, lambda t: Energy(t))
    cands = by_score[min(by_score)]
    
    # sanity checks
    if len(cands) == 0:
        raise RuntimeError('Zero candidates from FindMostUni')
    if any([sum(cand) != min(number, sum(month_counts)) for cand in cands]):
        raise RuntimeError('Wrong sum of cand: expected %d, got %d'
                           % (min(number, sum(month_counts)), sum(cand)))
    if any([Energy(cand) != min(by_score) for cand in cands]):
        raise RuntimeError('Returned score (%d) '
                           'does not equal calculated score (%s) for %s'
                           % (score, [Energy(cand) for cand in cands],
                              cands))

    
    if len(cands) > 1:
        # pick those requiring the smallest number of incomplete dates
        inc_counts = []
        for cand in cands:
            tmp = [cand[mo] - month_counts_w_date[mo] for mo in range(12)]
            tmp = sum([x for x in tmp if x > 0])
            inc_counts.append(tmp)
        mn = min(inc_counts)
        cands = [cand for i, cand in enumerate(cands) if inc_counts[i] == mn]
        
    yr = int(l[0].split(delim)[date_index].split('-')[0])  # rely all same year
    cands_scores = {}
    for cand in cands:
        chosen = []
        for mo in range(1, 13):
            if cand[mo - 1] == 0:
                continue
            ids_mo = [id for i, id in enumerate(l)
                      if dc[i] == 2
                      and int(id.split(delim)[date_index].split('-')[1]) == mo]
            dates = {id : int(id.split(delim)[date_index].split('-')[2])
                     for id in ids_mo}
            # critical that order within each list is order in l
            ids_by_date = GroupBy(ids_mo, lambda id: dates[id])
            date_counts = {date : len(l) for date, l in ids_by_date.items()}
            if yr % 4:
                days_in_mo = days_in_mo_common[mo]
            else:
                days_in_mo = days_in_mo_leap[mo]
            date_counts = [date_counts.get(i + 1, 0)
                           for i in range(days_in_mo)]
            
            if sum(date_counts) == 0:
                mo_cands = [days_in_mo * (0,)]
            else:
                for i in itertools.count(1):
                    # at most i from each day
                    dc_red = [min(c, i) for c in date_counts]
                    if sum(dc_red) >= cand[mo - 1] or dc_red == date_counts:
                        break
                if verbose:
                    print('Computing distribution within month %d; ' % mo,
                          end='')
                    print('Choosing %d from %d isolates with full dates'
                          % (cand[mo - 1], sum(date_counts)))
                    tick()

                # ress = [MinEnergyHeur(min(cand[mo - 1],
                #                           sum(date_counts)),
                #                       dc_red)
                #         for i in range(10)]
                # print('%d unique' % len(set(ress)))

                # en = [Energy(t) for t in ress]
                # mn = min(en)
                # mo_cands = [t for i, t in enumerate(ress)
                #             if en[i] == mn]

                mo_cands = list(sorted(MinEnergy3(min(cand[mo - 1],
                                                      sum(date_counts)),
                                                  dc_red)))
                                                        
                # sanity checks
                if (len(mo_cands) == 0):
                    raise RuntimeError('0 month candidates')
                if any([sum(mo_cand) != min(cand[mo - 1],
                                            sum(date_counts))
                        for mo_cand in mo_cands]):
                    raise RuntimeError('Wrong number chosen: expected %d; '
                                       'mo_cands = %s'
                                       % (min(cand[mo - 1], sum(date_counts)),
                                          mo_cands))
            

            mo_cand = mo_cands[0]
            to_add = []
            for i, num in enumerate(mo_cand):
                if num:
                    to_add.extend(ids_by_date[i + 1][:num])

            if verbose:
                print('Cat. maxes: %s' % (dc_red,))
                print('Result:     %s' % (mo_cand,))
                print('Minimum energy %d' % Energy(mo_cand))

                tock()
                
            if len(to_add) < cand[mo - 1]:
                nmissing = cand[mo - 1] - len(to_add)
                ids_mo_dateless = [id for idx, id in enumerate(l)
                                   if dc[idx] == 1
                                   and int(id.split(delim)[date_index]
                                   .split('-')[1]) == mo]
                to_add.extend(ids_mo_dateless[:nmissing])
            chosen.extend(to_add)
                        
        # we go with the first candidate
        nmissing = number - len(chosen)
        if nmissing:
            # fill remaining requirement with monthless IDs
            nadded = 0
            for idx, c in enumerate(dc):
                if not c:
                    chosen.append(l[idx])
                    nadded += 1
                    if nadded == nmissing:
                        break

        if len(chosen) < min(number, len(ids)):
            raise RuntimeError('Fewer chosen than expected; ids[0] = %s; '
                               'cand = %s; mo_cand = %s; '
                               'number expected = %d; number chosen = %d'
                               % (ids[0], cand, mo_cand,
                                  min(number, len(ids)), len(chosen)))

        if verbose:
            print()
        
        return chosen
        
    
def IndexToColName(idx):
    '''
    Given a zero-based index, return the corresponding column name
    in a spreadsheet.
    
    ex.:
      0 -> 'A'
      1 -> 'B'
      25 -> 'Z'
      26 -> 'AA'
      27 -> 'AB'
      52 -> 'BA'
    '''
    n = idx + 1  # zero- to one-based
    res = ''
    while n:
        rem = n % 26
        res = chr(64 + rem) + res  # 'A' is 65, 'B' is 66, etc.
        n = n // 26
    return res
        

empty_field_label = '(empty field)'

def WriteEditableCSV(groups, nsamp, fname):
    '''
    Write a file with two counts for each category:
    the total number of isolates in set, and the
    default number to be sampled (min(nsamp, total)).
    Meant to allow editing to express different preferences
    for number to be sampled from each category.

    groups  keys are tuples specifying categories;
            lenghts of values equal total counts
    nsamp   Requested number of samples per category
    fname   Name of file to write
    '''
    import csv
    with open(fname, 'w', newline='') as fid:
        csv_writer = csv.writer(fid)
        nfields = len(list(groups)[0])
        if True: #nfields > 1:
            if nfields > 1:
                nv = (nfields + 1) // 2  # number of fields to run vertically
                nh = nfields - nv  # number of fields to run horizontally
            else:
                # if just one, make it horizontal
                nv = 0
                nh = 1

            fvs = sorted(set(t[:nv] for t in groups))
            fhs = sorted(set(t[nv:] for t in groups))
            row_len = nh + 2 * len(fvs) + 2
            rows_written = 0
            
            # header: nv rows
            rows = [row_len * [''] for i in range(nv)]
            for i, t in enumerate(fvs):
                for j in range(nv):
                    if i == 0 or (t[:j + 1] != fvs[i - 1][:j + 1]):
                        val = fvs[i][j]
                        if val == '':
                            val = empty_field_label
                        rows[j][nh + 2 * i] = val
            csv_writer.writerows(rows)
            rows_written += len(rows)

            # data rows
            for i, fh in enumerate(fhs):
                row = row_len * ['']
                # field values in initial cols
                for j in range(nh):
                    if i == 0 or (fh[:j + 1] != fhs[i - 1][:j + 1]):
                        row[j] = fh[j]
                # counts
                for j, fv in enumerate(fvs):
                    cat = fv + fh  # category (tuple)
                    if cat in groups:
                        tot_count = len(groups[cat])
                        samp_count = min([nsamp, tot_count])  # sample <= tot
                        row[nh + 2 * j] = str(tot_count)
                        row[nh + 2 * j + 1] = str(samp_count)
                
                # Total samples formula
                if i == 0:
                    row[-1] = 'Total Samples'
                    row.append('Force Inclusion of:')
                if i == 1:
                    sum_arg = ','.join(['%s%d:%s%d'
                                        % (IndexToColName(col_idx), nv + 1,
                                           IndexToColName(col_idx), nv + len(fhs))
                                        for col_idx in
                                        range(nh + 1, nh + 1 + 2 * max([len(fvs), 1]), 2)])
                    sum_formula = 'sum(%s)' % sum_arg
                    row[-1] = '=' + sum_formula
                csv_writer.writerow(row)
                rows_written += 1

            
def ReadEditableCSV(fname):
    '''
    Read sample count info for categories from CSV file
    produced by WriteEditableCSV and possibly edited.
    
    Returns requested sample counts by category.
    '''
    # Read rows
    import csv
    with open(fname, newline='') as fid:
        rdr = csv.reader(fid)
        rows = [r for r in rdr]

    # First data row is first to have nonempty col 1
    nv = [i for i in range(len(rows)) if rows[i][0] != ''][0]  # vert fields

    hdr = rows[:nv]
    data = rows[nv:]

    if len(hdr) == 0:
        # Given rules for writing, this should only happen when
        # there is a single field, with horizontal orientation,
        # meaning (confusingly?) that each row has
        # category value in first col, total count in second,
        # and sampled count in third.
        res = {(row[0],) : int(row[2])
               for row in data
               if row[2].strip() and int(row[2]) != 0}
        return res

    # Last row of header tells us number of vert cats and number of hdr fields
    nvcats = len([s for s in hdr[-1] if s != ''])  # num vert cats
    nh = [i for i in range(len(hdr[-1])) if hdr[-1][i] != ''][0]
        
    nhcats = len(data)
        
    # Fill in header blanks by repeating values to their right into blanks
    for row in hdr:
        for i in range(nh, nh + 2 * nvcats):
            if row[i] == '':
                row[i] = row[i - 1]

    # Any originally empty fields have been given a nonempty label;
    # replace them with empty string.
    for row in hdr:
        if empty_field_label in row:
            for i in range(len(row)):
                if row[i] == empty_field_label:
                    row[i] = ''
                
    # Vertical cats, in left to right order
    tp = [x for x in zip(*hdr)]
    vcats = [tp[i] for i in range(nh + 1, nh + 2 * nvcats + 1, 2)]

    # Get the sample counts
    counts = {}
    hcat_prev = None
    for rnum, row in enumerate(data):
        hcat = row[:nh]
        for i in range(len(hcat)):
            if hcat[i] == '':
                hcat[i] = hcat_prev[i]
        
        # Any originally empty fields have been given a nonempty label;
        # replace them with empty string.
        if empty_field_label in hcat:
            for i in range(len(hcat)):
                if hcat[i] == empty_field_label:
                    hcat[i] = ''

        # Get counts for categories with nonzero counts
        hcat = tuple(hcat)
        hcat_prev = hcat
        for i in range(nvcats):
            cnum = nh + 2 * i + 1
            if row[cnum].strip():
                count = int(row[cnum])
                if count:
                    counts[vcats[i] + hcat] = count

    # Forced includes
    cnum = nh + 2 * nvcats + 2
    row = data[1]
    if len(row) >= cnum + 1:
        # Allow accessions or identifiers to be separated by ',' or whitespace
        forced_includes = row[cnum].replace(',', ' ').split()
    
    return counts, forced_includes


def LaunchFile(fname):
    '''
    Open a file with appropriate (hopefully) application.
    Return 0 on success, return (*not* raise) exception on failure.

    Uses OS-dependent mechanisms.
    For other than Windows and MacOS (e.g., linux), best we can do
    is through web browser.
    '''
    try:
        plat = sys.platform
        if plat.startswith('win32'):
            import subprocess
            # Extra parameter ('') necessary if
            # file name has spaces
            subprocess.call(['start', '', fname],
                            shell=True)
        elif plat.startswith('darwin'):
            import subprocess
            subprocess.call(['open', fname])
        else:
            webbrowser.open(fname)
        return 0  # 0 on success
    except Exception as e:
        # Do not allow the program to fail because of failed launch
        return e  # *return* exception, which tests as True, on failure

                    
                    

################### Find largest min distance among removed #################

# Necessary do-nothing functions

def tic():
    pass
    
    
def toc(s=None):
    pass


def profile(f):
    return f
    

# The class that does the computations
    
class MismatchCount(object):
    '''
    The purpose of this code is to determine the largest minimum distance
    of the sequences in one set to any of the sequences in another set.
    With the first set referred to as the queries and the second as the
    subjects, this is max for q in Q of (min for s in S of d(q, s)).
    This minimum, along with the corresponding query indices, is returned
    by MaxMinDist().  For simplest use of the class, the only other method 
    that might be called by user is convenience function EndGapsToN().

    Computing all pairwise dists for large sets would be slow, but is
    unnecessary.  Comparing some queries to all subjects gives a lower
    bound on max.  Comparing all queries to some subjects gives an upper
    bound for each min.  Queries for which the latter is smaller than the
    former can then be eliminated from consideration.  Additional comparisons
    and eliminations follow.

    Additional speed-ups:

    1. Calculate min differences on partitions of the sequence positions,
       after uniquifying (MMCountPieces).
       I.e., use seq[0::10], seq[1::10] ... seq[9::10], and uniquify before
       calculating distances.  These smaller sequences are more commonly
       identical between isolates, and number to compare can go down by_min
       a large factor.

    2. All query against all subject comparison of seq[0::10]
       (with uniquification) is fast and immediately provides good candidates
       for large minimum distance.
    '''

    nucs = 'AGCTRYSWMKBHDVN'
    
    bits = [0b0001, 0b0010, 0b0100, 0b1000,                  # AGCT
            0b0011, 0b1100, 0b0110, 0b1001, 0b0101, 0b1010,  # RYSWMK
            0b1110, 0b1101, 0b1011, 0b0111,                  # BHDV
            0b1111]                                          # N

    nucs_match = {}
    for j, nuc2 in enumerate(nucs):
        for i, nuc1 in enumerate(nucs):
            nucs_match[(nuc1, nuc2)] = bool(bits[i] & bits[j])   

    nucs_match['-', '-'] = True
    for nuc in nucs:
        nucs_match[nuc, '-'] = False
        nucs_match['-', nuc] = False

    verbose = False
    
    
    @staticmethod
    def EndGapsToN(seq):
        orig_len = len(seq)
        s1 = seq.lstrip('-')
        lgaps = orig_len - len(s1)
        seq = s1.rstrip('-')
        rgaps = len(s1) - len(seq)
        return lgaps * 'N' + seq + rgaps * 'N'

    
    @staticmethod
    def MakeColDescs(seqs):
        '''
        For each column of alignment, make a dict mapping characters to
        rows in which they occur.
        '''
        cols = zip(*seqs)
        groups = [GroupBy(range(len(seqs)), lambda i: col[i])
                  for col in cols]
        return groups


    @classmethod
    def MMCount(cls, queries, subjects, scol_descs=None):
        nucs_match = cls.nucs_match
        if scol_descs is None:
            nsubjects = len(subjects)
            scol_descs = cls.MakeColDescs(subjects)
        else:
            # do not use arg 'subjects', so that it can be, e.g., None 
            nsubjects = sum([len(l) for l in scol_descs[0].values()])
            
        qcol_descs = cls.MakeColDescs(queries)
        all_counts = [nsubjects * [0] for i in range(len(queries))]
        all_common_counts = [0 for i in range(len(queries))]
        for pos, qcol_desc in enumerate(qcol_descs):
            non_matching = {qnuc : [snuc for snuc in scol_descs[pos]
                                    if not nucs_match[qnuc, snuc]]
                            for qnuc in qcol_descs[pos]}
            for qnuc in qcol_descs[pos]:
                num_non_matching = sum(len(scol_descs[pos][snuc])
                                       for snuc in non_matching[qnuc])
                if num_non_matching <= nsubjects - num_non_matching:
                    # Straightforward: increment appropriate mismatch counts
                    snums = list(flatten(scol_descs[pos][snuc]
                                         for snuc in non_matching[qnuc]))
                    for qnum in qcol_descs[pos][qnuc]:
                        #all_counts_qnum = all_counts[qnum]  # just for speed
                        for snum in snums:
                            all_counts[qnum][snum] += 1
                else:
                    # Time-saving trick:
                    # increment overall count, decrement matchers
                    snums = list(flatten(scol_descs[pos][snuc]
                                                for snuc in scol_descs[pos]
                                                if snuc not in non_matching[qnuc]))
                    for qnum in qcol_descs[pos][qnuc]:
                        all_common_counts[qnum] += 1
                        for snum in snums:
                            all_counts[qnum][snum] -= 1

        # Finish the time-saving trick
        for qnum in range(len(queries)):
            for snum in range(nsubjects):
                all_counts[qnum][snum] += all_common_counts[qnum]
                
        return all_counts


    @classmethod
    def MMCountPieces(cls, queries, subjects, pieces, col_descs=None):
        '''
        Computes min distance over subjects for each query,
        analyzing subsets of positions separately for speed.
        
        This is fast for typical sequence sets because number of
        unique values of, e.g., seq[i::10], for any 0 <= i < 10, is
        much smaller than the number of sequences for a large set.
        
        The function returns only the minimimum for each query.
        Computing the distances for every query/subject pair would
        be too time-consuming.
        
        queries, subjects  are iterables of sequences
        pieces             is the number of subsets of positions to use
        col_descs          optional; results of MakeColDescs for each
                           position subset of subjects.  Can speed things
                           up if caller already has it.
        '''
        
        tic()
        qpieces = [GroupBy(range(len(queries)),
                                  lambda i: queries[i][start::pieces])
                   for start in range(pieces)]
        spieces = [GroupBy(range(len(subjects)),
                                  lambda i: subjects[i][start::pieces])
                   for start in range(pieces)]
        if col_descs is None:
            res = [cls.MMCount(qpieces[i].keys(), spieces[i].keys())
                   for i in range(pieces)]
        else:
            res = [cls.MMCount(qpieces[i].keys(), None, col_descs[i])
                   for i in range(pieces)]
        toc('counts of pieces')
        qmap = [{seq : i for i, seq in enumerate(qp)} for qp in qpieces]
        smap = [{seq : i for i, seq in enumerate(sp)} for sp in spieces]
        all_sindices = [[smap[p][s[p::pieces]] for p in range(pieces)]
                        for s in subjects]
        mins = []
        for i, q in enumerate(queries):
            qindices = [qmap[p][q[p::pieces]] for p in range(pieces)]
            mn = len(q)  # dist can't be any larger than this
            # Find the min for this query over all subjects
            for j, s in enumerate(subjects):
                sindices = all_sindices[j]
                sm = 0  # sum of dists for pieces
                for p in range(pieces):
                    sm += res[p][qindices[p]][sindices[p]]
                    if sm >= mn:
                        # Total cannot be smaller than min so far
                        break
                else:
                    mn = sm
            mins.append(mn)
        return mins

        
    class TooManyItError(RuntimeError):
        pass


    @profile
    @classmethod
    def MaxMinDist(cls, queries, subjects, end_gaps_to_n):
        '''
        Find the max over all q of the min over all s of the distance.
        '''
        skip = max([len(queries) // 65, 1])
        skip2 = max([skip // 2, 1])

        cands = list(range(len(queries)))  # All are initially candidates

        queries = [s.upper() for s in queries]
        subjects = [s.upper() for s in subjects]
        
        if end_gaps_to_n:
            queries = [cls.EndGapsToN(s) for s in queries]
            subjects = [cls.EndGapsToN(s) for s in subjects]            
        
        q_part = 0
        s_part = 0

        cur_max = 0
        cur_max_queries = []

        # Upper bound for min dist for each query
        ubs = [len(queries[0]) for q in queries]

        tic()
        pieces = 10
        spieces = [GroupBy(range(len(subjects)),
                           lambda i: subjects[i][start::pieces])
                   for start in range(pieces)]
        scol_descs10 = [cls.MakeColDescs(x) for x in spieces]
        searched_against_all = set()
        toc('Col descs for subjects in 10 pieces')
        for itnum in itertools.count():
            if cls.verbose:
                print('%d candidates' % len(cands))

            if len(cands) <= 20 and itnum > 0:
                # With a small number of candidates remaining,
                # cut to the chase.
                tic()
                mins = cls.MMCountPieces([queries[i] for i in cands], subjects,
                                         10, scol_descs10)
                toc('MMCount 3')
                max_min = max(mins)
                return max_min, [qnum for i, qnum in enumerate(cands)
                                 if mins[i] == max_min]
            
            # A fraction of query candidates against all subjects
            # to get lower bound on max of min
            if itnum == 0:
                # Comparison of all queries to all subjects
                # using just a fraction of positions.
                # This is a fast way to find good candidates (among queries)
                # for high minimum distance (to subjects).
                q10 = [q[::10] for q in queries]
                s10 = [s[::10] for s in subjects]
                qset = list(set(q10))
                sset = list(set(s10))
                res = cls.MMCount(qset, sset)
                mins = [min(r) for r in res]
                min_for_seq = {qset[i] : mins[i] for i in range(len(qset))}
                by_min = sorted(range(len(q10)),
                                key=lambda i: min_for_seq[q10[i]],
                                reverse=True)
                query_subset = by_min[:20]  # Likely to have high min. distance
            elif itnum % 2 == 0:
                query_subset = [i for i in cands if i % skip2 == q_part]
                q_part += 1
                cands_not_searched = [i for i in cands
                                      if i not in searched_against_all]
                cs = sorted(cands_not_searched,
                            key=lambda c: ubs[c],
                            reverse=True)
                query_subset += cs[:10]

                # Remove any that have already been searched against all
                query_subset = [qnum for qnum in query_subset            
                                if qnum not in searched_against_all] 
            
            if query_subset and itnum % 2 == 0:
                tic()
                mins = cls.MMCountPieces([queries[i] for i in query_subset],
                                         subjects,
                                         10, scol_descs10)
                toc('MMCount 1')
                searched_against_all.update(query_subset)
                for i, idx in enumerate(query_subset):
                    ubs[idx] = mins[i]
                max_min = max(mins)
                if max_min >= cur_max:
                    if cls.verbose:
                        print('max %d' % max_min)
                    if max_min > cur_max:
                        cur_max = max_min
                        cur_max_queries = [qnum for i, qnum
                                           in enumerate(query_subset)
                                           if mins[i] == max_min]
                        cands = [i for i in cands if ubs[i] >= cur_max]
                    else:
                        # equal to cur_max; append indices
                        cur_max_queries += [qnum for i, qnum
                                            in enumerate(query_subset)
                                            if mins[i] == max_min]
                    if len(set(cur_max_queries)) < len(cur_max_queries):
                        raise RuntimeError('Duplicate values in '
                                           'cur_max_queries: %s'
                                           % cur_max_queries)
            

            # All query candidates against a fraction of subjects
            # to get upper bound on min for each cand
            subject_subset = list(range(s_part, len(subjects), skip))
            subject_subset_seqs = [subjects[i] for i in subject_subset]
            s_part += 1
            tic()
            mins = cls.MMCountPieces([queries[i] for i in cands],
                                      subject_subset_seqs, 10)
            toc('MMCount 2')
            
            to_rm = set()  # candidates whose ub now falls below max
            for i, qnum in enumerate(cands):
                if mins[i] < ubs[qnum]:
                    ubs[qnum] = mins[i]
                if ubs[qnum] < cur_max:
                    to_rm.add(qnum)
            if to_rm:
                cands = [i for i in cands if i not in to_rm]

                
            # Are we done?
            # One way:
            if len(cands) == len(cur_max_queries):
                if set(cands) != set(cur_max_queries):
                    raise RuntimeError('set(cands) != set(cur_max_queries)',
                                       cands, cur_max_queries)
                return cur_max, cur_max_queries
            # Related sanity check:
            if len(cands) < len(cur_max_queries):
                raise RuntimeError('len(cands) < len(cur_max_queries)',
                                   cands, cur_max_queries)
            # Another way to be done:
            if s_part == skip or searched_against_all.issuperset(cands):
                # Upper bounds for all remaining candidates are actual mins
                mx = max(ubs[i] for i in cands)
                arg_max = [i for i in cands if ubs[i] == mx]
                return mx, arg_max
                
##############################################################################


#################### GUI ####################

class GuiError(RuntimeError):
   pass

class Gui(object):
    initialized = False
    dir = ''
    
    @staticmethod
    def InitGui():
        if Gui.initialized:
            return
        try:
            import wx
        except ImportError as e:
            raise GuiError('wxPython absent or incorrectly installed; '
                           'try "pip install wxPython"')
        Gui.app = wx.App()
        Gui.frame = wx.Frame()
        Gui.initialized = True
        
    @staticmethod
    def FileDialog(mode, prompt, type_text):
        Gui.InitGui()
        import wx
        if mode == 'r':
            flags = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
        elif mode == 'w':
            flags = wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT
        else:
            raise ValueError("'%s' is not a valid mode; "
                             "should be 'r' for input or 'w' for output"
                            % mode)
        fileDialog = wx.FileDialog(Gui.frame, prompt, Gui.dir, "", type_text, flags)
        res = fileDialog.ShowModal()
        if res == wx.ID_CANCEL:
            return None
        path = fileDialog.GetPath()
        if path.strip():
            Gui.dir = fileDialog.GetDirectory()
        fileDialog.Destroy()
        return path
            
#################### End GUI ################ 


start_time = time.localtime()

try:
    optlist, args = getopt.getopt(sys.argv[1:], 's:vc:h',
                                  ['seed=',
                                   'verbose',
                                   'comment=',
                                   'help'])

    optdict = dict(optlist)
    if len(args) not in [0, 1, 2]:
        raise getopt.GetoptError('Too many command-line parameters; '
                                 'maximum is 2')
    verbose = '-v' in optdict or '--verbose' in optdict
    if '-s' in optdict or '--seed' in optdict:
        if '-s' in optdict and '--seed' in optdict:
            raise getopt.GetoptError('Illegal to give both -s and --seed')
        if '-s' in optdict:
            seed = int(optdict['-s'])
        else:
            seed = int(optdict['--seed'])
    else:
        seed = None

    if '-c' in optdict or '--comment' in optdict:
        if '-c' in optdict and '--comment' in optdict:
            raise getopt.GetoptError('Illegal to give both -c and --comment')
        if '-c' in optdict:
            comment = optdict['-c']
        else:
            comment = (optdict['--comment'])
    else:
        comment = None
    if '-h' in optdict or '--help' in optdict:
        raise getopt.GetoptError('')
except getopt.GetoptError as e:
    sys.stderr.write(e.msg)

    sys.stderr.write('''
Usage: %s [OPTIONS] [input_fname [output_fname]]

  OPTIONS
    -v,--verbose   Print a lot of information
    -s,--seed      With an integer, sets random number seed
    -c,--comment   With a string, sets comment in top of output fasta
    -h,--help      This help
''' % sys.argv[0])
    Quit(1)

print()
print('** Modified program for COVID: changed meanings of A and G criteria **')
print()

if seed is not None:
    random.seed(seed)

summaries = []

for itnum in itertools.count():
    if itnum == 0 and len(args) > 0:
        in_fname = args[1]
    elif (itnum > 0 and
          GetYN('Use same input file (%s)?  y(es)/N(o) (RETURN for No): '
                 % in_fname, False)):
        # User wants to use same input file as in last iteration;
        # leave value as it is
        pass
    else:
        in_fname = GetInputFname('Enter FASTA-format file for input: ',
                                 'fasta files (*.fa; *.fas; *.fasta)'
                                 '|*.fa;*.fas;*.fasta|all files (*.*)|*.*')

    if itnum == 0 and len(args) > 1:
        out_fname = args[1]
    else:
        out_fname = GetOutputFname('Enter name for FASTA-format output file: ',
                                   'fasta files (*.fa)|*.fa'
                                   '|all files (*.*)|*.*')

    st_input = os.stat(in_fname)                          
    seqs, names = fasta_read(in_fname)

    aln = {name : seqs[i] for i, name in enumerate(names)}

    default_delim = '|'
    delim = GetSingleCharacter("Enter delimiter, or return for '%s': "
                                % default_delim,
                                  default_delim)

    ShowExampleTagVals(names, delim)
    print('Note: only the year will be used if the date field is included')
    print()

    date_index = FindDateField(aln.keys(), delim)

    using_ss_input = False
    fields = GetFieldChoices(len(names[0].split(delim)))

    if fields == 'S':
        # Request to load spreadsheet
        using_ss_input = True
        csv_input_fname = GetInputFname('Enter CSV-format file for input: ',
                                        'CSV files (*.csv)|*.csv'
                                        '|all files (*.*)|*.*')
        target_counts, forced_includes = ReadEditableCSV(csv_input_fname)

        # Figure out which fields these correspond to
        values = [set(x) for x in zip(*[n.split(delim) for n in names])]
        values[date_index] = set(s[:4] for s in values[date_index])
        ss_values = [set(x) for x in zip(*target_counts.keys())]
        fields = []
        for ss in ss_values:
            fnum, = [i for i, s in enumerate(values) if s.issuperset(ss)]
            fields.append(fnum)
            
        
    # Group by chosen fields
    groups = GroupBy(aln.keys(),
                     lambda s: FieldTuple(s, delim, fields, date_index))

    print()
    print('%d categories' % len(groups))
    print('Mean category size: %.2f'
          % (sum(len(l) for l in groups.values()) / float(len(groups))))
    print('Minimum: %d' % min(len(l) for l in groups.values()))
    print('Maximum: %d' % max(len(l) for l in groups.values()))
    print()

    if len(fields) == 1 or len(fields) == 2:
        do_ss = GetYN('Save a spreadsheet showing '
                      'category sizes? y(es)/N(o) (RETURN for No): ', False)
        if do_ss:
        
            csv_out_fname = GetOutputFname('Enter name for output spreadsheet (CSV) file: ',
                                           'CSV files (*.csv)|*.csv')
            import csv
            with open(csv_out_fname, 'w', newline='') as fid:
                csv_writer = csv.writer(fid)
                if len(fields) == 2:
                    f1s = sorted(set(t[0] for t in groups))
                    f2s = sorted(set(t[1] for t in groups))
                    csv_writer.writerow([''] + f2s + ['', 'Total'])  # Col names
                    for f1 in f1s:
                        row = ([f1]
                               + [len(groups[f1, f2]) if (f1, f2) in groups
                                  else ''  # or 0 here
                                 for f2 in f2s]
                               + ['']
                               # Row sum
                               + [sum(len(l) for t, l in groups.items()
                                  if t[0] == f1)])
                        csv_writer.writerow(row)
                    csv_writer.writerow([])
                    # Column sums
                    csv_writer.writerow(['Total'] +                      
                                        [sum([len(l) for t, l in groups.items()
                                         if t[1] == f2]) for f2 in f2s] +
                                         [''] +
                                         [sum(len(l) for l in groups.values())]
                                         )
                else:  # len(fields) must equal 1
                    f1s = sorted(set(t[0] for t in groups))
                    for f1 in f1s:
                        row = [f1, len(groups[f1,])]
                        csv_writer.writerow(row)
                    csv_writer.writerow([])
                    csv_writer.writerow(['Total',
                                         sum(len(l) for l in groups.values())])
            st_csv = os.stat(csv_out_fname)
            print()
            print('Wrote CSV file %s' % csv_out_fname)
            print()
            if GetYN('Launch CSV file %s now? y(es)/N(o) (RETURN for No): '
                     % csv_out_fname, False):
                try:
                    plat = sys.platform
                    if plat.startswith('win32'):
                        import subprocess
                        # Extra parameter ('') necessary if
                        # file name has spaces
                        subprocess.call(['start', '', csv_out_fname],
                                        shell=True)
                    elif plat.startswith('darwin'):
                        import subprocess
                        subprocess.call(['open', csv_out_fname])
                    else:
                        webbrowser.open(csv_out_fname)
                    print()
                    print('Launched CSV file')
                    print()
                except Exception as e:  # Do not allow program to fail
                    print()
                    print('Problem launching CSV file: %s' % e)
                    print()
            else:
                print()
        else:
            # No ss
            print()


    if not using_ss_input:
        # Show user number of samples that will result from various
        # choices of samples per category.
        cat_sizes = [len(l) for l in groups.values()]
        max_cat_size = max(cat_sizes)
        print('Total samples for some choices of samples per category:')
        print()
        pairs = []
        for num_samps in itertools.chain(range(1, 10),
                                         range(10, 100, 10),
                                         range(100, 1001, 100)):
            if num_samps > max_cat_size:
                num_samps = max_cat_size
            total = sum(min([num_samps, sz]) for sz in cat_sizes)
            pairs.append('%4d  %-10d' % (num_samps, total))
            if num_samps == max_cat_size:
                pairs[-1] = pairs[-1].rstrip() + ' (all)'
                break

        # Make two columns
        if len(pairs) % 2:
            pairs.append('')
        nlines = len(pairs) // 2
        for i in range(nlines):
            print(pairs[i] + 10 * ' ' + pairs[i + nlines].rstrip())
            
        print()
        
        num_samps = GetPosInt('Number to sample from each category: ')

        total = sum(min([num_samps, sz]) for sz in cat_sizes)
        print('%d sequences will be sampled' % total)
        print()
        
        do_ss2 = GetYN('Save a spreadsheet showing '
                       'category sizes and number sampled '
                       'for editing and reloading? '
                       'y(es)/N(o) (RETURN for No): ', False)
        if do_ss2:
            csv_out_fname2 = GetOutputFname('Enter name for output spreadsheet '
                                            '(CSV) file: ',
                                            'CSV files (*.csv)|*.csv')
            WriteEditableCSV(groups, num_samps, csv_out_fname2)
            print()
            print('Wrote CSV file %s' % csv_out_fname2)
            print()
            if GetYN('Launch CSV file %s now? y(es)/N(o) (RETURN for No): '
                     % csv_out_fname2, False):
                status = LaunchFile(csv_out_fname2)
                print()
                if not status:
                    print('Launched CSV file')
                else:
                    # Print Exception object
                    print('Problem launching CSV file: %s' % status)  
                print()
        else:
            # No ss2
            print()
        

    while True:
        priorities = GetPrefsWithinCat()
        if priorities[0] == 'U' and date_index not in fields:
            print('\n** U cannot be used unless date field is part of category **\n')
        else:
            break

    if not using_ss_input:
        target_counts = {key : min([num_samps, len(value)])
                         for key, value in groups.items()}
        forced_includes = []

    # Deal with forced includes: group by category, remove from category,
    # and adjust target count downward.
    forced_includes_for_cat = {}
    if forced_includes:
        print()
        print('Forcing inclusion of %d sequences' % len(forced_includes))
        print()
        acc_to_id = {id.split(delim)[0] : id for id in aln}
        id_to_cat = {}
        for cat, ids in groups.items():
            for id in ids:
                id_to_cat[id] = cat
        for fi in forced_includes:
            if len(fi.split(delim)) == 1:
                try:
                    id = acc_to_id[fi]
                except KeyError as e:
                    raise ValueError('Accession %s (forced include) not found'
                                     % fi)
            else:
                id = fi

            try:
                cat = id_to_cat[id]
            except KeyError as e:
                raise ValueError('ID %s (forced include) not found'
                                 % id)

            if cat not in forced_includes_for_cat:
                forced_includes_for_cat[cat] = []
            forced_includes_for_cat[cat].append(id)

        for cat, ids in forced_includes_for_cat.items():
            if len(ids) > target_counts[cat]:
                print('WARNING: Only %d sequences requested for %s,'
                      % (target_counts[cat], cat))
                print('         but inclusion of %d forced:' % len(ids))
                print('         %s' % ' '.join([id.split('|')[0] for id in ids]))
                print('         All %d will be included, increasing total by %d.'
                      % (len(ids), len(ids) - target_counts[cat]))
                print()
                target_counts[cat] = 0
            else:
                target_counts[cat] -= len(ids)
                
            for id in ids:
                groups[cat].remove(id)
            
    sampled = {key : (SubSample(value, target_counts[key], priorities,
                                delim, date_index, aln)
                      + forced_includes_for_cat.get(key, []))
               for key, value in groups.items()
               if key in target_counts}

    sampled_ids = list(flatten(sampled.values()))

    # A comment to add to first line of fasta
    if comment is None or comment.startswith('+'):
        if using_ss_input:
            num_samps_str = 'from %s' % csv_input_fname
        else:
            num_samps_str = str(num_samps)

        try:
            login = os.getlogin() + ' '  # os.getlogin() fails on some linux
        except Exception:
            login = ''

        tmp = ('%s %s fields: %s prefs: %s samps_per_cat: %s (%s %s)'
               % (os.path.basename(sys.argv[0]), in_fname,
                  fields, priorities,
                  num_samps_str, login, StrFTime('%d %b %Y %X %Z')))
        if comment is not None:
            comment = tmp + '; ' + comment[1:]
        else:
            comment = tmp
        
    writefasta([aln[id] for id in sampled_ids], sampled_ids,
               out_fname, [comment])
    
    # Description of actions taken
    in_mtime = datetime.datetime.fromtimestamp(st_input.st_mtime)
    digest = hashlib.md5(open(in_fname, 'rb').read()).hexdigest()
    line0 = ('Input file: %s\n(modified %s; size %d; MD5 %s)'
             % (in_fname, in_mtime.strftime('%d %b %Y %X %Z'),
                st_input.st_size, digest)
            )
    if seed is not None:
        line_seed = 'Random seed set to %d' % seed
    if not using_ss_input:
        line1 = ('Chose %d (or all) from each category '
                 'using fields %s and preferences %s'
                 % (num_samps, fields, priorities))
    else:
        line1 = ('Chose according to %s from each category '
                 'using fields %s and preferences %s'
                 % (csv_input_fname, fields, priorities))

    line2 = '%d sequences written to %s' % (len(sampled_ids), out_fname)
    if seed is None:
        summaries.append([line0, line1, line2])
    else:
        summaries.append([line0, line_seed, line1, line2])
    print()
    print(line1)
    if seed is not None:
        print(line_seed)
    print(line2)

    # Report largest min distance of non-sampled isolate to sampled isolates
    ssampled_ids = set(sampled_ids)
    non_sampled_ids = [id for id in names if id not in ssampled_ids]
    try:
        max_dist, indices = MismatchCount.MaxMinDist(
                                                [aln[id]
                                                 for id in non_sampled_ids],
                                                [aln[id]
                                                 for id in sampled_ids],
                                                True)
    except Exception as e:
        print()
        print('Problem calculating largest min distance:')
        sys.excepthook(type(e), e, e.__traceback__)
        print()
        print('This will not affect anything else.')
        print()
    else:
        print()
        print('Largest distance from non-sampled sequence '
              'to nearest sampled sequence:')
        print()
        print('          %d mismatches' % max_dist)
        print()
        print('%d sequence%s with this minimum distance:'
              % (len(indices), 's' if len(indices) > 1 else ''))
        # Print ID(s) for seqs with max min dist
        for idx in indices:
            print(non_sampled_ids[idx])
        print()
    
    print()
    if not GetYN('Create another subsample? y(es)/N(o) (RETURN for No): ',
                 False):
        break  # Done with iterations at user's request

end_time = time.localtime()
 
if len(summaries) > 0:
    try:
        login = os.getlogin()  # os.getlogin() fails on some linux
    except Exception:
        login = 'unknown'
    print()
    print()
    print('Summary of session by %s %s:'
          % (login, TimeRange(start_time, end_time)))
    print()
    for l in summaries:
        for line in l:
            print(line)
        print()

Quit(0)


      

