"""Immutable definition-level verifier. Search programs must not edit this file.

Uses literal retained-position combinations, independently of bit-level search.
Validates generalized n/t instances and cross-checks the frozen arena verifier at
n=16,t=2. No score from an infeasible state is a code-size result.
"""
import hashlib
import itertools
import json
from pathlib import Path


def verify(words, n=16, deletions=2):
    if not isinstance(words, list) or not words or len(set(words)) != len(words):
        raise ValueError('Nonempty distinct words required')
    owner = {}
    for word in words:
        if not isinstance(word, str) or len(word) != n or set(word) - {'0', '1'}:
            raise ValueError('Malformed binary word')
        received = {''.join(word[i] for i in keep)
                    for keep in itertools.combinations(range(n), n-deletions)}
        for value in received:
            if value in owner:
                raise ValueError(f'Collision: {owner[value]} and {word} -> {value}')
            owner[value] = word
    return {'size': len(words), 'n': n, 'deletions': deletions,
            'used_received_words': len(owner),
            'canonical_sha256': hashlib.sha256(('\n'.join(sorted(words))+'\n').encode()).hexdigest()}


def verify_file(path, n=16, deletions=2):
    data=json.loads(Path(path).read_text())
    result=verify(data['words'],n,deletions)
    if n==16 and deletions==2:
        import sys
        sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
        from verify import official
        result['arena_score']=official('two-deletion-code-16',data)
        assert result['arena_score']==result['size']
    return result


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--n',type=int,default=16)
    a=p.parse_args();print(json.dumps(verify_file(a.candidate,a.n),indent=2))
