import numpy as np

def forward(A, B, pi, O):
    """
    Calculates the probability of an observation sequence O given the model(A, B, pi).
    :param A: state transition probabilities (NxN)
    :param B: observation probabilites (NxM)
    :param pi: initial state probabilities (N)
    :param O: sequence of observations(T) where observations are just indices for the columns of B (0-indexed)
        N is the number of states,
        M is the number of possible observations, and
        T is the sequence length.
    :return: The probability of the observation sequence and the calculated alphas in the Trellis diagram with shape
             (N, T) which should be a numpy array.
    """
    n, m, t, result, init_state = B.shape[0], B.shape[1], len(O), 0, O[0]
    trellis = np.zeros((n,t))

    for i in range(n):
        trellis[i][0] = pi[i] * B[i][init_state]

    for obs in range(1,t):
        for state in range(n):
            for prev_state in range(n):
                trellis[state][obs] += trellis[prev_state][obs-1] * A[prev_state][state] * B[state][O[obs]]

    for i in range(n):
        result += trellis[i][t-1]

    return result, trellis

def viterbi(A, B, pi, O):
    """
    Calculates the most likely state sequence given model(A, B, pi) and observation sequence.
    :param A: state transition probabilities (NxN)
    :param B: observation probabilites (NxM)
    :param pi: initial state probabilities(N)
    :param O: sequence of observations(T) where observations are just indices for the columns of B (0-indexed)
        N is the number of states,
        M is the number of possible observations, and
        T is the sequence length.
    :return: The most likely state sequence with shape (T,) and the calculated deltas in the Trellis diagram with shape
             (N, T). They should be numpy arrays.
    """
    n, m, t, init_state, max_vit = B.shape[0], B.shape[1], len(O), O[0], 0
    trellis = np.zeros((n,t))
    result = np.zeros((t,),dtype=int)

    for i in range(n):
        trellis[i][0] = pi[i] * B[i][init_state]

    for obs in range(1,t):
        for state in range(n):
            pmax = 0
            for prev_state in range(n):
                p = trellis[prev_state][obs-1] * A[prev_state][state] * B[state][O[obs]]
                if p>pmax:
                    pmax = p
                    result[obs-1] = prev_state                
            trellis[state][obs] = pmax

    for i in range(n):
        if trellis[i][t-1] > trellis[i][max_vit]:
            max_vit = i
    result[-1] = max_vit

    return result, trellis
