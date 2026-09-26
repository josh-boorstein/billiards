"""Formal development of a word: the translation of a return word as a function of alpha.

Identify the plane with C.  The reflections in the three sides of the triangle of
`triangle.py` are the affine maps

    r_L2(w) = conj(w),     r_H(w) = e^{2 i alpha} conj(w),     r_L1(w) = -conj(w) + 2L,

with L = cot(alpha).  Composing them along a word gives a map
    w -> e^{i psi(a)} conj^k(w) + L * sum_{(h,q)} b[(h,q)] e^{i psi(h,q)},
    psi(h,q) = 2 h alpha + q pi,
with integer data (a, k, b) that do not depend on alpha.  The composition is done
formally on that data -- no angle is substituted and nothing is reduced modulo any
cyclotomic polynomial -- so the result is valid at every alpha at once.

A return word whose linear part is the identity (k = 0, h = 0, q even) develops to a
translation w -> w + T.  For the return to L1 in direction theta the offset along L1 is

    c = -Im(T e^{-i theta}) / cos(theta),

so the difference of two such words' offsets is  -(L / cos theta) * sum db sin(psi - theta),
an integer combination of sines: `sine_sum` returns it with theta = u alpha + v pi.
"""
from mpmath import mp, mpc, exp, im, cos, sin, pi


def _neg(hq):
    return (-hq[0], -hq[1])


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1])


REFL = {'2': ((0, 0), 1, {}),                   # conj
        'h': ((1, 0), 1, {}),                   # e^{2 i alpha} conj
        '1': ((0, 1), 1, {(0, 0): 2})}          # -conj + 2L   (e^{i pi} = -1)


def compose(G, R):
    """(G o R)(w) = aG conj^kG(aR conj^kR(w) + bR) + bG, on the formal data."""
    aG, kG, bG = G
    aR, kR, bR = R
    f = _neg if kG % 2 else (lambda hq: hq)
    a = _add(aG, f(aR))
    b = dict(bG)
    for hq, c in bR.items():
        key = _add(aG, f(hq))
        b[key] = b.get(key, 0) + c
    return (a, (kG + kR) % 2, {k: v for k, v in b.items() if v})


def develop(word):
    G = ((0, 0), 0, {})
    for s in word:
        G = compose(G, REFL[s])
    return G


def is_translation(word):
    a, k, _ = develop(word)
    return k == 0 and a[0] == 0 and a[1] % 2 == 0


def offset(word, alpha, theta):
    """The return offset c(alpha) of a translation word, numerically (mpmath)."""
    a, k, b = develop(word)
    assert k == 0 and a[0] == 0 and a[1] % 2 == 0, 'not a translation word'
    L = cos(alpha) / sin(alpha)
    T = L * sum(c * exp(mpc(0, 1) * (2 * h * alpha + q * pi)) for (h, q), c in b.items())
    return -im(T * exp(mpc(0, -1) * theta)) / cos(theta)


def sine_sum(word_a, word_b, u, v):
    """offset(word_b) - offset(word_a) = -(L / cos theta) * sum_n C_n sin(n alpha), for
    theta = u alpha + v pi.  Returns {n: C_n} with n > 0, exact integers."""
    acc = {}
    for word, sg in ((word_b, 1), (word_a, -1)):
        for (h, q), c in develop(word)[2].items():
            # sin(2h alpha + q pi - u alpha - v pi) = (-1)^(q - v) sin((2h - u) alpha)
            n = 2 * h - u
            coef = sg * c * (-1) ** ((q - v) % 2)
            if n < 0:
                n, coef = -n, -coef
            if n:
                acc[n] = acc.get(n, 0) + coef
    return {n: c for n, c in sorted(acc.items()) if c}
