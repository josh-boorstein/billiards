#!/usr/bin/env python3
"""
exact_state.py — Exact triangle state using CosPoly arithmetic.

Each vertex coordinate is a CosPoly (polynomial in cos(2α) with integer
coefficients and power-of-2 denominator). Reflections are computed exactly.

The reflection formula is:
  P' = 2·proj(P, line AB) - P
     = 2·((P-A)·(B-A)/|B-A|²)·(B-A) + 2A - P

Since |B-A|² is constant (1 for L₁, cot²α for L₂, csc²α for H), and
cot²α = (1+cos2α)/(1-cos2α), csc²α = 2/(1-cos2α), dividing by these
amounts to multiplying by (1-cos2α)/2 = sin²α or its powers.

Strategy: track each coordinate as a CosPoly. When we need to divide by
cot²α or csc²α, we multiply numerator by sin²α (or appropriate power)
and adjust.
"""

from cos_poly import CosPoly


class ExactTriangleState:
    """Track vertex positions as exact polynomials in cos(2α).

    Each of Ox, Oy, Rx, Ry, Ax, Ay is a CosPoly.

    Initial values (for right triangle with angle α at O):
      O = (0, 0)
      R = (cot α, 0)
      A = (cot α, 1)

    cot α is NOT a polynomial in cos(2α) — it involves sin α in the denominator.
    So we work with SCALED coordinates: multiply all x-coordinates by sin α.
    Then:
      O_x = 0, O_y = 0
      R_x = cos α, R_y = 0    (since cot α · sin α = cos α)
      A_x = cos α, A_y = 1

    Wait — but the y-coordinates (which give the split points) are NOT scaled.
    The split point is the RAW y-coordinate of a vertex. If we scale x by sin α,
    the y-coordinates stay as-is, which is what we want.

    Actually, let's think more carefully. The reflection formula:
      P' = P + 2·((A-P)·n̂)·n̂  where n̂ ⊥ line AB, |n̂|=1

    Or equivalently:
      P' = 2·A + 2·((P-A)·d̂/|d|²)·d - P  where d = B-A

    The issue is that d = B-A involves mixed x and y coordinates, and
    dividing by |d|² introduces the denominator.

    KEY REALIZATION: We don't need to scale x-coordinates separately.
    Instead, let's track coordinates with a COMMON denominator that's
    a power of sin²α.

    After n reflections off H or L₂, all coordinates are of the form:
      (integer polynomial in cos(2α)) / (2^a · sin^{2b}(α))

    Since sin²α = (1-cos2α)/2, we have sin^{2b}α = ((1-cos2α)/2)^b.
    So the denominator is 2^a · ((1-cos2α)/2)^b = 2^{a-b} · (1-cos2α)^b.

    Rather than tracking rational functions, let's just keep everything
    as a CosPoly and absorb the sin²α factors into the numerator by
    clearing denominators at each step.

    SIMPLEST APPROACH: Just do the reflection formula directly in CosPoly
    arithmetic, clearing the |d|² denominator by multiplying through.
    The split-point comparison (is vertex_y > lo? < hi?) uses the
    evaluated float anyway — the exact coefficients let us evaluate at
    arbitrary precision if needed.

    Let me try the direct approach: store 6 CosPolys + a normalizing
    denominator (CosPoly). Split point = vertex_y_num / denom, evaluated.
    """

    def __init__(self, alpha=None):
        if alpha is None:
            return
        # Store alpha for float evaluation
        self.alpha = alpha
        # Initial coordinates:
        # O = (0, 0), R = (cot α, 0), A = (cot α, 1)
        # cot α = cos α / sin α. In the cos(2α) basis:
        #   cos α = ? Not directly in Z[cos(2α)].
        #   cos²α = (1 + cos(2α)) / 2
        #
        # We'll store everything multiplied by a common factor.
        # Let's use a different normalization:
        #   Scale all coordinates by sin(α). Then:
        #   O = (0, 0)  →  (0, 0)
        #   R = (cot α, 0)  →  (cos α, 0) ... but sin α · 0 = 0
        #
        # NO — scaling x by sin α changes only x, not y. The split point
        # is a y-coordinate. Let's NOT scale.
        #
        # Instead: track a single scalar denominator D (a CosPoly) such that
        # actual coordinate = stored_poly / D. Initially D = 1.
        # After each reflection, we may multiply D by |AB|² (a CosPoly)
        # and adjust numerators accordingly.

        # cot α = cos α / sin α. We need this as initial x-coordinate.
        # cos α is NOT in Z[cos(2α)] — it's in Z[cos(α)].
        # But cos²α = (1 + cos(2α))/2, so cot²α = (1+cos2α)/(1-cos2α).
        #
        # WORKAROUND: Store x-coords and y-coords with DIFFERENT denominators?
        # That's messy.
        #
        # BETTER WORKAROUND: The y-coordinates are what matter for split points.
        # The x-coordinates only matter for computing the NEXT y-coordinate
        # via the reflection formula. Let's track everything in a way that
        # the y-coordinates come out as pure CosPolys.
        #
        # OBSERVATION: In the original TriangleState, the ONLY thing we
        # extract is vertex_y (the split point). The x-coordinates are
        # intermediate. So we need y-coordinates to be exact; x-coordinates
        # just need to be correct.
        #
        # PRACTICAL APPROACH: Store all 6 coordinates as CosPolys with a
        # common denominator. Use sin²α = (1-cos2α)/2 to clear denominators.
        #
        # Initial state: O=(0,0), R=(cotα, 0), A=(cotα, 1).
        # cotα = cosα/sinα. This is NOT algebraic in cos(2α) alone.
        #
        # KEY INSIGHT FROM THE Y-COORDINATE STRUCTURE:
        # The y-coordinates of vertices in the unfolded frame ARE in Z[cos(2α)]
        # (we proved this — Theorem 8, integer lattice). But x-coordinates
        # involve cotα which is cos(α)/sin(α).
        #
        # Since we only need y-coordinates for split points, and y only
        # mixes with y through projections involving x... let's just
        # store x and y separately and handle the algebra.
        #
        # Actually, let's look at what happens in the reflection formula.
        # Reflect P=(px,py) across line through A=(ax,ay) and B=(bx,by):
        #   dx = bx-ax, dy = by-ay, d2 = dx²+dy²
        #   t = ((px-ax)*dx + (py-ay)*dy) / d2
        #   P'x = 2*(ax + t*dx) - px
        #   P'y = 2*(ay + t*dy) - py
        #
        # Expanding:
        #   t = ((px-ax)*dx + (py-ay)*dy) / d2
        #   P'y = 2*ay + 2*t*dy - py
        #       = 2*ay + 2*dy*((px-ax)*dx + (py-ay)*dy)/d2 - py
        #
        # So P'y depends on px through the cross-term (px-ax)*dx*dy/d2.
        # Unless dx*dy = 0 (which happens for L₁ and L₂ but NOT for H).
        #
        # For L₁ (vertical, A=(cotα,1) to R=(cotα,0)): dx=0, dy=-1, d2=1.
        #   P'y = 2*1 + 2*(-1)*((px-cotα)*0 + (py-1)*(-1))/1 - py
        #       = 2 + 2*(py-1) - py = py. Wait that's wrong...
        #   Let me redo: L₁ goes from R=(cotα,0) to A=(cotα,1).
        #   dx=0, dy=1, d2=1.
        #   t = (py - 0) * 1 / 1 = py  ... no, A is one endpoint.
        #   Actually reflect across line R-A: a=(cotα,0), b=(cotα,1).
        #   dx=0, dy=1, d2=1.
        #   t = ((px-cotα)*0 + (py-0)*1) / 1 = py
        #   P'x = 2*(cotα + py*0) - px = 2cotα - px
        #   P'y = 2*(0 + py*1) - py = py
        #   So reflecting across L₁: x flips, y unchanged.
        #   y-coordinate doesn't change! And x just flips around cotα.
        #
        # For L₂ (horizontal, O=(0,0) to R=(cotα,0)): dy=0.
        #   a=(0,0), b=(cotα,0). dx=cotα, dy=0, d2=cot²α.
        #   t = ((px-0)*cotα + (py-0)*0) / cot²α = px/cotα
        #   P'x = 2*(px/cotα * cotα) - px = px  (unchanged)
        #   P'y = 2*0 - py = -py
        #   So reflecting across L₂: y flips sign, x unchanged.
        #   y-coordinate just negates!
        #
        # For H (hypotenuse, O=(0,0) to A=(cotα,1)):
        #   a=(0,0), b=(cotα,1). dx=cotα, dy=1, d2=cot²α+1=csc²α.
        #   t = (px*cotα + py*1) / csc²α = (px*cotα + py) * sin²α
        #   P'x = 2*(t*cotα) - px = 2*cotα*(px*cotα+py)*sin²α - px
        #       = 2*cos²α*px + 2*cotα*sin²α*py - px
        #       = (2cos²α - 1)*px + 2*cosα*sinα*py
        #       = cos(2α)*px + sin(2α)*py
        #   P'y = 2*(t*1) - py = 2*(px*cotα+py)*sin²α - py
        #       = 2*cosα*sinα*px + 2*sin²α*py - py
        #       = sin(2α)*px + (2sin²α-1)*py
        #       = sin(2α)*px - cos(2α)*py
        #
        # BEAUTIFUL! The H reflection is:
        #   P'x = cos(2α)·px + sin(2α)·py
        #   P'y = sin(2α)·px - cos(2α)·py
        #
        # And L₂ is: P'x = px, P'y = -py
        # And L₁ is: P'x = 2cotα - px, P'y = py
        #
        # For the y-coordinates (which we care about):
        # - L₁ reflection: y unchanged
        # - L₂ reflection: y negates
        # - H reflection: y' = sin(2α)·px - cos(2α)·py
        #
        # The H reflection mixes x into y! So we DO need x-coordinates.
        # But note: sin(2α) = 2·sin(α)·cos(α). So the x-contribution
        # to y is sin(2α)·px. If px involves cotα = cosα/sinα, then
        # sin(2α)·cotα = 2cos²α = 1 + cos(2α). Clean!
        #
        # So let's track: each coordinate as a CosPoly times possibly
        # a factor of cotα (for x-coordinates only).
        #
        # REFINED REPRESENTATION:
        #   y-coordinates: pure CosPoly (in Z[1/2][cos(2α)])
        #   x-coordinates: CosPoly × cotα (i.e., store the coefficient
        #     of cotα, which is itself a CosPoly)
        #
        # Then the reflection rules become:
        #   For H: y' = sin(2α)·(x_coeff·cotα) - cos(2α)·y
        #            = 2sinα·cosα·cosα/sinα · x_coeff - cos(2α)·y
        #            = 2cos²α · x_coeff - cos(2α)·y
        #            = (1+cos(2α)) · x_coeff - cos(2α)·y
        #         x_coeff' = cos(2α)·x_coeff + sin(2α)·y / cotα
        #                  ... wait, sin(2α)/cotα = sin(2α)·sinα/cosα = 2sin²α
        #                  = (1-cos(2α))
        #         So: x_coeff' = cos(2α)·x_coeff + (1-cos(2α))·y
        #
        #   For L₁: y unchanged, x_coeff' = 2 - x_coeff
        #     (since x' = 2cotα - x = (2-x_coeff)·cotα)
        #
        #   For L₂: y' = -y, x_coeff unchanged
        #
        # ALL OPERATIONS ARE IN Z[1/2][cos(2α)]! No denominators beyond powers of 2!

        # Initialize: O=(0,0), R=(cotα,0), A=(cotα,1)
        # x stored as coefficient of cotα:
        #   Ox_c = 0, Rx_c = 1, Ax_c = 1
        # y stored directly:
        #   Oy = 0, Ry = 0, Ay = 1

        self.Ox_c = CosPoly([0], 0)  # O.x / cotα = 0
        self.Oy = CosPoly([0], 0)    # O.y = 0
        self.Rx_c = CosPoly([1], 0)  # R.x / cotα = 1
        self.Ry = CosPoly([0], 0)    # R.y = 0
        self.Ax_c = CosPoly([1], 0)  # A.x / cotα = 1
        self.Ay = CosPoly([1], 0)    # A.y = 1

    def copy(self):
        s = ExactTriangleState.__new__(ExactTriangleState)
        s.alpha = self.alpha
        s.Ox_c = self.Ox_c.copy()
        s.Oy = self.Oy.copy()
        s.Rx_c = self.Rx_c.copy()
        s.Ry = self.Ry.copy()
        s.Ax_c = self.Ax_c.copy()
        s.Ay = self.Ay.copy()
        return s

    def reflect(self, side):
        """Reflect the appropriate vertex across the given side.

        H reflects R: R' = reflect(R, line O-A)
        L₂ reflects A: A' = reflect(A, line O-R)
        L₁ reflects O: O' = reflect(O, line R-A)

        Using the derived formulas:
          H: y' = (1+cos2α)·x_c - cos(2α)·y
             x_c' = cos(2α)·x_c + (1-cos(2α))·y
          L₁: y' = y (unchanged)
              x_c' = 2·(line_x_c) - x_c  ... but line is R-A
              Actually for L₁ (reflect across R-A which is vertical at x=cotα):
              x' = 2cotα - x → x_c' = 2 - x_c...
              WAIT: this only works if the line is at x=cotα.
              But after reflections, R and A have moved! The line R-A
              is no longer vertical.

              This is the fundamental issue: the "sides" move.
              The reflection is across the CURRENT side, not the original.
        """
        # The CURRENT sides are defined by the current vertex positions.
        # L₁ is always R-A (current R and A positions).
        # L₂ is always O-R (current O and R positions).
        # H is always O-A (current O and A positions).
        #
        # The clean decomposition (H: x'=cos2α·x+sin2α·y, etc.) only
        # works for the ORIGINAL triangle. After the first reflection,
        # the sides have moved.
        #
        # However! There's a key structure: reflections are isometries.
        # The triangle is always congruent to the original. The sides
        # have the same lengths. But the DIRECTION of H changes.
        #
        # We need to do the general reflection formula:
        #   Reflect P across line through A, B:
        #   d = B - A, d2 = |d|², t = (P-A)·d / d2
        #   P' = 2A + 2t·d - P ... wait no:
        #   P' = 2(A + t·d) - P
        #
        # With the (x_c, y) representation where x = x_c · cotα:
        #   Actual P = (Px_c · cotα, Py)
        #   Actual A = (Ax_c · cotα, Ay)
        #   Actual B = (Bx_c · cotα, By)
        #
        #   d = ((Bx_c - Ax_c)·cotα, By - Ay)
        #   d2 = (Bx_c - Ax_c)²·cot²α + (By - Ay)²
        #
        # d2 involves cot²α = (1+cos2α)/(1-cos2α)... which is a rational
        # function, not a polynomial. This brings us back to the denominator
        # problem.
        #
        # HOWEVER: we know d2 is a CONSTANT (doesn't depend on which
        # reflection step we're at, only on WHICH side). Specifically:
        #   |RA|² = 1, |OR|² = cot²α, |OA|² = csc²α = 1+cot²α
        #
        # So d2 is one of: 1, cot²α, or csc²α.
        # We need t = numerator / d2, and then P' = 2(A+t·d) - P.
        # Expanding P'y = 2·Ay + 2·t·dy - Py.
        #
        # To avoid dividing by cot²α or csc²α, multiply EVERYTHING
        # by d2. Track a common denominator that accumulates the d2 factors.
        #
        # But d2 = cot²α for L₂, and cot²α is not a polynomial...
        #
        # ALTERNATIVE: use (sin α · x, y) as coordinates instead of (x, y).
        # Let u = sin(α)·x, v = y. Then:
        #   u_O = 0, v_O = 0
        #   u_R = sin(α)·cot(α) = cos(α), v_R = 0
        #   u_A = cos(α), v_A = 1
        #
        # cos(α) is NOT in Z[cos(2α)] either (it's in Z[cos(α)]).
        # But cos²α = (1+cos2α)/2 IS.
        #
        # Let me try yet another scaling: use (sin²α · x, sin α · y)?
        # No, that makes y messy.
        #
        # FINAL APPROACH: Use the reflection formula with the CONSTANT d2,
        # and represent d2 as a CosPoly (rational in cos2α), keeping a
        # separate denominator polynomial.
        #
        # Actually, the cleanest thing: just use the formulas from the
        # FIRST reflection (which IS at the original sides) and note that
        # the word_tree's TriangleState does exactly one reflection per
        # `reflect` call, always across the CURRENT side.
        #
        # Let me just implement the GENERAL reflection in CosPoly
        # arithmetic with the (x_c · cotα, y) representation, and
        # handle the d2 denominators by clearing them.

        if side == 'H':
            # Reflect R across line O-A (current positions)
            # Line from (Ox_c·cotα, Oy) to (Ax_c·cotα, Ay)
            # dx = (Ax_c - Ox_c)·cotα, dy = Ay - Oy
            # d2 = (Ax_c-Ox_c)²·cot²α + (Ay-Oy)²
            # But |OA| is always csc²α? NO — only in the original.
            # After reflections, the sides can have different orientations
            # but SAME LENGTH. So |OA|² = csc²α always.
            #
            # Hmm wait — actually no. When we reflect R across O-A,
            # O and A don't move. So O-A stays the same as before.
            # When we reflect A across O-R, O and R don't move.
            # When we reflect O across R-A, R and A don't move.
            #
            # So the LENGTH of the reflection line is always one of the
            # three side lengths of the ORIGINAL triangle! The fixed pair
            # of vertices hasn't moved (or moved together preserving length).
            #
            # Wait, that's not right either. After the first reflection
            # (H, reflecting R), R moves but O and A don't. So |OR| and |RA|
            # change, but |OA| stays the same.
            # After the second reflection (say L₂, reflecting A across O-R),
            # O and R don't move, so |OR| stays the same.
            # But R has ALREADY moved from step 1! So |OR| is the distance
            # from O to the NEW R.
            #
            # Actually NO. The point is that each vertex that moves gets
            # reflected across the OPPOSITE side. The two endpoints of
            # that side DON'T move in that step. But they may have moved
            # in PREVIOUS steps.
            #
            # The key invariant: the triangle is always congruent to the
            # original (reflections preserve distances). So |OA| = csc α,
            # |OR| = cot α, |RA| = 1 AT ALL TIMES. The side lengths
            # never change!
            #
            # This means d2 is ALWAYS:
            #   reflecting across O-A (H): d2 = |OA|² = csc²α
            #   reflecting across O-R (L₂): d2 = |OR|² = cot²α
            #   reflecting across R-A (L₁): d2 = |RA|² = 1

            # Great! So for H reflection (reflect R across O-A):
            # d2 = csc²α = 1/sin²α
            #
            # General formula: P'y = Py + 2·dy·((Px-Ax)·dx+(Py-Ay)·dy)/d2 ...
            # No wait: P' = 2·proj - P where proj is projection of P onto line.
            #
            # Let me use: P' = 2A + 2t(B-A) - P where t = (P-A)·(B-A)/|B-A|²
            # No that's wrong too. Standard reflection of P across line AB:
            #   foot = A + ((P-A)·(B-A)/|B-A|²)(B-A)
            #   P' = 2·foot - P
            #
            # P'y = 2·foot_y - Py
            # foot_y = Ay + t·(By-Ay) where t = ((Px-Ax)(Bx-Ax)+(Py-Ay)(By-Ay))/|BA|²
            # Here A=O, B=A (the vertex), reflecting R.
            # Wait, reflecting R across side H = line from O to A(vertex).
            # Let's call the line endpoints (to avoid confusion) as P1=O, P2=A(vertex).

            # Reflect R across line O-to-A:
            # In (x_c·cotα, y) coordinates:
            # R = (Rx_c·cotα, Ry), O = (Ox_c·cotα, Oy), A = (Ax_c·cotα, Ay)
            # d = A - O = ((Ax_c-Ox_c)·cotα, Ay-Oy)
            # d2 = |d|² = (Ax_c-Ox_c)²·cot²α + (Ay-Oy)² = csc²α (constant)
            #
            # t = ((R-O)·d) / d2
            #   = ((Rx_c-Ox_c)·cotα · (Ax_c-Ox_c)·cotα + (Ry-Oy)·(Ay-Oy)) / csc²α
            #   = ((Rx_c-Ox_c)(Ax_c-Ox_c)·cot²α + (Ry-Oy)(Ay-Oy)) · sin²α
            #
            # R'y = Oy + t·(Ay-Oy) + (Oy + t·(Ay-Oy) - Ry)
            #     = 2·(Oy + t·(Ay-Oy)) - Ry
            #
            # This is getting complex but tractable. Let me implement it step by step.

            self._reflect_across(
                # point to reflect: R
                'R',
                # line endpoints: O, A
                'O', 'A'
            )

        elif side == 'L2':
            # Reflect A across line O-R
            self._reflect_across('A', 'O', 'R')

        elif side == 'L1':
            # Reflect O across line R-A
            self._reflect_across('O', 'R', 'A')

    def _reflect_across(self, point_label, line_start, line_end):
        """Reflect point across line defined by two other vertices.

        Uses the representation (x_c · cotα, y) for all points.

        The formula: P' = 2·foot - P, where foot is the projection of P
        onto the line through A,B.

        foot = A + t·(B-A), t = (P-A)·(B-A) / |B-A|²

        Since |B-A|² is constant:
          |RA|² = 1
          |OR|² = cot²α
          |OA|² = csc²α = 1 + cot²α

        We multiply through by |B-A|² (expressed in our basis) to avoid
        division. Then divide at the end by the same factor.

        In the (x_c·cotα, y) representation:
          dot product (P-A)·(B-A) = (Px_c-Ax_c)(Bx_c-Ax_c)·cot²α + (Py-Ay)(By-Ay)

        And |B-A|² = (Bx_c-Ax_c)²·cot²α + (By-Ay)²

        Since cot²α = (1+cos2α)/(1-cos2α), multiplying by (1-cos2α) = 2sin²α
        clears the denominator. Let's define:
          D2_cleared = |B-A|² · (1-cos2α) = (Bx_c-Ax_c)²·(1+cos2α) + (By-Ay)²·(1-cos2α)

        And the dot product cleared:
          dot_cleared = (Px_c-Ax_c)(Bx_c-Ax_c)·(1+cos2α) + (Py-Ay)(By-Ay)·(1-cos2α)

        Then t = dot_cleared / D2_cleared.

        foot_y = Ay + t·(By-Ay) = (Ay·D2_cleared + dot_cleared·(By-Ay)) / D2_cleared
        P'y = 2·foot_y - Py = (2Ay·D2_cleared + 2·dot_cleared·(By-Ay) - Py·D2_cleared) / D2_cleared

        Similarly for x_c:
        foot_xc = Ax_c + t·(Bx_c-Ax_c) = (Ax_c·D2_cleared + dot_cleared·(Bx_c-Ax_c)) / D2_cleared
        P'x_c = 2·foot_xc - Px_c = (2Ax_c·D2 + 2·dot·(Bx_c-Ax_c) - Px_c·D2) / D2_cleared

        Since D2_cleared is a CONSTANT (only depends on which side, not on
        the coordinates), and the triangle is always congruent to the original,
        D2_cleared is the same every time we reflect across the same side type.

        For |RA|²=1: D2_cleared = 0²·(1+c2α) + 1²·(1-c2α) = 1-cos2α = 2sin²α
        For |OR|²=cot²α: D2_cleared = 1²·(1+c2α) + 0²·(1-c2α) = 1+cos2α = 2cos²α
        For |OA|²=csc²α: D2_cleared = 1²·(1+c2α) + 1²·(1-c2α) = 2

        Wait but that assumes Bx_c-Ax_c and By-Ay are the ORIGINAL differences.
        After reflections the vertices move! The side R-A in the original triangle
        has dx_c=0, dy=1, but after reflecting R, it changes.

        Hmm, this is the crux issue. Let me think differently.

        The triangle is CONGRUENT to the original, so |R-A|²=1, |O-R|²=cot²α,
        |O-A|²=csc²α always. But the DIRECTION of the sides changes.

        For the general formula, we need to compute:
          d2 = (Bx_c-Ax_c)²·cot²α + (By-Ay)² = CONSTANT (one of 1, cot²α, csc²α)

        So: (Bx_c-Ax_c)²·cot²α + (By-Ay)² = known constant.
        We can use this to simplify, but the intermediate terms still depend
        on the actual coordinates.

        OK, I think the cleanest implementation is:
        1. Compute the dot product and d2 in terms of our CosPoly coordinates
        2. Use the known d2 constant to simplify
        3. Express the result as CosPoly / scalar

        Let me just implement it numerically with CosPolys, clearing the
        cot²α factors by multiplying by sin²α as needed.
        """
        # Get coordinates
        Px_c, Py = self._get_vertex(point_label)
        Ax_c, Ay = self._get_vertex(line_start)
        Bx_c, By = self._get_vertex(line_end)

        # Differences
        dAx = Bx_c - Ax_c  # (B-A) x-coefficient of cotα
        dAy = By - Ay       # (B-A) y-component
        dPx = Px_c - Ax_c   # (P-A) x-coefficient of cotα
        dPy = Py - Ay       # (P-A) y-component

        # dot product (P-A)·(B-A) = dPx·dAx·cot²α + dPy·dAy
        # d2 = |B-A|² = dAx²·cot²α + dAy²
        #
        # Since cot²α = (1+cos2α)/(1-cos2α), let's multiply everything
        # by (1-cos2α) = 2sin²α to clear it.
        #
        # dot_cleared = dPx·dAx·(1+cos2α) + dPy·dAy·(1-cos2α)
        # d2_cleared = dAx²·(1+cos2α) + dAy²·(1-cos2α)
        #
        # Note: (1+cos2α) and (1-cos2α) are CosPolys:
        one_plus_c2 = CosPoly([1, 1], 0)   # 1 + cos(2α)
        one_minus_c2 = CosPoly([1, -1], 0)  # 1 - cos(2α)

        # dot_cleared = dPx·dAx·(1+c2) + dPy·dAy·(1-c2)
        dot_cleared = dPx * dAx * one_plus_c2 + dPy * dAy * one_minus_c2

        # d2_cleared = dAx²·(1+c2) + dAy²·(1-c2)
        d2_cleared = dAx * dAx * one_plus_c2 + dAy * dAy * one_minus_c2

        # P' = 2·(A + t·(B-A)) - P where t = dot_cleared / d2_cleared
        # P'y = (2·Ay·d2_cleared + 2·dot_cleared·dAy - Py·d2_cleared) / d2_cleared
        # P'x_c = (2·Ax_c·d2_cleared + 2·dot_cleared·dAx - Px_c·d2_cleared) / d2_cleared

        # Numerators:
        new_Py_num = 2 * Ay * d2_cleared + 2 * dot_cleared * dAy - Py * d2_cleared
        new_Px_c_num = 2 * Ax_c * d2_cleared + 2 * dot_cleared * dAx - Px_c * d2_cleared

        # Divide by d2_cleared. Since d2_cleared is a known constant (polynomial
        # that evaluates to 2sin²α, 2cos²α, or 2 depending on side), we can
        # divide exactly.
        #
        # But dividing one CosPoly by another is NOT trivial — it's polynomial
        # division that may not be exact in general.
        #
        # HOWEVER: since the triangle is always congruent to the original,
        # the result MUST be exact (the new coordinates are also in Z[1/2][cos2α]).
        # So d2_cleared divides the numerator exactly.
        #
        # For now: evaluate d2_cleared at our specific alpha to get a scalar,
        # and divide the numerator coefficients by it. If it's a simple
        # constant (like 2), this is trivial.
        #
        # Actually, d2_cleared is NOT constant in general! It depends on
        # dAx and dAy which are themselves CosPolys (current vertex positions).
        #
        # The LENGTH |B-A|² is constant, but d2_cleared = |B-A|²·(1-cos2α)
        # = |B-A|²·2sin²α is also constant (a fixed number for each side type).
        # But dAx²·(1+c2) + dAy²·(1-c2) = 2sin²α·|B-A|² should always equal
        # the same value... is that true?
        #
        # |B-A|² = dAx²·cot²α + dAy² (in original coords)
        # d2_cleared = dAx²·(1+c2) + dAy²·(1-c2)
        #            = dAx²·2cos²α + dAy²·2sin²α  (since 1+c2=2cos²α, 1-c2=2sin²α)
        #            ... no, 1+cos2α = 2cos²α and 1-cos2α = 2sin²α. Yes!
        # So d2_cleared = 2·(dAx²·cos²α + dAy²·sin²α)
        # But |B-A|² = dAx²·cot²α + dAy² = dAx²·cos²α/sin²α + dAy²
        #            = (dAx²·cos²α + dAy²·sin²α) / sin²α
        # So d2_cleared = 2·|B-A|²·sin²α
        #
        # Since |B-A|² is constant, d2_cleared = 2·C·sin²α where C is:
        #   H side (O-A): C = csc²α → d2_cleared = 2·csc²α·sin²α = 2
        #   L₂ side (O-R): C = cot²α → d2_cleared = 2·cot²α·sin²α = 2cos²α = 1+cos2α
        #   L₁ side (R-A): C = 1 → d2_cleared = 2sin²α = 1-cos2α
        #
        # These ARE simple CosPolys! And they're the same regardless of
        # how many reflections we've done. So we CAN divide exactly by
        # multiplying by the inverse (which for these simple polynomials
        # requires polynomial division).
        #
        # But general polynomial division in Z[cos2α] is not straightforward.
        # Division by (1+cos2α) or (1-cos2α) requires the numerator to be
        # divisible by these — which it must be by the congruence argument.
        #
        # Let me implement exact division by these simple polynomials.

        # Determine which side we're reflecting across
        side_pair = {line_start, line_end}
        if side_pair == {'O', 'A'}:
            # H side: d2_cleared = 2 (constant)
            # Just divide all coefficients by 2 (add 1 to shift)
            new_Py_num.shift += 1
            new_Py_num._reduce()
            new_Px_c_num.shift += 1
            new_Px_c_num._reduce()
            self._set_vertex(point_label, new_Px_c_num, new_Py_num)
        elif side_pair == {'O', 'R'}:
            # L₂ side: d2_cleared = 1 + cos(2α)
            new_Py = _divide_by_1_plus_cos2(new_Py_num)
            new_Px_c = _divide_by_1_plus_cos2(new_Px_c_num)
            self._set_vertex(point_label, new_Px_c, new_Py)
        elif side_pair == {'R', 'A'}:
            # L₁ side: d2_cleared = 1 - cos(2α)
            new_Py = _divide_by_1_minus_cos2(new_Py_num)
            new_Px_c = _divide_by_1_minus_cos2(new_Px_c_num)
            self._set_vertex(point_label, new_Px_c, new_Py)

    def _get_vertex(self, label):
        if label == 'O':
            return self.Ox_c, self.Oy
        elif label == 'R':
            return self.Rx_c, self.Ry
        elif label == 'A':
            return self.Ax_c, self.Ay

    def _set_vertex(self, label, x_c, y):
        if label == 'O':
            self.Ox_c, self.Oy = x_c, y
        elif label == 'R':
            self.Rx_c, self.Ry = x_c, y
        elif label == 'A':
            self.Ax_c, self.Ay = x_c, y

    def vertex_y(self, vtype):
        """Get the y-coordinate CosPoly for a vertex."""
        if vtype == 'O':
            return self.Oy
        elif vtype == 'R':
            return self.Ry
        elif vtype == 'A':
            return self.Ay

    def vertex_y_float(self, vtype):
        """Evaluate vertex y-coordinate at the stored alpha."""
        return self.vertex_y(vtype).evaluate(self.alpha)


def _divide_by_1_plus_cos2(poly):
    """Exact division by (1 + cos(2α)) using Chebyshev recurrence.

    If f = (1+T₁)·g where T₁ = cos(2α), and f = Σ fₖ Tₖ, g = Σ gₖ Tₖ:
      f_n = g_{n-1}/2  (for n≥2)
      f_k = g_k + (g_{k-1} + g_{k+1})/2  for k≥2
      f_1 = g_1 + g_0 + g_2/2
      f_0 = g_0 + g_1/2

    Solving top-down: g_{n-1} = 2f_n, then recur downward.
    All arithmetic stays in integers (with shift tracking).
    """
    return _divide_by_1_pm_cos2(poly, +1)


def _divide_by_1_minus_cos2(poly):
    """Exact division by (1 - cos(2α)) using Chebyshev recurrence."""
    return _divide_by_1_pm_cos2(poly, -1)


def _divide_by_1_pm_cos2(poly, sign):
    """Divide poly by (1 + sign·cos(2α)) exactly in Chebyshev basis.

    sign=+1: divide by (1+cos2α). sign=-1: divide by (1-cos2α).

    Uses the recurrence from (1+s·T₁)·g = f:
      f_n = s·g_{n-1}/2
      f_k = g_k + s·(g_{k-1}+g_{k+1})/2  for k≥2
      f_1 = g_1 + s·g_0 + s·g_2/2
      f_0 = g_0 + s·g_1/2

    Solved top-down: g_{n-1}=2s·f_n, g_{k-1}=2s(f_k-g_k)-g_{k+1}, then g_0.
    """
    if poly.is_zero():
        return CosPoly([0], 0)

    f = poly.coeffs[:]
    n = len(f) - 1
    while n > 0 and f[n] == 0:
        n -= 1

    if n == 0:
        if f[0] == 0:
            return CosPoly([0], 0)
        assert False, f"Cannot divide constant {f[0]} by (1{'+' if sign>0 else '-'}cos2α)"

    # g has degree n-1. All g values are integers at shift = poly.shift.
    # The ×2 in g_{n-1}=2s·f_n offsets the /2 implicit in the recurrence.
    g = [0] * n
    result_shift = poly.shift

    if n == 1:
        # f_1 = s·g_0 → g_0 = s·f_1 (from the k=1 eq with g_1=0, g_2=0...
        # Actually f_1 = s·g_0 (since g_1=0), so g_0 = f_1·s = f_1/s = s·f_1
        g[0] = sign * f[1]
    else:
        # Top-down: g_{n-1} = 2·sign·f_n
        g[n-1] = 2 * sign * f[n]

        # For k from n-1 down to 2: g_{k-1} = 2s(f_k - g_k) - g_{k+1}
        for k in range(n-1, 1, -1):
            g_kp1 = g[k+1] if k + 1 <= n - 1 else 0
            g[k-1] = 2 * sign * (f[k] - g[k]) - g_kp1

        # g_0: from f_1 = g_1 + s·g_0 + s·g_2/2
        # → g_0 = s·(f_1 - g_1) - g_2/2
        g2 = g[2] if n > 2 else 0
        if g2 % 2 == 0:
            g[0] = sign * (f[1] - g[1]) - g2 // 2
        else:
            g = [c * 2 for c in g]
            result_shift += 1
            g[0] = 2 * sign * (f[1] - g[1] // 2) - (g2 * 2) // 2
            # Wait: after doubling, g[1] is the doubled value. g_1_orig = g[1]//2.
            # Formula: 2·g_0 = 2·sign·(f_1 - g_1_orig) - g_2_orig
            #                 = 2·sign·f_1 - 2·sign·(g[1]//2) - g2
            # Simpler: just scale f too
            g[0] = sign * (2 * f[1] - g[1]) - g2

    result = CosPoly(g, result_shift)
    result._reduce()
    return result


def test_exact_state():
    """Verify ExactTriangleState matches TriangleState floats."""
    import math
    from word_tree import TriangleState

    for p, q in [(1, 5), (2, 7), (1, 7), (3, 7), (1, 11), (2, 11)]:
        alpha = p * math.pi / q
        if alpha >= math.pi / 2:
            continue

        # Test a sequence of reflections
        words = [
            ['H'],
            ['H', 'L2'],
            ['H', 'L1'],
            ['H', 'L2', 'H'],
            ['H', 'L2', 'L1'],
            ['H', 'L1', 'H'],
            ['H', 'L1', 'L2'],
            ['H', 'L2', 'H', 'L2'],
            ['H', 'L2', 'H', 'L1'],
            ['H', 'L1', 'H', 'L2', 'H'],
            ['H', 'L1', 'H', 'L2', 'L1', 'H', 'L2'],
        ]

        errors = []
        for word in words:
            float_state = TriangleState(alpha)
            exact_state = ExactTriangleState(alpha)

            for side in word:
                float_state.reflect(side)
                exact_state.reflect(side)

            # Compare y-coordinates
            for v in ['O', 'R', 'A']:
                float_y = float_state.vertex_y(v)
                exact_y = exact_state.vertex_y_float(v)
                err = abs(float_y - exact_y)
                if err > 1e-10:
                    errors.append((p, q, ''.join(s[0] for s in word), v, float_y, exact_y, err))

        if errors:
            print(f"p={p} q={q}: ERRORS")
            for e in errors[:5]:
                print(f"  word={e[2]} vertex={e[3]}: float={e[4]:.10f} exact={e[5]:.10f} err={e[6]:.2e}")
        else:
            print(f"p={p} q={q}: all {len(words)} words match (max word len {max(len(w) for w in words)})")


if __name__ == "__main__":
    test_exact_state()
